

random_sampling <- function(aln, max_reads, seed = 42){
    set.seed(seed)
    aln <- aln[sample(seq_along(aln), min(max_reads, length(aln)))]
    
    return(aln)
}

by_start_read <- function(aln, max_reads, strand){
    
    if (strand == "+") {
        aln <- aln[order(start(aln), decreasing = FALSE)][1:max_reads]
    } else {
        aln <- aln[order(end(aln), decreasing = TRUE)][1:max_reads]
    }
    
    return(aln)
    
}

by_read_length <- function(aln, max_reads, strand){
    
    aln <- aln[order(width(aln), decreasing = TRUE)][1:max_reads]
    
    return(aln)
    
}

by_start_te_island <- function(aln, max_reads, strand, te){
    
    message("Filtering reads by start position relative to TE island...")
    message(paste0("Total n reads before filtering: ", length(aln)))
    read_starts <- granges(aln)
    
    if (as.character(strand(te)) == "+"){
        
        start_coord = start(te) - 1000
        end_coord = end(te) + 1000
        window <- GenomicRanges::GRanges(seqnames = as.character(seqnames(te)),
                                         ranges = IRanges::IRanges(start = start_coord,
                                                                  end = end_coord),
                                         strand = "+")
        
        start(read_starts) <- start(aln)
        end(read_starts) <- start(aln) + 1
        
    } else {
        
        start_coord = start(te) - 1000
        end_coord = end(te) + 1000
        
        window <- GenomicRanges::GRanges(seqnames = as.character(seqnames(te)),
                                         ranges = IRanges::IRanges(start = start_coord,
                                                                  end = end_coord))
        
        start(read_starts) <- end(aln) - 1
        end(read_starts) <- end(aln) 
        
    }
    
   
   message(paste0("Filtering reads starting within: ", as.character(window)))
   message(paste0("Window start: ", start_coord, " end: ", end_coord))
   message(paste0("Total n reads before filtering: ", length(read_starts)))
    
   # browser()
    
   hits <- findOverlaps(read_starts, window)
   
   message(paste0("n reads overlapping the window: ", length(hits)))
   
   if (length(hits) < max_reads){
       max_reads = length(hits)
   }
   
   if (max_reads > 0){
       filtered_reads <- aln[queryHits(hits)][1:max_reads]
   }else{
       filtered_reads <- NULL
   }
   
   message(paste0("n reads after filtering: ", length(filtered_reads)))
   
   return(filtered_reads)
    
}

sub_sampling_reads <- function(bam_file, chromosome, start_coord, end_coord, strand, max_reads = 30, sampling_strat = 'by_pos'){
    
    message(paste0("Reading alignments from BAM file: ", bam_file))
    which <- GenomicRanges::GRanges(seqnames = chromosome,
                                    ranges = IRanges::IRanges(start = start_coord, end = end_coord),
                                    strand = strand)
    
    param <- Rsamtools::ScanBamParam(which = which)
    
    aln <- GenomicAlignments::readGAlignments(bam_file, param = param)
    
    message(paste0("n alignments: ", length(aln)))
    # Sub sample to avoid plotting too many reads
    if (length(aln) > max_reads) {
        if (sampling_strat == 'by_pos') {
            print("Sampling by position...")
            aln <- by_start_read(aln, max_reads, strand = strand)
        }
        else if (sampling_strat == 'by_length') {
            print("Sampling by read length...")
            aln <- by_read_length(aln, max_reads, strand = strand)
        }
        #} else if(sampling_strat == 'by_te'){
        #    print("Sampling by TE...")
        #    aln <- by_start_te_island(aln, max_reads, strand, which)
        
        else {
            print("Random sampling...")
            aln <- random_sampling(aln, max_reads)
        }
    }
    message("Filter for reads beginning nearby the TE island...")
    aln <- by_start_te_island(aln, max_reads, strand, which)
    if (!is.null(aln)) {
        aln <- granges(aln)
    }
    
    return(aln)
}


te_island_track <- function(te_island_id_list,
                            te_island_annotation_file){
    
    te_island_ranges <- genomation::readGeneric(te_island_annotation_file, 
                                              strand = 6, 
                                              meta.cols = list(names = 4))
    
    id_names_list <- str_replace(te_island_id_list, "_Cluster_", " island ")
    
    AnnotationTrack(te_island_ranges[te_island_ranges@elementMetadata$names %in% te_island_id_list],
                    name = 'TE island',
                    id = id_names_list,
                    featureAnnotation = "id",
                    fontcolor.item = "white",
                    fill = "black",
                    size = 0.5)
}

coloring_te_instances <- function(te_instance_class_list){
    # map the order.color to the te_instance_class_list to get a list of colors that should be returned by this
    # function
    colors <- order.color[te_instance_class_list]
    colors[is.na(colors)] <- "black"
    return(colors)
    
}

local_te_track <- function(te_ranges, te){
    
    # Subset to TEs in the region of interest
    local_tes = data.frame(te_ranges) %>% 
        filter(seqnames == as.character(seqnames(te)), 
               strand == as.character(strand(te)),
               dplyr::between(start, start(te), end(te))) %>% 
        blackRcloud::splitTEID("te.id") %>% 
        unique()
    
    color <- coloring_te_instances(local_tes$order)
        
    ranges_of_interest <- te_ranges[te_ranges@elementMetadata$te.id %in% local_tes$te.id & !duplicated(te_ranges@elementMetadata$te.id)]
    ids_of_ranges_of_interest <- ranges_of_interest %>% 
            as.data.frame() %>% 
            blackRcloud::splitTEID("te.id") %>%
            dplyr::pull(super_family)
    
    te_instance_track <- AnnotationTrack(ranges_of_interest,
                                  featureAnnotation = "id",
                                  name = 'TE instances',
                                  id = ids_of_ranges_of_interest,
                                  stacking = "dense",
                                  #col = 'white',
                                  fill = color,
                                  size = 0.5)
    
    return(te_instance_track)
}


nanopore_track <- function(te, bam_file, max_reads = 25, name = "nanopore reads", sampling_strat = 'by_te'){

    
    readRanges <- sub_sampling_reads(bam_file, 
                                     as.character(seqnames(te)), 
                                     start(te) - 1500, 
                                     end(te) + 5000,
                                     as.character(strand(te)),
                                     max_reads = max_reads,
                                     sampling_strat = sampling_strat)
    
    if (is.null(readRanges)) {
        return(NULL)
    }
    
    message(paste0("n reads to plot: ", length(readRanges)))
    
    nanopore_reads <- Gviz::AnnotationTrack(
        range = readRanges,
        genome = "mm10",
        name = "nanopore reads",
        fill = "lightgray",
        chromosome = as.character(seqnames(te)),
        id = c(),
        size = 0.5
    )
    
}


cage_seq_tracks <- function(te){
    if (as.character(strand(te)) == "+"){
        
        bw <- rtracklayer::import.bw(cage_fwd_bw) %>%
            plyranges::filter_by_overlaps(as(paste0(
                as.character(seqnames(te)), ':', start(te), '-', end(te)
            ), "GRanges"))
        
        track_col = "darkgreen"
    }else {
        
        bw <- rtracklayer::import.bw(cage_rev_bw) %>%
            plyranges::filter_by_overlaps(as(paste0(
                as.character(seqnames(te)), ':', start(te), '-', end(te)
            ), "GRanges"))
        
        track_col = "lightgreen"
    }
    
    
    cage_limits = c(0, plyr::round_any(max(bw$score), 10, f = ceiling))
    
    cage_track <-  DataTrack(
        range = bw,
        type = "h",
        name = 'CAGE',
        col = track_col,
        ylim = cage_limits
    )
    
    return(cage_track)
    
}

genomeBrowserNanoporeTracks <- function(te_island_id, 
                                  te, 
                                  te_island_annotation_file,
                                  bam_file){
    
    gtrack <- GenomeAxisTrack()
    te_island_track <- te_island_track(te_island_id_list = te_island_id, te_island_annotation_file)
    te_instance_track <- local_te_track(teRanges, te)
    nanopore_track <- nanopore_track(te, bam_file) 
    cage_tracks <- cage_seq_tracks(te)
    
    if (is.null(nanopore_track)) {
        return(NULL)
    }
    
    tracks <- list(
        #ideoTrack,
        gtrack,
        te_island_track,
        te_instance_track,
        cage_tracks,
        nanopore_track
    )
    
    return(tracks)
    
    
}

plotGenomeBrowserNanopore <- function(track_list, te){
    grid.grabExpr(Gviz::plotTracks(
        unname(track_list),
        #add35 = TRUE,
        shape = 'box',
        scale = 0.25,
        #chromosome = seqnames(te),
        from = start(te) - 1500,
        to = end(te) + 1500,
        col.axis = "black",
        col.title = "black",
        cex.feature = 0.7,
        background.title = "white",
        fontcolor.legend = "black"
    ))
}
