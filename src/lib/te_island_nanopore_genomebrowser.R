

random_sampling <- function(aln, max_reads, seed = 42){
    set.seed(seed)
    aln <- aln[sample(seq_along(aln), min(max_reads, length(aln)))]
    
    return(aln)
}

by_position <- function(aln, max_reads, strand){
    
    if (strand == "+") {
        aln <- aln[order(start(aln), decreasing = FALSE)][1:max_reads]
    } else {
        aln <- aln[order(end(aln), decreasing = TRUE)][1:max_reads]
    }
    
    return(aln)
    
}

sub_sampling_reads <- function(bam_file, chromosome, start_coord, end_coord, strand, max_reads = 25, sampling_strat = 'by_pos'){
    
    which <- GenomicRanges::GRanges(seqnames = chromosome,
                                    ranges = IRanges::IRanges(start = start_coord, end = end_coord))
    
    param <- Rsamtools::ScanBamParam(which = which)
    
    aln <- GenomicAlignments::readGAlignments(bam_file, param = param)
    # Sub sample to avoid plotting too many reads
    if (length(aln) > max_reads) {
        if (sampling_strat == 'by_pos') {
            print("Sampling by position...")
            aln <- by_position(aln, max_reads, strand = strand)
        } else {
            print("Random sampling...")
            aln <- random_sampling(aln, max_reads)
        }
    }
    
    return(granges(aln))
}


te_island_track <- function(te_island_id_list,
                            te_island_ranges,
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
                    size = 2)
}


local_te_track <- function(te_ranges, te){
    
    # Subset to TEs in the region of interest
    local_tes = data.frame(te_ranges) %>% 
        filter(seqnames == as.character(seqnames(te)), 
               strand == as.character(strand(te)),
               dplyr::between(start, start(te), end(te))) %>% 
        blackRcloud::splitTEID("te.id") %>% 
        unique()
        
    ranges_of_interest <- te_ranges[te_ranges@elementMetadata$te.id %in% local_tes$te.id & !duplicated(te_ranges@elementMetadata$te.id)]
    ids_of_ranges_of_interest <- ranges_of_interest %>% 
            as.data.frame() %>% 
            blackRcloud::splitTEID("te.id") %>%
            dplyr::pull(super_family)
    
    te_instance_track <- AnnotationTrack(ranges_of_interest,
                                  featureAnnotation = "id",
                                  id = ids_of_ranges_of_interest,
                                  stacking = "dense",
                                  col = 'white',
                                  size = 2)
    
    return(te_instance_track)
}


nanopore_track <- function(te, bam_file, max_reads = 25, name = "nanopore reads"){

    
    readRanges <- sub_sampling_reads(bam_file, 
                                     as.character(seqnames(te)), 
                                     start(te), 
                                     end(te),
                                     as.character(strand(te)), max_reads = max_reads)
    
    nanopore_reads <- Gviz::AnnotationTrack(
        range = readRanges,
        genome = "mm10",
        name = "reads",
        fill = "lightgray",
        chromosome = as.character(seqnames(te)),
        id = c()
    )
    
}


cage_seq_tracks <- function(te){
    
    bw <- rtracklayer::import.bw(cage_fwd_bw) %>%
        plyranges::filter_by_overlaps(as(paste0(
            as.character(seqnames(te)), ':', start(te), '-', end(te)
        ), "GRanges"))
    
    
    cage_limits = c(0, plyr::round_any(max(bw$score), 10, f = ceiling))
    
    cage_track <-  DataTrack(
        range = bw,
        type = "h",
        name = 'CAGE',
        col = "darkgreen",
        ylim = cage_limits
    )
    
    return(cage_track)
    
}

genomeBrowserNanopore <- function(te_island_id, 
                                  te, 
                                  te_island_annotation_file,
                                  bam_file){
    
    gtrack <- GenomeAxisTrack()
    te_island_track <- te_island_track(te_island_id_list = te_island_id, te_island_ranges, te_island_annotation_file)
    te_instance_track <- local_te_track(teRanges, te) #chromosome, start_coord, end_coord) 
    nanopore_track <- nanopore_track(te, bam_file) #start_coord, end_coord, strand, bam_file)
    cage_tracks <- cage_seq_tracks(te) #chromosome, start_coord, end_coord)
    
    Gviz::plotTracks(
        list(
            #ideoTrack,
            gtrack,
            te_island_track,
            te_instance_track,
            cage_tracks,
            nanopore_track
        ),
        #add35 = TRUE,
        shape = 'box',
        scale = 0.25,
        chromosome = chromosome,
        from = start_coord,
        to = end_coord,
        col.axis = "black",
        col.title = "black",
        cex.feature = 0.7,
        background.title = "white",
        fontcolor.legend = "black"
    )
    
}
