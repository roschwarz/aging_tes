
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


local_te_track <- function(te_ranges, chromosome, start_coord, end_coord){
    
    local_tes = data.frame(te_ranges) %>% 
        filter(seqnames == chromosome, 
               strand == strand,
               dplyr::between(start, start_coord, end_coord)) %>% 
        blackRcloud::splitTEID("te.id") %>% 
        unique()
        
        #pull(te.id) %>% 
        # unique()
    
    #te_ids <- blackRcloud::splitTEID()
    
    te_instance_track <- AnnotationTrack(te_ranges[te_ranges@elementMetadata$te.id %in% local_tes$te.id & !duplicated(te_ranges@elementMetadata$te.id)],
                                  featureAnnotation = "id",
                                  #fontcolor.item = "white",
                                  id = local_tes$super_family,
                                  stacking = "dense",
                                  col = 'white',
                                  size = 2)
    
    return(te_instance_track)
}


nanopore_track <- function(chromosome, bam_file, name = "nanopore reads"){
    
    nanopore_reads <- Gviz::AnnotationTrack(
        range = bam_file,
        genome = "mm10",
        name = "reads",
        fill = "lightgray",
        chromosome = chromosome,
        id = c()
    )
    
}


cage_seq_tracks <- function(chromosome, start_coord, end_coord){
    
    bw <- rtracklayer::import.bw(cage_fwd_bw) %>%
        plyranges::filter_by_overlaps(as(paste0(
            chromosome, ':', start_coord, '-', end_coord
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
                                  chromosome,
                                  start_coord, 
                                  end_coord, 
                                  te_island_annotation_file,
                                  bam_file){
    
    gtrack <- GenomeAxisTrack()
    
    te_island_track <- te_island_track(te_island_id_list = te_island_id, te_island_ranges, te_island_annotation_file)
    te_instance_track <- local_te_track(teRanges, chromosome, start_coord, end_coord) 
    nanopore_track <- nanopore_track(chromosome, bam_file)
    cage_tracks <- cage_seq_tracks(chromosome, start_coord, end_coord)
    
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
