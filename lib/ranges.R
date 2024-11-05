
addPeakCoordinate <- function(countTable, DESeqRes){
    
    # Coordinates of each peak are stored in the featureCount table
    # The coordinates are added to the DESeq results and the DEseq results are
    # returned.
    coordinates <- as.data.frame(fread(countTable))
    
    coordinates <- coordinates %>% dplyr::select(Geneid , Chr, Start, End, Strand)
    
    DESeqRes <- rownames_to_column(as.data.frame(DESeqRes), 'Geneid')
    DESeqRes <- merge(DESeqRes, coordinates, by='Geneid', all.x = T)
    
    return(DESeqRes)
    
}

createCageGRange <- function(df){
    
    cagePeaks <- GRanges(seqnames = df$Chr,
                         ranges = IRanges(df$Start,
                                          end = df$End,
                                          names = df$Geneid),
                         strand = df$Strand,
                         peak_name = df$Geneid)
                         #baseMean = df$baseMean)
    
    return(cagePeaks)
    
}


getQueryOverlapInformation <- function(query, subject, info.column){
    # Do an overlap between two GRange objects and returns the column of
    # interest as a vector
    overlap <- query[queryHits(findOverlaps(query, subject))]
    
    overlap <- as.data.frame(overlap, row.names = NULL)
    
    return(unique(overlap[[info.column]]))
}


extendTERanges <- function(teRanges, geneRanges, gene.exonRanges){
    # This function characterizes the position of the TEs with respect to
    # genes. Three groups are possible for position - exonic, intronic, intergenic
    
    te.gene.intersection <-  intersectGranger(teRanges,
                                              geneRanges,
                                              'all') %>% 
        dplyr::select(query.te.id, subject.ensembl_gene_id, subject.external_gene_name)
    
    names(te.gene.intersection) <- c('te.id', 'ensembl_gene_id','external_gene_name')
    
    te.exon.intersection <-  intersectGranger(teRanges,
                                              gene.exonRanges,
                                              'all') %>% 
        dplyr::pull(query.te.id)
    
    
    teRanges <- transGrange(teRanges)
    
    teRanges <- teRanges %>% 
        mutate(position = case_when(te.id %in% te.exon.intersection ~ 'exonic',
                                    !(te.id %in% te.gene.intersection$te.id) ~ 'intergenic',
                                    TRUE ~ 'intronic'))
    
    teRanges <- merge(teRanges, te.gene.intersection, by = "te.id", all.x = TRUE )
    
    teRanges <-  GRanges(seqnames = teRanges$seqnames,
                         ranges = IRanges(teRanges$start,
                                          end = teRanges$end,
                                          names = teRanges$te.id),
                         strand = teRanges$strand,
                         te.id = teRanges$te.id,
                         position = teRanges$position,
                         ensembl_gene_id = teRanges$ensembl_gene_id,
                         external_gene_name = teRanges$external_gene_name)
    
    return(teRanges)
}



# Returns the promoter region of an annotation
# the data frame needs a start and end column
# you can also set a size of your promoter (default = 500).
getPromoter <- function(df, size = 500){
    
    df %>%  
    mutate(promoter_start = case_when(strand == "+" ~ as.numeric(start) - size, .default = as.numeric(end)),
           promoter_end = case_when(strand == "+" ~ as.numeric(start), .default = as.numeric(end) + size))
    
    
}
