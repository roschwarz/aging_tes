
# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}



gene.exonRanges
geneRanges
genePromoter <- IRanges::promoters(geneRanges, upstream = 500, downstream = 0)

geneDown <- IRanges::flank(geneRanges, 500, start = FALSE)


blackRcloud::intersectGranger(teRanges, geneRanges) %>% data.frame() %>% nrow()


deseq_res <- read.csv(paste0(table_dir, "02_deseq_results_te_instances.csv"))

N_TES <- nrow(data.frame(teRanges))


getNumbers <- function(rangesOfInterest, tissue, rangeType){
    
    total_count <- nrow(data.frame(rangesOfInterest))
    
    promoter_count = blackRcloud::intersectGranger(rangesOfInterest, genePromoter) %>% 
        data.frame() %>% nrow()
    
    exon_count =  blackRcloud::intersectGranger(rangesOfInterest, gene.exonRanges) %>% 
        data.frame() %>% nrow()
    
    gene_count <- blackRcloud::intersectGranger(rangesOfInterest, geneRanges) %>% 
        data.frame() %>% nrow()
    
    intron_count = gene_count - exon_count
    
    down_count = blackRcloud::intersectGranger(rangesOfInterest, geneDown) %>% 
        data.frame() %>% nrow()
    
    intergenic_count = total_count - promoter_count - exon_count - intron_count - down_count
    
    counts <- data.frame(rangeType = rangeType, 
                         tissue = tissue,
                         total = total_count,
                         promoter = promoter_count,
                         exon =  exon_count,
                         intron = intron_count,
                         down = down_count,
                         intergenic = intergenic_count,
                         total_prop = total_count/N_TES*100,
                         promoter_prop = promoter_count/total_count*100,
                         exon_prop = exon_count/total_count*100,
                         intron_prop = intron_count/total_count*100,
                         down_prop = down_count/total_count*100,
                         intergenic_prop = intergenic_count/total_count *100
                         
    )
    
    return(counts)
    
    
}


############################
# Counts for expressed TEs #
############################

counts_ex_tes <- do.call('rbind', sapply(tissues, simplify = FALSE, function(x){
    
    df <- deseq_res %>% filter(tissue == x, !is.na(padj)) 
    
    ex_te_ranges <- GRanges(seqnames = df$chromosome,
                            ranges = IRanges(start = df$start,
                                             end = df$end,
                                             names = df$te_id),
                            strand = df$strand,
                            cluster.ID = df$te_id)
    
    
    counts <- getNumbers(rangeType = "expressed TEs", 
                         tissue = x, 
                         rangesOfInterest = ex_te_ranges)
    
    
    return(counts)
    
    
    
})
)

####################
# Counts for DETEs #
####################

counts_detes <- do.call('rbind', sapply(tissues, simplify = FALSE, function(x){
    
    df <- deseq_res %>% filter(tissue == x, padj <= 0.05) 
    
    ex_te_ranges <- GRanges(seqnames = df$chromosome,
                            ranges = IRanges(start = df$start,
                                             end = df$end,
                                             names = df$te_id),
                            strand = df$strand,
                            cluster.ID = df$te_id)
    
    
    counts <- getNumbers(rangeType = "DETEs", 
                         tissue = x, 
                         rangesOfInterest = ex_te_ranges)
    
    
    return(counts)
    
    
    
})
)

#######################
# Counts for CAGE_TEs #
#######################

cageRanges <- sapply(tissues, simplify = F, function(tissue){
    
    directory <- strsplit(cage.files[[tissue]], "/")[[1]][3]
    file <- paste0('results/cage_RS/',
                   directory,
                   '/raw_peaks/',
                   tissue,
                   '.cage_peaks.bed')
    
    #file = paste0('../data/cage/',tissue,'.peakAnnotation.bed')
    return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
    
})

te.Cages <- sapply(names(cageRanges), simplify = F, function(x) {
    
    df <- intersectGranger(cageRanges[[x]], teRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'cage',
                                              'subject' = 'te'))
    
    df <- df %>% filter(!duplicated(te.te.id))
    
    ranges <- GRanges(seqnames = df$te.seqnames,
                      ranges = IRanges(start = df$te.start,
                                       end = df$te.end,
                                       names = df$te.te.id),
                      strand = df$te.strand,
                      cluster.ID = df$te.te.id)
    
    return(ranges)
    
})

counts_cage_tes <- do.call('rbind', sapply(tissues, simplify = FALSE, function(x){
    
    te_cages <- te.Cages[[x]]
    
    
    counts <- getNumbers(rangeType = "CAGE TEs", 
                         tissue = x, 
                         rangesOfInterest = te_cages)
    
    
    return(counts)
    
    
    
})
)

########################
# Counts for QUANT_TEs #
########################

quantRanges <- sapply(tissues, simplify = F, function(tissue){
    
    directory <- strsplit(quant.files[[tissue]], "/")[[1]][3]
    file <- paste0('results/quant_RS/',
                   directory,
                   '/raw_peaks/',
                   tissue,
                   '.quant_peaks.bed')
    
    #file = paste0('../data/cage/',tissue,'.peakAnnotation.bed')
    return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
    
})

te.quants <- sapply(names(quantRanges), simplify = F, function(x) {
    
    df <- intersectGranger(quantRanges[[x]], teRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'cage',
                                              'subject' = 'te'))
    
    df <- df %>% filter(!duplicated(te.te.id))
    
    ranges <- GRanges(seqnames = df$te.seqnames,
                      ranges = IRanges(start = df$te.start,
                                       end = df$te.end,
                                       names = df$te.te.id),
                      strand = df$te.strand,
                      cluster.ID = df$te.te.id)
    
    return(ranges)
    
})

counts_quant_tes <- do.call('rbind', sapply(tissues, simplify = FALSE, function(x){
    
    te_quants <- te.quants[[x]]
    
    
    counts <- getNumbers(rangeType = "Quant TEs", 
                         tissue = x, 
                         rangesOfInterest = te_quants)
    
    
    return(counts)
    
    
    
})
)


count_data <- rbind(counts_ex_tes, counts_detes, rbind(counts_cage_tes, counts_quant_tes))


counts_for_word <- count_data %>% 
    mutate(total = paste0(total, " (", round(total_prop,2), ")"),
                      promoter = paste0(promoter, " (", round(promoter_prop,2), ")"),
                      exon = paste0(exon, " (", round(exon_prop,2), ")"),
                      intron = paste0(intron, " (", round(intron_prop,2), ")"),
                      down = paste0(down, " (", round(down_prop,2), ")"),
                      intergenic = paste0(intergenic, " (", round(intergenic_prop,2), ")")
                      ) %>% 
    dplyr::select(rangeType, tissue, total, promoter, exon, intron, down, intergenic) %>% view() 
