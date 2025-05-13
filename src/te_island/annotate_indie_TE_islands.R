# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}


aging_tes::load_rna_seq_env()
aging_tes::load_cage_peak_annotation()
aging_tes::load_annotations()
aging_tes::load_te_island_env()


load_te_ranges()

te_island_5_prime_extended <- readGeneric('data/processed/annotation/mm10_TE_island_5prime_extended.bed',
                                          strand = 6,
                                          meta.cols = list(names = 4))

# Collect TE instances of TE-regions to get information about the composition of the TE
# regions later in the analysis.
te_island_instances <- blackRcloud::intersectGranger(te_island_5_prime_extended, teRanges, tab = 'all')

# TEs can intersect with multiple genes, so that they occur multiple times in
# the teRange object. Therefore a unique is applied at the end of the pipeline
te_island_instances <- te_island_instances %>% 
    dplyr::select(query.names, subject.te.id, subject.position) %>% 
    filter(!duplicated(subject.te.id)) %>% 
    unique()

names(te_island_instances) <- c("te_region_id", "te_id", "te_position")


#####################################
# individually expressed TE islands #
#####################################

# Usage of CAGE- and RNA-seq to get individually expressed TE islands. 
#  1. Intersect the extended TE islands (te_island_5_prime_extended) with CAGE-peaks to get potential 
#     individually expressed TE islands.
#  2. Get individually expressed TE islands that also contain at least one expressed TE instance that 
#     is detected via RNA-Seq. An TE instance is considered as expressed when an adjusted p-value after
#     the DESeq analysis is available.

# Cage intersection
indie_te_islands <- sapply(names(cageRanges), 
                          simplify = F, 
                          USE.NAMES = T, 
                          function(x){
                              
                              blackRcloud::intersectGranger(cageRanges[[x]], te_island_5_prime_extended, tab = 'all')
                              
                          })

# RNA intersection
load_deseq_results() # loads rna_deseq_res_te

expressed_te_islands <- sapply(names(rna_deseq_res_te), simplify = F, function(x){
    
    expressed_TEs_rna <- rownames(rna_deseq_res_te[[x]]  %>% filter(!is.na(padj)))
    
    exp_islands <- te_island_instances %>% 
        filter(te_id %in% expressed_TEs_rna) %>% 
        dplyr::pull('te_region_id')
    
    # regions Ids occur multiple times because of the composition of the islands
    # However, only one id is needed to get the expressed TE island 
    return(unique(exp_islands))
})

# filter the cage intersecting TE islands for those that intersect with 
# expressed TEs detected via RNA-Seq
indie_te_islands <- sapply(names(indie_te_islands), simplify = F, function(x){
    
    te_island <- indie_te_islands[[x]]
    exp.te <- expressed_te_islands[[x]]
    
    te_island %>% 
        filter(subject.names %in% exp.te) %>% 
        filter(!duplicated(subject.names))
    
})


## Write bed-files of individually expressed TE islands
sapply(names(indie_te_islands), function(x){
    
    df <- indie_te_islands[[x]]
    
    bed_file <- GRanges(seqnames = df$subject.seqnames,
                        ranges = IRanges(start = df$subject.start,
                                         end = df$subject.end,
                                         names = df$subject.names),
                        strand = df$subject.strand,
                        te_island_id = df$subject.names,
                        value = 1)
    
    bed_file <- transGrange(bed_file)
    
    bed_file <- bed_file[c('seqnames', 'start', 'end', 'te_island_id',  'value', 'strand')]
    
    filename = indie_te_island_bed[[x]]
    
    write.table(bed_file, file = filename, sep = '\t', col.names = F, row.names = F, quote = F)
    
})
