# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_analysis_env()
aging_tes::load_cage_seq_env()
aging_tes::load_annotations()

cageRanges <- load_cage_peak_annotation()
geneRanges <- load_gene_ranges()
load_te_island_annotation()

# ---------------------------- CAGE-peaks in TE islands ------------------------

te_island_cages <- sapply(names(cageRanges), 
                          simplify = F, 
                          USE.NAMES = T, 
                          function(x){
                              
                              df <-  intersectGranger(cageRanges[[x]], te_island_5primeRanges, tab = 'all')
                              names(df) <- str_replace_all(names(df), c('query' = 'cage', 'subject' = 'te_island'))
                              
                              return(df)
                              
                          })

deseq_cage <- loadRdata(cage_deseq_res)

te_island_gene_relation <- te_island_gene_rel()

deseq_cage_te_island_merged <- do.call('rbind', sapply(names(deseq_cage), simplify = F, function(x){
    
    df <- deseq_cage[[x]] %>%
        rownames_to_column(var = "peak_id") %>% 
        mutate(tissue = x)
    
    
    df <- merge(df, te_island_cages[[x]], by.x = 'peak_id', by.y = 'cage.names')
    
    df <- merge(df, te_island_gene_relation, 
                by.x = 'te_island.names', 
                by.y = 'te_island_id',
                all.x = T)
    
    return(df)
    
}))


write.table(deseq_cage_te_island_merged, 
            file = paste0(table_dir, 'deseq_cage_te_island.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')
