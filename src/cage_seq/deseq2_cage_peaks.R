# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_analysis_env()
aging_tes::load_cage_seq_env()


# ------------------------------------------------------------------------------
# Quantification/DESeq2 of CAGE-peaks
# ------------------------------------------------------------------------------
if (!file.exists(paste0(cage_deseq_dir, cage_deseq_res))) {
    
    data <- sapply(names(cage_counts),
                   simplify = F,
                   function(x) getFeatureCountTab(cage_counts[[x]], filter = F))
    
    count_tables <- sapply(names(data),
                           function(x) updateHeader(data[[x]][['counts']]))
    
    conditions <- sapply(names(count_tables),
                         simplify = F,
                         function(x) getCondition_cage(names(count_tables[[x]])))
    
    
    # runDESeq
    dds.cage <- sapply(names(count_tables), simplify = F,  function(x){
        
        doDEseq(count.matrix = count_tables[[x]],
                col_data = conditions[[x]],
                reference = 'young',
                target = 'all'
        )
    })
    
    # collect DESeq results
    deseq.cage <- sapply(names(dds.cage), simplify = F, function(x){
        
        getDEseqResults(dds.cage[[x]], coefficient = 'condition_old_vs_young', FDR.filter = F)
        
    })
    
    save(dds.cage, file = cage_dds)
    save(deseq.cage, file = cage_deseq_res)
    
    rm(dds.cage)
    
}else{
    
    deseq.cage <- loadRdata(cage_deseq_res)
}

deseq_cage_merged <- do.call('rbind',
                             sapply(names(deseq.cage),
                                    simplify = F,
                                    function(x) {
                                        df <- deseq.cage[[x]] %>%
                                            rownames_to_column(var = "peak_id") %>%
                                            mutate(tissue = x)
                                        
                                        return(df)
                                        
                                    }))

write.table(deseq_cage_merged, 
            file = paste0(table_dir, cage_deseq_res_csv), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


# ------------------------------------------------------------------------------
# CAGE Peaks in genes
# ------------------------------------------------------------------------------
#
# If cage annotation is not available, run cage_store_annotation.R

cageRanges <- load_cage_peak_annotation()

aging_tes::load_annotations()

gene_cages <- sapply(names(cageRanges), simplify = F, function(x) {
    
    df <- intersectGranger(cageRanges[[x]], geneRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'cage',
                                              'subject' = 'gene'))
    
    return(df)
    
})

deseq_gene_cage_merged <- do.call('rbind', sapply(names(deseq.cage), simplify = F, function(x){
    
    df <- deseq.cage[[x]] %>%
        rownames_to_column(var = "peak_id") %>% 
        mutate(tissue = x)
    
    df <- merge(df, gene_cages[[x]], by.x = 'peak_id', by.y = 'cage.names') 
    
    return(df)
    
}))

write.table(deseq_gene_cage_merged, 
            file = paste0(table_dir, '07_deseq_cage_gene.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

# ------------------------------------------------------------------------------
# CAGE Peaks in transposable elements
# ------------------------------------------------------------------------------
#
# If cage annotation is not available, run cage_store_annotation.R

te_cages <- sapply(names(cageRanges), simplify = F, function(x) {
    
    df <- intersectGranger(cageRanges[[x]], teRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'cage',
                                              'subject' = 'te'))
    
    return(df)
    
})

deseq_te_cage_merged <- do.call('rbind',
                                sapply(names(deseq.cage),
                                       simplify = F,
                                       function(x) {
                                           df <- deseq.cage[[x]] %>%
                                               rownames_to_column(var = "peak_id") %>%
                                               mutate(tissue = x)
                                           
                                           df <-
                                               merge(df, te_cages[[x]], by.x = 'peak_id', by.y = 'cage.names')
                                           
                                       }))

write.table(deseq_te_cage_merged, 
            file = paste0(table_dir, '07_deseq_cage_te.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

