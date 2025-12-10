# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_analysis_env()
aging_tes::load_quant_seq_env()

# ------------------------------------------------------------------------------
# Quant-Seq DESeq2 of Peaks
# ------------------------------------------------------------------------------

if (!file.exists(paste0(quant_deseq_dir, quant_deseq_res))) {
    
    data <- sapply(names(quant_counts),
                   simplify = F,
                   function(x) getFeatureCountTab(quant_counts[[x]], filter = F))
    
    count_tables <- sapply(names(data),
                           function(x) updateHeader(data[[x]][['counts']]))
    
    conditions <- sapply(names(count_tables),
                         simplify = F,
                         function(x) getCondition_quant(names(count_tables[[x]])))
    
    
    # runDESeq
    dds_quant <- sapply(names(count_tables), simplify = F,  function(x){
        
        doDEseq(count.matrix = count_tables[[x]],
                col_data = conditions[[x]],
                reference = 'young',
                target = 'all'
        )
    })
    
    # collect DESeq results
    deseq_quant <- sapply(names(dds_quant), simplify = F, function(x){
        
        getDEseqResults(dds_quant[[x]], coefficient = 'condition_old_vs_young', FDR.filter = F)
        
    })
    
    save(dds_quant, file = quant_dds)
    save(deseq_quant, file = quant_deseq_res)
    
    rm(dds_quant)
}else{
    
    deseq_quant <- loadRdata(quant_deseq_res)
}

deseq_quant_merged <- do.call('rbind',
                             sapply(names(deseq_quant),
                                    simplify = F,
                                    function(x) {
                                        df <- deseq_quant[[x]] %>%
                                            rownames_to_column(var = "peak_id") %>%
                                            mutate(tissue = x)
                                        return(df)
                                    }))

write.csv(deseq_quant_merged,
          file = paste0(quant_deseq_dir, quant_deseq_res_csv),
          row.names = F,
          col.names = T,
          quote = F,
          sep = ",")


