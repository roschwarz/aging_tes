
# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_analysis_env()
aging_tes::load_rna_seq_env()

if (!dir.exists(rna_seq_deseq_dir)) {
    dir.create(rna_seq_deseq_dir, recursive = TRUE)
}


if (!file.exists(paste0(rna_seq_deseq_dir, deseq_results_mixed))) {
    # Load required libraries for parallel processing
    library(BiocParallel)
    library(future.apply)
    
    register(MulticoreParam(workers = 4)) # BiocParallel
    plan(multicore, workers = 3)  # future apply
    
    dds_list <- future_lapply(tissues, function(tissue) {
        message("Running DESeq2 for ", tissue)
        
        counts <- fread(counts_rna[[tissue]]) %>%
            tibble::column_to_rownames("TE")
        
        condition = getConditions(names(counts))
        
        dds <- doDEseq(
            counts,
            condition,
            paral = TRUE,
            reference = 'young',
            target = 'all'
        )
        
        return(dds)
    })
    
    names(dds_list) <- tissues
    
    # Increase memory limit for large DESeq2 results
    options(future.globals.maxSize = 4 * 1024^3)
    
    res_list <- future_lapply(tissues, function(tissue) {
        res <-  getDEseqResults(
            dds_list[[tissue]],
            coefficient = "condition_old_vs_young",
            parallel = TRUE,
            FDR.filter = FALSE
        )
        
        return(res)
    })
    
    names(res_list) <- tissues
    
    save(dds_list,
         file = paste0(rna_seq_deseq_dir, deseq_dds_mixed))
    save(res_list,
         file = paste0(rna_seq_deseq_dir, deseq_results_mixed))
    
}else{
    
    res_list <- loadRdata(paste0(rna_seq_deseq_dir, deseq_results_mixed))
}

# ------------------------------------------------------------------------------
# Step 2: Merge DESeq2 results across tissues into one table
# ------------------------------------------------------------------------------

deseq.mixed.merged <- do.call("rbind",
                           sapply(names(res_list), simplify = F, function(x) {
                               res_list[[x]] %>%
                                   rownames_to_column(var = "te_id") %>%
                                   mutate(tissue = x)
                               
                           }))

# Annotate merged results with TE metadata
write.table(deseq.mixed.merged,
            file = paste0(table_dir, deseq_results_mixed_csv), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

