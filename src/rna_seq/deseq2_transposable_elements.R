# --------------------------------- Notes--- -----------------------------------
# This script is used for the RNA-Seq quantification
#
# Input:
#
# - Count tables of the respective tissue (EXPR.csv) in ./results/rna_seq/<tissue>/alignment_SalmonTE/EXPR.csv
# - te region file for deseq analysis for TE islands 
#
# Output:
#   - ./results/rna_seq/deseq2/
#

# TE instances
#
# dds.TE.Salmon.Rdata - DESeq2 dds object for TEs base on SalmonTE
# deseq.TE.SalmonTE.Rdata - DESeq2 results for TEs based on SalmonTE counts
# 02_deseq_results_te.csv - tissue merged csv file of DESeq results of TE
#   instances
#


if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}


if (!dir.exists(rna_seq_deseq_dir)) {
    dir.create(rna_seq_deseq_dir, recursive = TRUE)
}

# ----- Quantification of TE expression during aging in different tissues ------

# read the count tables and store them in a list
# path to count tables can be found in 01_load_environment
if (!file.exists(paste0(rna_seq_deseq_dir, deseq_results))) {
    
    # Setup for Parallel
    library(BiocParallel)
    register(MulticoreParam(workers = 4))
    
    library(future.apply)
    plan(multicore, workers = 3) 
    
    dds_list <- future_lapply(tissues, function(tissue) {
        message("Running DESeq2 for ", tissue)
        
        counts <- load_salmonTE_counts(counts_rna[[tissue]], "instance")[["instance"]] %>%
            tibble::column_to_rownames("TE")
        
        condition = getConditions(names(counts))
        
        dds <- doDEseq(
            counts,
            condition,
            paral = TRUE,
            reference = 'young',
            target = 'te'
        )
        
        return(dds)
    })
    
    names(dds_list) <- tissues
    
    
    options(future.globals.maxSize = 4 * 1024^3) # need to increase memory because it is limited to 500 mb
    
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
         file = paste0(rna_seq_deseq_dir, deseq_dds))
    save(res_list,
         file = paste0(rna_seq_deseq_dir, deseq_results))
    
}else{
    
    res_list <- loadRdata(paste0(rna_seq_deseq_dir, deseq_results))
}

#####
#
# Notiz
#
# update der Ordnerstruktur
# Letzter schritt die DESeq analyse parallelisiert
#
# Next steps:
#   - mache skripte fur gene un TE islands
#   - check ob du Volcanos hier mit ins skript intigierst
#   - schliesse Kapitel eins von deinem Manuskript ab mit allen benoetigten bilder und skripten
#   - commit to git




# Merge the results of the different tissues to get one huge table with all results
deseq.te.merged <- do.call("rbind",
                           sapply(names(res_list), simplify = F, function(x) {
                               res_list[[x]] %>%
                                   rownames_to_column(var = "te_id") %>%
                                   mutate(tissue = x)

                           }))

deseq.te.merged <- merge(deseq.te.merged, te.annotation, by = 'te_id')

write.table(deseq.te.merged,
            file = paste0(table_dir, '02_deseq_results_te_instances.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

rm(deseq.te.merged, deseq.te)


