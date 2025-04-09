# ==============================================================================
# RNA-Seq Expression Analysis for Transposable Elements (TEs)
# ==============================================================================
# This script performs differential expression analysis of transposable elements
# (TEs) across various tissues using RNA-Seq count data processed via SalmonTE.
# It utilizes DESeq2 for normalization and statistical testing, supports parallel
# execution for scalability, and merges results across tissues for downstream use.
#
# ------------------------------------------------------------------------------
# INPUT:
# - Count tables (SalmonTE output) for each tissue:
#   ./results/rna_seq/<tissue>/alignment_SalmonTE/EXPR.csv
#
# OUTPUT:
# - DESeq2 results and objects:
#   ./results/rna_seq/deseq2/dds_TE_instances_salmonTE.Rdata         : DESeq2 dds objects
#   ./results/rna_seq/deseq2/deseq_TE_instances_salmonTE.Rdata     : DESeq2 result objects
#   ./results/tables/02_deseq_results_te.csv     : Merged results (CSV)
#
# DEPENDENCIES:
# - Custom package: `aging_tes` (loaded via `devtools`)
# - Bioconductor: `BiocParallel`, `DESeq2`
# - Others: `future.apply`, `dplyr`, `tibble`
#
# ==============================================================================

# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_analysis_env()
aging_tes::load_rna_seq_env()

if (!dir.exists(rna_seq_deseq_dir)) {
    dir.create(rna_seq_deseq_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Step 1: Run DESeq2 for each tissue if results do not exist
# ------------------------------------------------------------------------------

if (!file.exists(paste0(rna_seq_deseq_dir, deseq_results_te))) {
    
    # Load required libraries for parallel processing
    library(BiocParallel)
    library(future.apply)
    
    register(MulticoreParam(workers = 4)) # BiocParallel
    plan(multicore, workers = 3)  # future apply
    
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
         file = paste0(rna_seq_deseq_dir, deseq_dds_te))
    save(res_list,
         file = paste0(rna_seq_deseq_dir, deseq_results_te))
    
}else{
    
    res_list <- loadRdata(paste0(rna_seq_deseq_dir, deseq_results_te))
}

# ------------------------------------------------------------------------------
# Step 2: Merge DESeq2 results across tissues into one table
# ------------------------------------------------------------------------------

deseq.te.merged <- do.call("rbind",
                           sapply(names(res_list), simplify = F, function(x) {
                               res_list[[x]] %>%
                                   rownames_to_column(var = "te_id") %>%
                                   mutate(tissue = x)

                           }))

# Annotate merged results with TE metadata
deseq.te.merged <- merge(deseq.te.merged, te.annotation, by = 'te_id')

write.table(deseq.te.merged,
            file = paste0(table_dir, deseq_results_te_csv), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')
