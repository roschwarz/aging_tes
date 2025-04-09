# ==============================================================================
# RNA-Seq Gene-Level Differential Expression Analysis
# ==============================================================================
# This script performs differential expression analysis of genes across multiple
# tissues using RNA-Seq data. It leverages DESeq2 for statistical inference and
# supports parallel processing for performance optimization.
#
# ------------------------------------------------------------------------------
# INPUT:
# - Gene-level count tables (SalmonTE output):
#   ./results/rna_seq/<tissue>/alignment_SalmonTE/EXPR.csv
#
# OUTPUT:
# - DESeq2 results and objects:
#   ./results/rna_seq/deseq2/dds.TE.Gene.Rdata            : DESeq2 dds objects
#   ./results/rna_seq/deseq2/deseq.TE.GeneResults.Rdata   : DESeq2 result objects
#   ./results/rna_seq/deseq2/02_deseq_results_gene.csv    : Merged results (CSV)
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
aging_tes::load_annotations()
aging_tes::load_rna_seq_env()

# ------------------------------------------------------------------------------
# Step 1: Run DESeq2 for each tissue if results do not exist
# ------------------------------------------------------------------------------

if (!file.exists(paste0(rna_seq_deseq_dir, deseq_results_gene))) {
    
    # Load required libraries for parallel processing
    library(BiocParallel)
    library(future.apply)
    
    register(MulticoreParam(workers = 4)) # BiocParallel
    plan(multicore, workers = 3)  # future apply
    
    dds_list <- future_lapply(tissues, function(tissue) {
        message("Running DESeq2 for ", tissue)
        
        counts <- load_salmonTE_counts(counts_rna[[tissue]], "gene")[["gene"]]
        
        condition = getConditions(names(counts))
        
        dds <- doDEseq(
            counts,
            condition,
            paral = TRUE,
            reference = 'young',
            target = 'gene'
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
         file = paste0(rna_seq_deseq_dir, deseq_dds_gene))
    save(res_list,
         file = paste0(rna_seq_deseq_dir, deseq_results_gene))
    
}else{
    
    res_list <- loadRdata(paste0(rna_seq_deseq_dir, deseq_results_gene))
}

# ------------------------------------------------------------------------------
# Step 2: Merge DESeq2 results across tissues into one table
# ------------------------------------------------------------------------------

deseq.gene.merged <- do.call("rbind",
                           sapply(names(res_list), simplify = F, function(x) {
                               res_list[[x]] %>%
                                   rownames_to_column(var = "ensembl_gene_id") %>%
                                   mutate(tissue = x)
                               
                           }))

deseq.gene.merged <- merge(deseq.gene.merged,
                            transGrange(geneRanges),
                            by = 'ensembl_gene_id')

write.table(deseq.gene.merged,
            file = paste0(table_dir, deseq_results_gene_csv), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')
