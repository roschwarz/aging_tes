# ==============================================================================
# RNA-Seq Expression Analysis for Transposable Elements (TEs)
# ==============================================================================
#
# This script performs a differential expression analysis of transposable elements
# (TEs) in public data using RNA-Seq count data mapped via salmon (see rna_seq_public.sh).
# It utilizes DESeq2 for normalization and statistical testing and supports parallel
# execution to improve scalability.
#
# The public dataset contains samples of liver, gastroncnemious muscle, and white adipose tissue
# from female and male mice of different ages (6 and 24 months).
#
# The count tables are loaded and pre-filtered by removing transposable elements with less than 11 reads
# in sum across all samples.
#
# ------------------------------------------------------------------------------
# INPUT:
# - Count tables (SalmonTE output) for female and male
# 
#
# OUTPUT:
# - DESeq2 results and objects:
#   -   dds_TE_instances_salmonTE_<tissue>_<sex>.Rdata
#   -   deseq_TE_instances_salmonTE_<tissue>_<sex>.Rdata
# - csv-file of DESeq2 results (all tissues and sexes merged)
#    - deseq_results_te_instances_public.csv

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
aging_tes::load_rna_seq_public_data_env()

if (!dir.exists(rna_seq_deseq_dir)) {
    dir.create(rna_seq_deseq_dir, recursive = TRUE)
}

public_tissues = c("Liver", "Gastrocnemius muscle", "White adipose tissue")

# Load required libraries for parallel processing
library(BiocParallel)
library(future.apply)

register(MulticoreParam(workers = 4)) # BiocParallel
plan(multicore, workers = 6)  # future apply

# Increase memory limit for large DESeq2 results
options(future.globals.maxSize = 4 * 1024^3)
# double loop for sex (female and male) and tissue (liver, Gastroncnemious muscle, White adipose tissue)
deseq_sex_results <- future_lapply(names(counts_rna), function(sex){
    message("Running DESeq2 for ", sex)
    
    counts <- data.table::fread(counts_rna[[sex]])
    
    file_names <- c("TE", 
                    str_replace(basename(dirname(names(counts))), "_1.fastq.gz", "")[-1])
    
    names(counts) <- file_names
    
    counts <- counts %>% 
        column_to_rownames("TE")
    
    deseq_tissue_results <- future_lapply(public_tissues, function(t){
        
        message("Running DESeq2 for ", sex, " and ", t)
        
        coldata <- read.csv(meta_data[[sex]]) %>% 
            dplyr::select(Run, AGE, tissue) %>% 
            filter(tissue == t) %>% 
            tibble::column_to_rownames('Run')
        
        tissue_counts <- counts[, rownames(coldata)]
        
        dds <- doDEseq(count.matrix = tissue_counts,
                col_data = coldata,
                design_formula = ~AGE,
                paral = TRUE,
                reference = '6m',
                count_threshold = 10,
                target = 'te',
                relevel_condition = "AGE")
        
        res <- results(dds)
        
        dds_file_name = paste0("dds_TE_instances_salmonTE_", t, "_", sex, ".Rdata")
        res_file_name = paste0("deseq_TE_instances_salmonTE_", t, "_", sex, ".Rdata")
        
        save(dds,
             file = paste0(rna_seq_deseq_dir, dds_file_name))
        save(res,
             file = paste0(rna_seq_deseq_dir, res_file_name))
        
        res <- data.frame(res) %>% 
            mutate(sex = sex,
                   tissue = t)
        
        return(res)
    })
    
    names(deseq_tissue_results) <- public_tissues
    
    # Combine the different tissues in one data frame
    deseq_tissue_results_merged <- do.call("rbind",
            lapply(names(deseq_tissue_results), function(x) {
                
                deseq_tissue_results[[x]] %>% 
                    rownames_to_column("te_id")
            })
            ) %>% 
        blackRcloud::splitTEID("te_id")
        
    return(deseq_tissue_results_merged)
})

names(deseq_sex_results) <- names(counts_rna)

# Combine male in female in one table and store it as a csv file.
public_deseq_res <- do.call("rbind", 
                            lapply(names(deseq_sex_results), function(x) return(deseq_sex_results[[x]]))
                            )

write.csv(public_deseq_res,
          file = paste0(rna_seq_deseq_dir, "deseq_results_te_instances_public.csv"))
