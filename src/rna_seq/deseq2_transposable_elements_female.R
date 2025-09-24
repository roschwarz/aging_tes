# ==============================================================================
# RNA-Seq Expression Analysis for Transposable Elements (TEs)
# ==============================================================================
# This script performs differential expression analysis of transposable elements
# (TEs) across various tissues in female using RNA-Seq count data processed via 
# salmon (see rna_seq_female.sh).
# It utilizes DESeq2 for normalization and statistical testing, supports parallel
# execution for scaleability, and merges results across tissues for downstream use.
#
# ------------------------------------------------------------------------------
# INPUT:
# - Count tables (SalmonTE output) for each tissue:
#   ./results/rna_seq/female/detector/EXPR.csv
#
# OUTPUT:
# - DESeq2 results and objects:
#   ./results/rna_seq/female/deseq2/dds_TE_instances_salmonTE.Rdata         : DESeq2 dds objects
#   ./results/rna_seq/femal/deseq2/deseq_TE_instances_salmonTE.Rdata     : DESeq2 result objects
#   ./results/tables/02_deseq_results_te_female.csv     : Merged results (CSV)
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
aging_tes::load_rna_seq_female_env()

if (!dir.exists(rna_seq_deseq_dir)) {
    dir.create(rna_seq_deseq_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------------------------

file_assignment <- read.csv('/misc/paras/data/rschwarz/projects/aging_tes/data/raw/rna_seq/female/FILENAMES', header = TRUE)

meta_data <- file_assignment %>% 
    tidyr::separate(new_name, into = c('id', 'tissue', 'age'), sep = "_", remove = FALSE) %>% 
    dplyr::select(-id) %>% 
    mutate(age_group = case_when(age == "18w" ~ 'young', .default = 'old'))

counts <- data.table::fread(counts_rna) %>% 
    filter(grepl("^chr", Name)) %>% 
    column_to_rownames('Name')

# ------------------------------------------------------------------------------
# Step 1: Run DESeq2 for each tissue if results do not exist
# ------------------------------------------------------------------------------

if (!file.exists(paste0(rna_seq_deseq_dir, deseq_results_te))) {
# Load required libraries for parallel processing
    library(BiocParallel)
    library(future.apply)
    
    register(MulticoreParam(workers = 4)) # BiocParallel
    plan(multicore, workers = 2)  # future apply
    
    dds_list <- future_lapply(c('brain', 'skin'), function(tis) {
        message("Running DESeq2 for ", tis)
        meta <- meta_data %>% 
            filter(tissue == tis, grepl("R1$", new_name))
        
        tissue_counts <- counts[,meta$new_name]
        
        col_data <- meta %>% dplyr::select(new_name, age_group) %>% column_to_rownames("new_name")
        
        dds <-  dds <- doDEseq(
            tissue_counts,
            col_data,
            design_formula = formula(~age_group),
            paral = TRUE,
            reference = 'young',
            target = 'te',
            relevel_condition = 'age_group'
            
        )
        
        return(dds)
        
    })
    
    names(dds_list) <- c('brain', 'skin')
    
    # Increase memory limit for large DESeq2 results
    options(future.globals.maxSize = 4 * 1024^3)
    
    res_list <- future_lapply(c('brain', 'skin'), function(tissue) {
        res <-  getDEseqResults(
            dds_list[[tissue]],
            coefficient = "age_group_old_vs_young",
            parallel = TRUE,
            FDR.filter = FALSE
        )
        
        return(res)
    })
    
    names(res_list) <- c('brain', 'skin')
        
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

aging_tes::load_annotations()
load_te_annotation()

deseq.te.merged <- do.call("rbind",
                           sapply(names(res_list), simplify = F, function(x) {
                               res_list[[x]] %>%
                                   rownames_to_column(var = "te_id") %>%
                                   mutate(tissue = x)
                               
                           }))

# Annotate merged results with TE metadata
deseq.te.merged <- merge(deseq.te.merged, te_annotation, by = 'te_id')

write.table(deseq.te.merged,
            file = paste0(table_dir, deseq_results_te_csv), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

