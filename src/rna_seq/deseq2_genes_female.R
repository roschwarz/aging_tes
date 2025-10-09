# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_analysis_env()
aging_tes::load_annotations()
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


counts <- load_salmonTE_counts(counts_rna, request = "gene")[["gene"]]

# ------------------------------------------------------------------------------
# Step 1: Run DESeq2 for each tissue if results do not exist
# ------------------------------------------------------------------------------

if (!file.exists(paste0(rna_seq_deseq_dir, deseq_results_gene))) {
    
    # Load required libraries for parallel processing
    library(BiocParallel)
    library(future.apply)
    
    register(MulticoreParam(workers = 4)) # BiocParallel
    plan(multicore, workers = 3)  # future apply
    
    dds_list <- future_lapply(c('brain', 'skin'), function(tis) {
        message("Running DESeq2 for ", tis)
        
        meta <- meta_data %>% 
            filter(tissue == tis, grepl("R1$", new_name))
        
        tissue_counts <- counts[,meta$new_name]
        
        col_data <- meta %>% dplyr::select(new_name, age_group) %>% column_to_rownames("new_name")
        
        dds <- doDEseq(
            counts,
            col_data = col_data,
            design_formula = formula(~age_group),
            paral = TRUE,
            reference = 'young',
            target = 'gene',
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

load_gene_ranges()

deseq.gene.merged <- merge(deseq.gene.merged,
                           transGrange(geneRanges),
                           by = 'ensembl_gene_id')

write.table(deseq.gene.merged,
            file = paste0(table_dir, deseq_results_gene_csv), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')
