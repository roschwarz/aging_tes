# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    source('./01_load_environment.R')
} else if (!ENVIRONMENT_LOADED) {
    
    source('./01_load_environment.R')
    
}

tpm.own.rna.exp <- loadRdata(paste0(data_dir, "rna_RS_tpms_TE_instances_SalmonTE.Rdata"))

deseq.own.rna.merged <- read.csv(paste0(table_dir,
                                        "02_deseq_results_te_instances.csv"))

chrom = 'chr1'
strand = '-'
start_coord = 21322152
end_coord = 21326292
#gtrack <- GenomeAxisTrack()

#teRanges <- blackRcloud::getRanges("mmusculus", 102, "te", 1)

TEs_of_interest = data.frame(teRanges) %>%
    filter(seqnames == chrom,
           strand == strand,
           dplyr::between(start, start_coord, end_coord)) %>%
    pull(te.id) %>%
    unique()


norm <- tpm.own.rna.exp$blood %>% filter(te_id %in% TEs_of_interest)


norm <- norm %>% 
    separate('sample', c("age.name", "age", "sample.id", "tissue"), sep = "[_]", remove = F) %>% 
    dplyr::select(-age.name)

norm <- norm %>% 
    mutate(te.annotation = paste0(order, '|', super_family, '|', family, '|', start)) %>% 
    dplyr::select(te.annotation, sample, TPM) %>% 
    spread(sample, TPM) %>%
    column_to_rownames(var = 'te.annotation')


norm <- as.matrix(norm)

norm <- t(apply(norm, 1, function(x) calc_z_score(x, capped = 2)))

# removes rows full of NA --> Why such rows are there?
norm <- norm[rowSums(is.na(norm)) != ncol(norm), ]

ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points") 
col_fun = colorRamp2(c(-2, 0, 2), c('#457b9d','white','#e63946'))

row_text_size = 6


Heatmap(norm,
        col = col_fun)
