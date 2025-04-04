if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}

getConditions <- function(file_names) {
    # Project specific. It uses the naming convention of the fastq file to
    # determine the tissue and age of each sample.
    
    
    samples <- data.frame(
        SampleID = file_names,
        condition = file_names)
    
    # define the conditions for each sample
    if (grepl("skin", file_names[1], ignore.case = T)) {
        samples <- samples %>%
            separate(condition,
                     c("ID", "age"),
                     sep = "[_]") %>%
            mutate(age = str_replace(age, "[0-9]*Skin", "")) %>%
            mutate(condition = as.factor(case_when(
                age == "MO" ~ "old",
                age == "MY" ~ "young"
            ))) %>%
            dplyr::select(SampleID, condition)
        
        return(samples)
    }
    
    samples <- samples %>%
        separate(condition,
                 c("ID", "org", "age", "tissue", "number"),
                 sep = "[_]") %>%
        mutate(condition = as.factor(case_when(age == "o" ~ "old",
                                               age == "y" ~ "young"))) %>%
        dplyr::select(SampleID, condition)
    
    return(samples)
    
    
}

counts <- fread('results/TEItx/teitx_salmonTE_mapped/EXPR.csv')
counts <- fread('results/TEItx/teitx_salmonTE_stranded_mapped/EXPR.csv')

teitx_counts <- counts %>% 
    filter(grepl('^TE', TE)) %>% 
    column_to_rownames(var = 'TE')
    

skin_samples <- names(teitx_counts)[grepl("skin", names(teitx_counts), ignore.case = T)]
brain_samples <- names(teitx_counts)[grepl("brain", names(teitx_counts), ignore.case = T)]
blood_samples <- names(teitx_counts)[grepl("blood", names(teitx_counts), ignore.case = T)]

#########
# Brain #
#########

brain_counts <- teitx_counts[brain_samples]

condition <- getConditions(names(brain_counts))

dds <- DESeqDataSetFromMatrix(countData = round(brain_counts),
                              colData = condition,
                              design = ~condition)


dds <- dds[rowSums(counts(dds)) >= 10, ]

dds$condition <- relevel(dds$condition, ref = 'young')
dds <- DESeq(dds)

deseq_res <- data.frame(results(dds))


brain.volcano <- volcanoPlot(deseq_res, FDR = 0.05) +
    theme_rob(base_size = 12, base_family = 'Arial') +
    theme(legend.position = 'None')




# Get Promoter region of top 500 down-regulated TE regions


brain_promoter_top500_down <- deseq_res %>% 
    filter(!is.na(padj), log2FoldChange < 0) %>% 
    slice_min(order_by = padj, n = 500) %>% 
    rownames_to_column("id") %>% 
    separate(id, 
             into = c('te_region_id', 
                      'NA',
                      'chromosome',
                      'coord',
                      'strand'),
             sep = "[:(]") %>%
    separate(coord, into = c('start', 'end'), sep = "[-]") %>%
    mutate(across('strand', str_replace, "[)]", "")) %>% 
    getPromoter() %>% 
    dplyr::select(chromosome, promoter_start, promoter_end, te_region_id, baseMean, strand)
    

writeBedFile(brain_promoter_top500_down, file.name = "results/TEItx/homer/brain_promoter_top500_down.bed")


########
# Skin #
########

skin_counts <- teitx_counts[skin_samples]

skin_condition <- getConditions(names(skin_counts))

skin_dds <- DESeqDataSetFromMatrix(countData = round(skin_counts),
                              colData = skin_condition,
                              design = ~condition)


skin_dds <- skin_dds[rowSums(counts(skin_dds)) >= 10, ]

skin_dds$condition <- relevel(skin_dds$condition, ref = 'young')
skin_dds <- DESeq(skin_dds)

skin_deseq_res <- data.frame(results(skin_dds))

skin.volcano <- volcanoPlot(skin_deseq_res, FDR = 0.05) +
    theme_rob(base_size = 12, base_family = 'Arial') +
    theme(legend.position = 'None')


# Get Promoter region of top 500 down-regulated TE regions

skin_promoter_top500_down <- skin_deseq_res %>% 
    filter(!is.na(padj), log2FoldChange < 0) %>% 
    slice_min(order_by = padj, n = 500) %>% 
    rownames_to_column("id") %>% 
    separate(id, 
             into = c('te_region_id', 
                      'NA',
                      'chromosome',
                      'coord',
                      'strand'),
             sep = "[:(]") %>%
    separate(coord, into = c('start', 'end'), sep = "[-]") %>%
    mutate(across('strand', str_replace, "[)]", "")) %>% 
    getPromoter() %>% 
    dplyr::select(chromosome, promoter_start, promoter_end, te_region_id, baseMean, strand)


writeBedFile(skin_promoter_top500_down, file.name = "results/TEItx/homer/skin_promoter_top500_down.bed")

########
# blood #
########

blood_counts <- teitx_counts[blood_samples]

blood_condition <- getConditions(names(blood_counts))

blood_dds <- DESeqDataSetFromMatrix(countData = round(blood_counts),
                              colData = blood_condition,
                              design = ~condition)


blood_dds <- blood_dds[rowSums(counts(blood_dds)) >= 10, ]

blood_dds$condition <- relevel(blood_dds$condition, ref = 'young')
blood_dds <- DESeq(blood_dds)

blood_deseq_res <- data.frame(results(blood_dds))


blood.volcano <- volcanoPlot(blood_deseq_res, FDR = 0.05) +
    theme_rob(base_size = 12, base_family = 'Arial') +
    theme(legend.position = 'None')


# Get Promoter region of top 500 down-regulated TE regions


blood_promoter_top500_down <- blood_deseq_res %>% 
    filter(!is.na(padj), log2FoldChange < 0) %>% 
    slice_min(order_by = padj, n = 500) %>% 
    rownames_to_column("id") %>% 
    separate(id, 
             into = c('te_region_id', 
                      'NA',
                      'chromosome',
                      'coord',
                      'strand'),
             sep = "[:(]") %>%
    separate(coord, into = c('start', 'end'), sep = "[-]") %>%
    mutate(across('strand', str_replace, "[)]", "")) %>% 
    getPromoter() %>% 
    dplyr::select(chromosome, promoter_start, promoter_end, te_region_id, baseMean, strand)


writeBedFile(brain_promoter_top500_down, file.name = "results/TEItx/homer/blood_promoter_top500_down.bed")

ggarrange(brain.volcano, skin.volcano, blood.volcano, nrow =1)


