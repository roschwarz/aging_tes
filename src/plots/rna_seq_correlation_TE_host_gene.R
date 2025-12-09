# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_annotations()
aging_tes::load_plotting_env()

te_annotation <- load_te_annotation()

intragenic_tes <- te_annotation %>% 
    filter(position == 'intronic') %>% 
    dplyr::select(te_id, ensembl_gene_id)

gene_deseq <- fread('results/tables/02_deseq_results_genes.csv')
te_instances_deseq <- fread('results/tables/02_deseq_results_te_instances.csv')

intra_combination <-  do.call('rbind', sapply(tissues, simplify = F, function(x){
    
    genes <- gene_deseq %>% filter(tissue == x, 
                                         ensembl_gene_id %in% intragenic_tes$ensembl_gene_id) %>% 
        dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
        dplyr::rename(gene_log2FoldChange = log2FoldChange)
    
    tes <- te_instances_deseq %>% filter(tissue == x,
                                               te_id %in% intragenic_tes$te_id) %>% 
        dplyr::select(te_id, log2FoldChange, ensembl_gene_id) %>% 
        dplyr::rename(te_log2FoldChange = log2FoldChange)


    combination <- merge(genes, tes, by = 'ensembl_gene_id')
    combination$position <- 'intra'
    combination$distance <- 0
    
    return(combination)
}))

intra_combination$tissue <- factor(intra_combination$tissue, level = tissues)


#########################################################################
##### 
### Correlation between intergenic TEs and their associated (closest)
##  genes
#
#########################################################################


te_gene_distance <- read.csv(paste0(table_dir,
                                    "intergenic_te_gene_distance.csv"))

inter_combination <-  do.call('rbind', sapply(tissues, simplify = F, function(x){
    
    inter_te_deseq <-
        te_instances_deseq %>% filter(tissue == x,
                                      te_id %in% te_gene_distance$te_id,
        ) %>%
        dplyr::select(te_id, log2FoldChange, tissue) %>%
        dplyr::rename(te_log2FoldChange = log2FoldChange)
    
    inter_te_deseq <-
        merge(inter_te_deseq, te_gene_distance, by = 'te_id')
    
    inter_gene_deseq <- gene_deseq %>% filter(tissue == x,
                                                    ensembl_gene_id %in% te_gene_distance$ensembl_gene_id,
    ) %>%
        dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name) %>%
        dplyr::rename(gene_log2FoldChange = log2FoldChange)
    
    inter_combination <-
        merge(inter_te_deseq, inter_gene_deseq, by = 'ensembl_gene_id')


    inter_combination$position <- 'inter'
    inter_combination$distance <- 0
    
    return(inter_combination)
}))

new_order <- names(intra_combination)

inter_combination <- inter_combination[,..new_order]

data <- rbind(inter_combination, intra_combination)

data$tissue <- factor(data$tissue, levels = c('brain', 'skin', 'blood'))
data$position <- factor(data$position, levels = c('intra', 'inter'))

write.table(data, 
            file = 'results/tables/te_host_gene_correlation_data.csv',
            sep = ',', quote = F, row.names = F, col.names = T)

