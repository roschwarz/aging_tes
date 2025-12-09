# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()
aging_tes::load_annotations()

load_gene_ranges()

indie_te_islands <- sapply(names(indie_te_island_bed), simplify = F, function(x){
    
    indie_df <- readGeneric(indie_te_island_bed[[x]], strand = 6, 
                           meta.cols = list(names = 4))
    
    inde_te_islands_genes_df <- intersectGranger(geneRanges, indie_df, 'all') %>% 
        dplyr::select(query.ensembl_gene_id, query.external_gene_name, subject.names) %>% 
        dplyr::rename(ensembl_gene_id = query.ensembl_gene_id, 
                      external_gene_name = query.external_gene_name, 
                      te_island_id = subject.names) 
    
    return(inde_te_islands_genes_df)
    
})

# ================================= Load TE island rna =========================

te_island_deseq_results <- read.csv(paste0(table_dir, "02_deseq_results_te_island.csv"))

# indie_te_island_rna <- sapply(tissues, simplify = F, function(x){
#     
#     df <- te_island_deseq_results %>% filter(tissue == x,
#                                    te_island_id %in% indie_te_islands[[x]][['te_island_id']]) %>% 
#         dplyr::select(te_island_id, log2FoldChange, padj, baseMean, pvalue)
#     
#     names(df) <- c('te_island_id', 'rna_L2FC', 'rna_padj', 'te_island_baseMean', 'te_island_pvalue')
#     
#     return(df)
# })

indie_te_island_rna <- sapply(tissues, simplify = F, function(x){
    
    df <- te_island_deseq_results %>% filter(tissue == x)
    
    df <- merge(df, indie_te_islands[[x]]) %>% 
                                   #te_island_id %in% indie_te_islands[[x]][['te_island_id']]) %>% 
        dplyr::select(te_island_id, log2FoldChange, padj, baseMean, pvalue, ensembl_gene_id)
    
    names(df) <- c('te_island_id', 'rna_L2FC', 'rna_padj', 'te_island_baseMean', 'te_island_pvalue', 'ensembl_gene_id')
    
    return(df)
})

# ================================= Load TE island cage ========================
#
# used to identify individually expressed TE islands in the rna-seq data
# cage_te_island <- read.csv(paste0(table_dir, "deseq_cage_te_island.csv"))
# 
# indie_te_island_cages <- sapply(tissues, simplify = F, function(x){
#     
#     df <- cage_te_island %>% 
#         filter(tissue == x, te_island.names %in% indie_te_islands[[x]][['te_island_id']]) %>% 
#         dplyr::select(peak_id, te_island.names, ensembl_gene_id, log2FoldChange, padj, baseMean, pvalue)
#     
#     names(df) <- c('peak_id', 'te_island_id', 'ensembl_gene_id', 'cage_L2FC', 'cage_padj', 'cage_baseMean', 'cage_pvalue')
#     
#     return(df)
# })

# ================================= Load Gene data =============================

rna_gene <- read.csv(paste0(table_dir, "02_deseq_results_genes.csv"))

genes_with_indie_te_islands <- sapply(tissues, simplify = FALSE, function(x){
    
    df <- rna_gene %>% 
        filter(tissue == x, ensembl_gene_id %in% indie_te_islands[[x]][['ensembl_gene_id']]) %>% 
        dplyr::select(ensembl_gene_id, log2FoldChange, padj, external_gene_name, baseMean, pvalue)
    
    names(df) <- c('ensembl_gene_id', 'gene_L2FC', 'gene_padj', 'external_gene_name', 'gene_baseMean', 'gene_pvalue')
    
    return(df)
})


# ================================= Merge Tables ===============================

gene_te_island_association <- sapply(tissues, simplify = F, function(x){
    
    df <- indie_te_island_rna[[x]]
    #df <- merge(indie_te_island_cages[[x]], indie_te_island_rna[[x]], by = 'te_island_id', all = TRUE)
    df <- merge(df, genes_with_indie_te_islands[[x]], by = 'ensembl_gene_id')
    
    return(df)
    
})

gene_te_island_association_merged <- do.call('rbind', sapply(tissues, simplify = F, function(x){
    gene_te_island_association[[x]] %>% 
        mutate(tissue = x) %>% 
        mutate(rna_quadrant = ifelse(gene_L2FC > 0 & rna_L2FC > 0, 1,
                                     ifelse(gene_L2FC > 0 & rna_L2FC < 0, 2,
                                            ifelse(gene_L2FC < 0 & rna_L2FC < 0, 3, 4))))
}))

write.table(gene_te_island_association_merged,
            file = paste0(table_dir, 'indie_te_island_host_gene_correlation.csv'),
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

# ================================= Brain specific ===============================
#
# The protocadherine cluster popped up as an interesting locus, therefore the te 
# islands within this cluster should be considered as well. That requires an updated
# data frame.

load_te_island_annotation()

pcdhb <- rna_gene %>% 
    filter(tissue == 'brain',
           grepl(pattern = '^Pcdhb', x = external_gene_name))

pcdhb_bed <- transGrange(geneRanges) %>% 
    filter(ensembl_gene_id %in% pcdhb$ensembl_gene_id) %>% 
    dplyr::select(seqnames, start ,end, external_gene_name, width, strand) %>% 
    arrange(start)


write.table(pcdhb_bed,
            file = paste0(table_dir, 'protocadherine_in_pcdhb_cluster.bed'),
            col.names = F, row.names = F, quote = F, 
            sep = '\t')

# handpicked in the genome browser
cluster_of_interest <- c('TE_Cluster_728951', 
                         'TE_Cluster_728962', 
                         'TE_Cluster_728983', #
                         'TE_Cluster_728986', #
                         'TE_Cluster_728991') #

te_island_bed <- transGrange(te_island_5primeRanges) %>% 
    filter(names %in% cluster_of_interest) %>% 
    dplyr::select(seqnames, start ,end, names, width, strand)


write.table(te_island_bed,
            file = paste0(table_dir, 'te_islands_in_pcdhb_cluster.bed'),
            col.names = F, row.names = F, quote = F, 
            sep = '\t')


df <- system(paste0('bedtools closest -a ', table_dir, 'te_islands_in_pcdhb_cluster.bed -b ', 
                    table_dir, 'protocadherine_in_pcdhb_cluster.bed -s'), wait = TRUE, intern = TRUE)


pcdhb_island_pairs <- read.table(text = df) %>% 
    dplyr::select(V4, V10) %>% 
    dplyr::rename(te_island_id = V4,
                  external_gene_name = V10)


pcdhb_island_pairs <- merge(pcdhb_island_pairs, 
                            te_island_deseq_results %>% filter(tissue == "brain"), 
                            by.x = 'te_island_id', by.y = 'te_island_id') %>% 
    dplyr::select(te_island_id, external_gene_name, log2FoldChange, padj) %>% 
    dplyr::rename(rna_L2FC = log2FoldChange,
                  rna_padj = padj)


pcdhb_island_pairs <- merge(pcdhb_island_pairs,
                            rna_gene %>% filter(tissue == 'brain'),
                            by = 'external_gene_name') %>% 
    mutate(peak_id = NA,
           cage_L2FC = NA,
           cage_padj = NA) %>% 
    dplyr::select(ensembl_gene_id, te_island_id, 
                  peak_id, cage_L2FC, cage_padj,
                  rna_L2FC, rna_padj, log2FoldChange, padj, external_gene_name) %>% 
    dplyr::rename(gene_L2FC = log2FoldChange,
                  gene_padj = padj)

brain_data_pcdhb <- rbind(brain_data %>% 
                              dplyr::select(ensembl_gene_id, te_island_id, 
                                            peak_id, cage_L2FC, cage_padj,
                                            rna_L2FC, rna_padj, gene_L2FC,
                                            gene_padj, external_gene_name),
                          pcdhb_island_pairs)




