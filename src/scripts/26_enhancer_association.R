# The purpose of that script is identify CAGE-Peaks that intersect with TEs and enhancers. The output
# of the script is a table that gives enhancer-TEs and their associated genes regulated in the same
# direction.

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}


# Get all intergenic TEs to store them in a bed to make a bedtools closest 
# to get the neighbored genes.
#te.annotation <- te.annotation %>% 
#    filter(chromosome %in% chr_of_interest) #%>% 

###############
# load tracks #
###############

geneHancer <- readGeneric('/misc/paras/data/rschwarz/common_data/mm10/geneHancer/geneHancer.bed',
                          strand = 6,
                          meta.cols = list(names = 4))

cage_peaks <- list()

cage_peaks[['blood']] <- readGeneric('results/cage_RS/blood_gclipped/raw_peaks/blood.cage_peaks.bed',
                                     strand = 6, 
                                     meta.cols = list(names = 4))

cage_peaks[['brain']] <- readGeneric('results/cage_RS/brain_downsampled_gclipped/raw_peaks/brain.cage_peaks.bed',
                                     strand = 6, 
                                     meta.cols = list(names = 4))

cage_peaks[['skin']] <- readGeneric('results/cage_RS/skinII_gclipped/raw_peaks/skin.cage_peaks.bed',
                                    strand = 6, 
                                    meta.cols = list(names = 4))

############################################
# Collect CAGE-TE-Enhancer-Gene connection #
############################################

te_CAGEs_enh_Gens <- list()

for (tissue in tissues) {
    
    te_cages <-  intersectGranger(cage_peaks[[tissue]], 
                                  teRanges, 'all') %>% 
        dplyr::select(query.names, subject.te.id) %>% 
        dplyr::rename(peak_id = query.names,
                      te_id = subject.te.id)
    
    tmp_te_CAGEs <- intersectGranger(cage_peaks[[tissue]], 
                                           teRanges, 'query')
    
    te_enhancer_Cages <- intersectGranger(tmp_te_CAGEs, 
                                           geneHancer, 'all') %>% 
        dplyr::select(query.names, subject.names) %>% 
        tidyr::separate(col = "subject.names", into = c("external_gene_name", "enhancer_id"),
                        sep = '/') %>% 
        dplyr::rename(peak_id = query.names)
    
    te_CAGEs_enh_Gens[[tissue]] <- merge(te_cages, te_enhancer_Cages, by = 'peak_id')
    
    
}

# load DESeq results and merge it to the tables
gene_deseq <- fread('results/tables/02_deseq_results_gene.csv')
te_instances_deseq <- fread('results/tables/02_deseq_results_te_instances.csv')
cage_deseq <- fread('results/tables/07_deseq_cage_all.csv')

gene_deseq <- gene_deseq %>% 
    dplyr::select(ensembl_gene_id, external_gene_name, log2FoldChange, padj, tissue) %>% 
    dplyr::rename(gene_log2FC = log2FoldChange,
                  gene_padj = padj)

te_instances_deseq <- te_instances_deseq %>% 
    dplyr::select(te_id, log2FoldChange, tissue, padj, position) %>% 
    dplyr::rename(te_log2FC = log2FoldChange,
                  te_padj = padj)

cage_deseq <- cage_deseq %>% 
    dplyr::select(peak_id, log2FoldChange, padj, tissue) %>% 
    dplyr::rename(cage_log2FC = log2FoldChange,
                  cage_padj = padj)

### Brain

te_CAGEs_enh_genes_brain <- merge(te_CAGEs_enh_Gens$brain, 
                                  gene_deseq %>% filter(tissue == 'brain'),
                                  by = 'external_gene_name') 

te_CAGEs_enh_genes_brain <- merge(te_CAGEs_enh_genes_brain,
                                  te_instances_deseq %>% filter(tissue == 'brain'),
                                  by = c('te_id', 'tissue'))

te_CAGEs_enh_genes_brain <- merge(te_CAGEs_enh_genes_brain,
                                  cage_deseq %>% filter(tissue == 'brain'),
                                  by = c('peak_id', 'tissue'))
### blood

te_CAGEs_enh_genes_blood <- merge(te_CAGEs_enh_Gens$blood, 
                                  gene_deseq %>% filter(tissue == 'blood'),
                                  by = 'external_gene_name') 

te_CAGEs_enh_genes_blood <- merge(te_CAGEs_enh_genes_blood,
                                  te_instances_deseq %>% filter(tissue == 'blood'),
                                  by = c('te_id', 'tissue'))

te_CAGEs_enh_genes_blood <- merge(te_CAGEs_enh_genes_blood,
                                  cage_deseq %>% filter(tissue == 'blood'),
                                  by = c('peak_id', 'tissue'))
### skin

te_CAGEs_enh_genes_skin <- merge(te_CAGEs_enh_Gens$skin, 
                                  gene_deseq %>% filter(tissue == 'skin'),
                                  by = 'external_gene_name') 

te_CAGEs_enh_genes_skin <- merge(te_CAGEs_enh_genes_skin,
                                  te_instances_deseq %>% filter(tissue == 'skin'),
                                  by = c('te_id', 'tissue'))

te_CAGEs_enh_genes_skin <- merge(te_CAGEs_enh_genes_skin,
                                  cage_deseq %>% filter(tissue == 'skin'),
                                  by = c('peak_id', 'tissue'))

# bind tissue specific tables to one comprehensive table
data_full <- rbind(te_CAGEs_enh_genes_blood, rbind(te_CAGEs_enh_genes_skin, te_CAGEs_enh_genes_brain))

# Filter those where CAGE and Genes are equally regulated
for (p in c(0.05, 0.1, 1)) {
    
    data_table_up <- data_full %>% filter(gene_log2FC > 0, cage_log2FC > 0, cage_padj <= p)
    data_table_down <- data_full %>% filter(gene_log2FC < 0, cage_log2FC < 0, cage_padj <= p)
    
    data_table_cage <- distinct(rbind(data_table_up, data_table_down))
    
    data_table_cage <- splitTEID(data_table_cage, 'te_id')
    
    write.table(data_table_cage, file = paste0(table_dir, "26_TE_enhancers_adjp_", p, ".csv"),
                row.names = FALSE, quote = FALSE, sep = ',')

}

# ============================ Filtered for DEGs =================================================
#
# Filter those where TE_CAGEs and Genes are regulated in the same direction and the gene is
# differentially expressed.

for (p in c(0.05, 0.1, 1)) {
    
    data_table_gene_up <- data_full %>% 
        filter(gene_log2FC > 0, cage_log2FC > 0, gene_padj <= p)
    data_table_gene_down <- data_full %>% 
        filter(gene_log2FC < 0, cage_log2FC < 0, gene_padj <= p)
    
    data_table_gene <- unique(rbind(data_table_gene_up, data_table_gene_down))
    
    data_table_gene <- splitTEID(data_table_gene, 'te_id')
    
    write.table(data_table_gene, file = paste0(table_dir, "26_TE_enhancers_gene_adjp_", p, ".csv"),
                row.names = FALSE, quote = FALSE, sep = ',')
}
