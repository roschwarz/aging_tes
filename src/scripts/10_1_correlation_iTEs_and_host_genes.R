# The purpose of that script is identify independently expressed intronic TEs (TEs that intersect with a 
# CAGE-peak) and check the correlation of the l2FC of these TEs and their host gene.

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}


tissue_names = c('blood_gclipped',
                 'brain_downsampled_gclipped',
                 'skinII_gclipped')

# =================================== TE-Annotations ===============================================

# Get all intergenic TEs to store them in a bed to make a bedtools closest 
# to get the neighbored genes.
te.annotation <- te.annotation %>% 
    filter(chromosome %in% chr_of_interest) #%>% 

intronic_TEs_bed <- te.annotation %>% 
        filter(position == 'intronic') %>% 
        dplyr::select(chromosome, start, end, te_id, Kimura, strand) %>% 
        mutate(Kimura = 1)
    
intronic_TEs_gr <- GRanges(seqnames = intronic_TEs_bed$chromosome,
                          ranges = IRanges(as.numeric(intronic_TEs_bed$start),
                                           end = as.numeric(intronic_TEs_bed$end),
                                           names = intronic_TEs_bed$te_id),
                          strand = intronic_TEs_bed$strand,
                          te_id = intronic_TEs_bed$te_id)

# Contains double TE-ids when multiple genes have the same distance
te_gene_distance <- read.csv(paste0(table_dir,
                                    "intragenic_te_gene_distance.csv"))

intergenic_TEs <- te.annotation %>% 
    filter(te_id %in% te_gene_distance$te_id) %>% 
    dplyr::select(-c(ensembl_gene_id, external_gene_name))

gene_anno <- data.frame(geneRanges) %>% 
    dplyr::select(ensembl_gene_id, external_gene_name)

gene_anno <- merge(te_gene_distance, gene_anno, by = 'ensembl_gene_id')

intergenic_TEs <- merge(intergenic_TEs, gene_anno, by = 'te_id' )

intergenic_TEs_gr <- GRanges(seqnames = intergenic_TEs$chromosome,
                          ranges = IRanges(as.numeric(intergenic_TEs$start),
                                           end = as.numeric(intergenic_TEs$end),
                                           names = intergenic_TEs$te_id),
                          strand = intergenic_TEs$strand,
                          te_id = intergenic_TEs$te_id)

# ======================================== CAGE-Peaks================================================
    
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


independent_intronic_tes <- list()

for (tissue in tissues) {
    
    independent_intronic_tes[[tissue]] <- intersectFeatures(intronic_TEs_gr, 
                                                            cage_peaks[[tissue]], 'te_id')
    
}
 

gene_deseq <- fread('results/tables/02_deseq_results_gene.csv')
te_instances_deseq <- fread('results/tables/02_deseq_results_te_instances.csv')

introList <- list()

tes_brain <- te_instances_deseq %>% 
        filter(tissue == 'brain',
               te_id %in% independent_intronic_tes$brain#,
               #baseMean > 5,
               #abs(log2FoldChange) > 0.1
               ) %>% 
        dplyr::select(te_id, log2FoldChange, ensembl_gene_id) %>% 
        dplyr::rename(te_log2FoldChange = log2FoldChange)
    
gene_brain <- gene_deseq %>% 
       filter(tissue == 'brain', ensembl_gene_id %in% tes_brain$ensembl_gene_id#,
              #baseMean > 5,
              #abs(log2FoldChange) > 0.1
              ) %>% 
       dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
       dplyr::rename(gene_log2FoldChange = log2FoldChange) 
   
merge_brain <- merge(gene_brain, tes_brain, by = 'ensembl_gene_id')

tes_skin <- te_instances_deseq %>% 
        filter(tissue == 'skin',
               te_id %in% independent_intronic_tes$skin,
               #baseMean > 5,
               #abs(log2FoldChange) > 0.1
               ) %>% 
        dplyr::select(te_id, log2FoldChange, ensembl_gene_id) %>% 
        dplyr::rename(te_log2FoldChange = log2FoldChange)
    
gene_skin <- gene_deseq %>% 
       filter(tissue == 'skin', ensembl_gene_id %in% tes_skin$ensembl_gene_id,
              #baseMean > 5,
              #abs(log2FoldChange) > 0.1
              ) %>% 
       dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
       dplyr::rename(gene_log2FoldChange = log2FoldChange) 
   
merge_skin <- merge(gene_skin, tes_skin, by = 'ensembl_gene_id')


tes_blood <- te_instances_deseq %>% 
        filter(tissue == 'blood',
               te_id %in% independent_intronic_tes$blood,
               #baseMean > 5,
               #abs(log2FoldChange) > 0.1
               ) %>% 
        dplyr::select(te_id, log2FoldChange, ensembl_gene_id) %>% 
        dplyr::rename(te_log2FoldChange = log2FoldChange)
    
gene_blood <- gene_deseq %>% 
       filter(tissue == 'blood', ensembl_gene_id %in% tes_blood$ensembl_gene_id,
              #baseMean > 5,
              #abs(log2FoldChange) > 0.1
              ) %>% 
       dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
       dplyr::rename(gene_log2FoldChange = log2FoldChange) 
   
merge_blood <- merge(gene_blood, tes_blood, by = 'ensembl_gene_id')

intro_combination <- rbind(merge_brain, rbind(merge_skin, merge_blood))

intro_combination$tissue <- factor(intro_combination$tissue, level = tissues)

pl <- ggplot(intro_combination, aes(gene_log2FoldChange, te_log2FoldChange)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(method = 'lm', color = '#028763') +
    stat_quadrant_counts(xintercept = 0, yintercept = 0) +
    coord_cartesian() +
    stat_cor(method = "spearman", label.y = 3, label.x = -2) +
    #stat_poly_eq(use_label= c("adj.R2", "p"),label.y = 0.7, label.x = 0.1) +
    facet_grid(.~tissue) +
    theme_rob() +
    theme(text = element_text(family = 'arial'))


pl <- color_strips(pl, 
                   bg_cols = c(tissue.color[2:4]), 
                   text_cols = c( "#ffffff", "#000000","#ffffff"))

pl <- grid.draw(pl)

show(pl)

#########################################################################
##### Intronic TE_CAGES-Gene-correlation
###
##
#
#########################################################################

te_cages <- list()

gene_deseq <- fread('results/tables/02_deseq_results_gene.csv')
cage_deseq <- fread('results/tables/07_deseq_cage_all.csv')
te_instances_deseq <- fread('results/tables/02_deseq_results_te_instances.csv')

for (tissue in tissues) {
    
    te_cages[[tissue]] <-  intersectGranger(cage_peaks[[tissue]], 
                                  intronic_TEs_gr, 'all') %>% 
        dplyr::select(query.names, subject.te_id) %>% 
        dplyr::rename(peak_id = query.names,
                      te_id = subject.te_id) %>% 
        dplyr::mutate(tissue = tissue)
    
    te_gene_asso <- te_instances_deseq %>% 
        filter(tissue == tissue) %>% 
        dplyr::select(te_id, ensembl_gene_id)
    
    te_cages[[tissue]] <- merge(te_cages[[tissue]], 
                                te_gene_asso, 
                                by = 'te_id')
    
    cage_deseq_temp <- cage_deseq %>% 
        filter(tissue == tissue) %>% 
        dplyr::select('peak_id', 'log2FoldChange', 'padj') %>% 
        dplyr::rename(cage_log2FC = log2FoldChange,
                      cage_padj = padj)
    
    gene_deseq_temp <- gene_deseq %>% 
        filter(tissue == tissue) %>% 
        dplyr::select('ensembl_gene_id', 'log2FoldChange', 'padj', 'external_gene_name') %>% 
        dplyr::rename(gene_log2FC = log2FoldChange,
                      gene_padj = padj)
    
    te_cages[[tissue]] <- merge(te_cages[[tissue]], 
                                cage_deseq_temp, 
                                by = 'peak_id')
    
    te_cages[[tissue]] <- merge(te_cages[[tissue]], 
                                gene_deseq_temp, 
                                by = 'ensembl_gene_id')
    
}

te_cages_comp <- do.call('rbind', te_cages)

te_cages_comp$tissue <- factor(te_cages_comp$tissue, level = tissues)

pl <- ggplot(te_cages_comp, aes(gene_log2FC, cage_log2FC)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(method = 'lm', color = '#028763') +
    stat_quadrant_counts(xintercept = 0, yintercept = 0) +
    coord_cartesian() +
    stat_cor(method = "spearman", label.y = 6, label.x = -3) +
    #stat_poly_eq(label.y = 0.7, label.x = 0.1) +
    facet_grid(.~tissue) +
    theme_rob() +
    theme(text = element_text(family = 'arial'))


pl <- color_strips(pl, 
                   bg_cols = c(tissue.color[2:4]), 
                   text_cols = c( "#ffffff", "#000000","#ffffff"))

pl <- grid.draw(pl)

show(pl)

#########################################################################
##### Intergenic TE_CAGES-Gene-correlation
###
##
#
#########################################################################

independent_intergenic_tes <- list()

for (tissue in tissues) {
    
    independent_intergenic_tes[[tissue]] <- intersectFeatures(intergenic_TEs_gr, 
                                                            cage_peaks[[tissue]], 'te_id')
    
}


tes_brain <- te_instances_deseq %>% 
    filter(tissue == 'brain',
           te_id %in% independent_intergenic_tes$brain#,
           #baseMean > 5,
           #abs(log2FoldChange) > 0.1
    ) %>% 
    dplyr::select(te_id, log2FoldChange) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange)

tes_brain <- merge(tes_brain, te_gene_distance, by='te_id')

gene_brain <- gene_deseq %>% 
    filter(tissue == 'brain', ensembl_gene_id %in% tes_brain$ensembl_gene_id#,
           #baseMean > 5,
           #abs(log2FoldChange) > 0.1
    ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange) 

merge_brain <- merge(gene_brain, tes_brain, by = 'ensembl_gene_id')

tes_skin <- te_instances_deseq %>% 
    filter(tissue == 'skin',
           te_id %in% independent_intergenic_tes$skin,
           #baseMean > 5,
           #abs(log2FoldChange) > 0.1
    ) %>% 
    dplyr::select(te_id, log2FoldChange) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange)

tes_skin <- merge(tes_skin, te_gene_distance, by='te_id')

gene_skin <- gene_deseq %>% 
    filter(tissue == 'skin', ensembl_gene_id %in% tes_skin$ensembl_gene_id,
           #baseMean > 5,
           #abs(log2FoldChange) > 0.1
    ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange) 

merge_skin <- merge(gene_skin, tes_skin, by = 'ensembl_gene_id')


tes_blood <- te_instances_deseq %>% 
    filter(tissue == 'blood',
           te_id %in% independent_intergenic_tes$blood,
           #baseMean > 5,
           #abs(log2FoldChange) > 0.1
    ) %>% 
    dplyr::select(te_id, log2FoldChange) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange)

tes_blood <- merge(tes_blood, te_gene_distance, by='te_id')

gene_blood <- gene_deseq %>% 
    filter(tissue == 'blood', ensembl_gene_id %in% tes_blood$ensembl_gene_id,
           #baseMean > 5,
           #abs(log2FoldChange) > 0.1
    ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange) 

merge_blood <- merge(gene_blood, tes_blood, by = 'ensembl_gene_id')

intergenic_combination <- rbind(merge_brain, rbind(merge_skin, merge_blood))

intergenic_combination$tissue <- factor(intergenic_combination$tissue, level = tissues)

pl <- ggplot(intergenic_combination, aes(gene_log2FoldChange, te_log2FoldChange)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(method = 'lm', color = '#028763') +
    stat_quadrant_counts(xintercept = 0, yintercept = 0) +
    coord_cartesian() +
    stat_cor(method = "spearman", label.y = 3, label.x = -2) +
    #stat_poly_eq(use_label= c("adj.R2", "p"),label.y = 0.7, label.x = 0.1) +
    facet_grid(.~tissue) +
    theme_rob() +
    theme(text = element_text(family = 'arial'))


pl <- color_strips(pl, 
                   bg_cols = c(tissue.color), 
                   text_cols = c( "#ffffff", "#000000","#ffffff"))

pl <- grid.draw(pl)

show(pl)

