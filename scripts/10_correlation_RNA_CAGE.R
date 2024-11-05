# Load Environment
if( !exists("ENVIRONMENT_LOADED") ){
    
    source('./01_load_environment.R')
    
} else if( !ENVIRONMENT_LOADED ){
    
    source('./01_load_environment.R')
    
}

library(extrafont)

font_import()

loadfonts(device = "postscript")

# ---------------------------- Prepare data ------------------------------------
#
# An autonomous TE region is a region that contains at least one CAGE-Peak and
# TE instance that is considered as expressed and was already defined in 
# TE.region.R and stored as a bed file. Currently, the coordinates of these 
# autonomous TE region are stored in ../data/densityPlot/, but that will changed 
# when the scripts are updated.
#
# Load coordinates of autonomously expressed TE regions as a Grange object
# These are all TE regions with at least one CAGE-Peak and RNA-Seq signal
# Subsequently, the autonomous TEs are intersected with genes to identify the
# host genes of autonomous TEs that are located in genes.
#
# Afterwards, the deseq results of CAGE and RNA for TE regions are loaded and
# filtered for the set of autonomous TEs that are located in genes. The same
# is done for the RNA deseq results of genes, whereas this table is filtered 
# for the host genes.
#
# Finally, the filtered tables are merged and can be used for the down stream
# analysis. The merged table is stored as 10_auto_TE_Host_correlation.csv in
# the table directory. The rna_quadrant column is the information in which 
# quadrant the data point is located considering only the correlation between
# the rna signal of the TE region and the host gene.


# =================== Load annotation & identify host genes ====================

auto.TE.cluster.files <- list(brain = "data/shared/brain_downsampled_independent_TE_regions.bed",
                              blood = "data/shared/blood_independent_TE_regions.bed",
                              skin = "data/shared/skinII_independent_TE_regions.bed"
)

auto_TE_regions <- sapply(names(auto.TE.cluster.files), simplify = F, function(x){
    
    auto.df <- readGeneric(auto.TE.cluster.files[[x]], strand = 6, 
                           meta.cols = list(names = 4))
    
    auto.gene.df <- intersectGranger(geneRanges, auto.df, 'all') %>% 
        dplyr::select(query.ensembl_gene_id, query.external_gene_name, subject.names) %>% 
        dplyr::rename(ensembl_gene_id = query.ensembl_gene_id, 
                      external_gene_name = query.external_gene_name, 
                      te_region_id = subject.names) 
    
    return(auto.gene.df)
    
})

# ================================= Load TE region rna =========================

te.region.rna <- read.csv(paste0(table_dir, "02_deseq_results_te_region.csv"))

auto_TE_regions_rna <- sapply(tissues, simplify = F, function(x){
    
    df <- te.region.rna %>% filter(tissue == x,
                                   te_region_id %in% auto_TE_regions[[x]][['te_region_id']]) %>% 
        dplyr::select(te_region_id, log2FoldChange, padj, baseMean, pvalue)
    
    names(df) <- c('te_region_id', 'rna_L2FC', 'rna_padj', 'te_region_baseMean', 'te_region_pvalue')
    
    return(df)
})

# ================================= Load TE region cage ========================

te.region.cage <- read.csv(paste0(table_dir, "07_deseq_cage_te_region.csv"))

auto_TE_regions_cage <- sapply(tissues, simplify = F, function(x){
    
    df <- te.region.cage %>% 
        filter(tissue == x, te_region.names %in% auto_TE_regions[[x]][['te_region_id']]) %>% 
        dplyr::select(peak_id, te_region.names, ensembl_gene_id, log2FoldChange, padj, baseMean, pvalue)
    
    names(df) <- c('peak_id', 'te_region_id', 'ensembl_gene_id', 'cage_L2FC', 'cage_padj', 'cage_baseMean', 'cage_pvalue')
    
    return(df)
})

# ================================= Load Gene data =============================


gene.rna <- read.csv(paste0(table_dir, "02_deseq_results_gene.csv"))

genes_with_auto_TE_regions <- sapply(tissues, simplify = FALSE, function(x){
   
    df <- gene.rna %>% 
        filter(tissue == x, ensembl_gene_id %in% auto_TE_regions[[x]][['ensembl_gene_id']]) %>% 
        dplyr::select(ensembl_gene_id, log2FoldChange, padj, external_gene_name, baseMean, pvalue)
    
    names(df) <- c('ensembl_gene_id', 'gene_L2FC', 'gene_padj', 'external_gene_name', 'gene_baseMean', 'gene_pvalue')
    
    return(df)
})

# ================================= Merge Tables ===============================

combined.df <- sapply(tissues, simplify = F, function(x){
    
    df <- merge(auto_TE_regions_cage[[x]], auto_TE_regions_rna[[x]], by = 'te_region_id')
    df <- merge(df, genes_with_auto_TE_regions[[x]], by = 'ensembl_gene_id')
    
    return(df)
    
})

combined.df.merged <- do.call('rbind', sapply(tissues, simplify = F, function(x){
    combined.df[[x]] %>% 
        mutate(tissue = x) %>% 
        mutate(rna_quadrant = ifelse(gene_L2FC > 0 & rna_L2FC > 0, 1,
                                     ifelse(gene_L2FC > 0 & rna_L2FC < 0, 2,
                                            ifelse(gene_L2FC < 0 & rna_L2FC < 0, 3, 4))))
}))

write.table(combined.df.merged,
            file = paste0(table_dir, '10_auto_TE_Host_correlation.csv'),
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

# combined.resis <- combined.df.merged %>% 
#     group_by(tissue) %>% 
#     filter(!duplicated(te_region_id))
# 
# 
# write.table(combined.resis,
#             file = '../../submission/Resis/figure4/figure_4b.csv',
#             col.names = T, 
#             row.names = F,
#             quote = F,
#             sep = ',')


# ------------------------------ Correlation Plots -----------------------------
#
# The tables may contain multiple entries per TE region, as there is a 
# possibility that more than one CAGE-peak may exists per TE region. 
# Consequently, a filtered table is used for the log2FC correlation plots.
# 
# ============================== Brain =========================================

brain_data <- combined.df$brain %>% 
    filter(!duplicated(te_region_id))

faded <- brain_data %>% filter(abs(rna_L2FC) < 0.1)
focus <- brain_data %>% filter(abs(rna_L2FC) > 0.1)


ggplot(brain_data, aes(rna_L2FC, gene_L2FC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point_rast(data = faded, alpha = 0.05, color = '#192E37') +
    geom_point(data = focus, color = '#192E37') +
    geom_smooth(method = 'lm', color = '#028763') +
    xlim(c(-2,2)) +
    ylim(c(-2,2)) +
    geom_text_repel(data = focus,
                     aes(label = external_gene_name), 
                     size = 10*0.36,
                     box.padding = 0.7,
                     #label.padding = 0.1,
                     max.overlaps = 50) +
    stat_quadrant_counts(xintercept = 0, yintercept = 0) +
    labs(x = expression(paste(log[2], "(fold ", change[TE_region], ")")),
         y = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    theme_rob() +
    theme(text = element_text(family = 'arial'))

summary(
    lm(rna_L2FC~gene_L2FC, brain_data %>% filter(rna_L2FC > 0.1, gene_L2FC > 0.1))
    )

cor.test(brain_data$gene_L2FC, brain_data$rna_L2FC, method = 'spearman')


ggsave(
    filename = paste0(figure_dir, '10_brain_correlation_host_gene_autonomous_TE_regions.svg'),
    plot = last_plot(),
    width = 20,
    height = 20,
    units = "cm",
    dpi = 300)

# ########################## Brain with protocaderhins #########################

pcdhb <- gene.rna %>% 
    filter(tissue == 'brain',
           grepl(pattern = '^Pcdhb', x = external_gene_name))


pcdhb.bed <- transGrange(geneRanges) %>% 
    filter(ensembl_gene_id %in% pcdhb$ensembl_gene_id) %>% 
    dplyr::select(seqnames, start ,end, external_gene_name, width, strand) %>% 
    arrange(start)


write.table(pcdhb.bed,
            file = paste0(table_dir, '10_pcdhb.bed'),
            col.names = F, row.names = F, quote = F, 
            sep = '\t')

# handpicked in the genome browser
cluster.of.interest <- c('TE_Cluster_728951', 
                         'TE_Cluster_728962', 
                         'TE_Cluster_728983', #
                         'TE_Cluster_728986', #
                         'TE_Cluster_728991') #

te.region.bed <- transGrange(teregionRanges) %>% 
    filter(names %in% cluster.of.interest) %>% 
    dplyr::select(seqnames, start ,end, names, width, strand)


write.table(te.region.bed,
            file = paste0(table_dir, '10_pcdhb.te.regions.bed'),
            col.names = F, row.names = F, quote = F, 
            sep = '\t')


df <- system(paste0('bedtools closest -a ', table_dir, '10_pcdhb.te.regions.bed -b ', 
             table_dir, '10_pcdhb.bed -s'), wait = TRUE, intern = TRUE)

pcdhb.region.pairs <- read.table(text = df) %>% 
    dplyr::select(V4, V10) %>% 
    dplyr::rename(te_region_id = V4,
                  external_gene_name = V10)


pcdhb.region.pairs <- merge(pcdhb.region.pairs, 
                            te.region.rna %>% filter(tissue == 'brain'), 
                            by.x = 'te_region_id', by.y = 'te_region_id') %>% 
    dplyr::select(te_region_id, external_gene_name, log2FoldChange, padj) %>% 
    dplyr::rename(rna_L2FC = log2FoldChange,
                  rna_padj = padj)


pcdhb.region.pairs <- merge(pcdhb.region.pairs,
                            gene.rna %>% filter(tissue == 'brain'),
                            by = 'external_gene_name') %>% 
    mutate(peak_id = NA,
           cage_L2FC = NA,
           cage_padj = NA) %>% 
    dplyr::select(ensembl_gene_id, te_region_id, 
                  peak_id, cage_L2FC, cage_padj,
                  rna_L2FC, rna_padj, log2FoldChange, padj, external_gene_name) %>% 
    dplyr::rename(gene_L2FC = log2FoldChange,
                  gene_padj = padj)

brain_data.pcdhb <- rbind(brain_data %>% 
                              dplyr::select(ensembl_gene_id, te_region_id, 
                                            peak_id, cage_L2FC, cage_padj,
                                            rna_L2FC, rna_padj, gene_L2FC,
                                            gene_padj, external_gene_name),
                          pcdhb.region.pairs)

faded <- brain_data.pcdhb %>% filter(abs(rna_L2FC) < 0.1)
focus <- brain_data.pcdhb %>% filter(abs(rna_L2FC) > 0.1)


brain_data.pcdhb$tissue <- 'brain'

corr_brain <- ggplot(brain_data.pcdhb, aes(rna_L2FC, gene_L2FC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point_rast(data = faded, 
                    alpha = 0.05, 
                    color = '#192E37',
                    size = 0.3) +
    geom_point(data = focus, 
               color = '#192E37',
               size = 0.3) +
    geom_smooth(method = 'lm', 
                color = '#028763',
                linewidth = 0.3) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'Rho',
             p.accuracy = .Machine$double.xmin,
             size = 6/.pt,
             label.sep = '\n',
             label.x = -2, 
             label.y = 1.5,
             family = 'Arial',) +
    xlim(c(-2,2)) +
    ylim(c(-2,2)) +
    geom_text_repel(data = focus,
                    aes(label = external_gene_name), 
                    family = 'Arial',
                    size = 6/.pt,
                    box.padding = 0.25,
                    #label.padding = 0.3,
                    max.overlaps = 80) +
    stat_quadrant_counts(xintercept = 0, 
                         yintercept = 0,  
                         size = 6/.pt,
                         family = 'Arial',) +
    facet_grid(.~tissue) +
    labs(x = expression(paste(log[2], "(fold ", change[TE_region], ")")),
         y = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          #text = element_text(family = 'arial'),
          strip.background = element_rect(fill = tissue.color['brain']),
          strip.text = element_text(size = 8, colour = 'white'))



# =============================== Skin =========================================

skin_data <- combined.df$skin %>% 
    filter(!duplicated(te_region_id))

skin_faded <- skin_data %>% filter(abs(rna_L2FC) < 0.3)
skin_focus <- skin_data %>% filter(abs(rna_L2FC) > 0.3)

skin_data$tissue <- 'skin'

corr_skin <- ggplot(skin_data, aes(rna_L2FC, gene_L2FC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point_rast(data = skin_faded, 
                    alpha = 0.05,
                    size = 0.3,
                    color = '#192E37') +
    geom_point(data = skin_focus,
               color = '#192E37',
               size = 0.3) +
    geom_smooth(method = 'lm', 
                color = '#028763',
                linewidth = 0.3) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'Rho',
             size = 6/.pt,
             family = "Arial",
             p.accuracy = .Machine$double.xmin,
             label.sep = '\n',
             label.x = -4, 
             label.y = 3) +
    xlim(c(-4,4)) +
    ylim(c(-4,4)) +
    geom_text_repel(data = skin_focus %>% filter(abs(rna_L2FC) > 1.5),
                    aes(label = external_gene_name), 
                    size = 6/.pt,
                    family = 'Arial', 
                    box.padding = 0.3,
                    #label.padding = 0.1,
                    max.overlaps = 50) +
    stat_quadrant_counts(xintercept = 0,
                         yintercept = 0,
                         size = 6/.pt,
                         family = 'Arial') +
    labs(x = expression(paste(log[2], "(fold ", change[TE_region], ")")),
         y = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    facet_grid(.~tissue) +
    #theme_rob() +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          strip.background = element_rect(fill = tissue.color['skin']),
          strip.text = element_text(size = 10, colour = 'white'))


# =============================== blood =========================================

blood_data <- combined.df$blood %>% 
    filter(!duplicated(te_region_id), !is.na(rna_padj))

blood_faded <- blood_data %>% filter(abs(rna_L2FC) < 0.2)
blood_focus <- blood_data %>% filter(abs(rna_L2FC) > 0.2)

blood_data$tissue <- 'blood'

corr_blood <- ggplot(blood_data, aes(rna_L2FC, gene_L2FC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point_rast(data = blood_faded,
                    alpha = 0.05, 
                    color = '#192E37',
                    size = 0.3) +
    geom_point(data = blood_focus,
               size = 0.3,
               color = '#192E37') +
    geom_smooth(method = 'lm',
                color = '#028763',
                linewidth = 0.3) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'Rho',
             size = 6/.pt,
             family = 'Arial',
             p.accuracy = .Machine$double.xmin,
             #p.digits = 300,
             label.sep = '\n',
             label.x = -3, 
             label.y = 2) +
    xlim(c(-3,3)) +
    ylim(c(-3,3)) +
    geom_text_repel(data = blood_focus %>% filter(abs(rna_L2FC) > 0.75),
                    aes(label = external_gene_name), 
                    size = 6/.pt,
                    family = 'Arial',
                    box.padding = 0.3,
                    #label.padding = 0.1,
                    max.overlaps = 50) +
    facet_grid(.~tissue) +
    stat_quadrant_counts(xintercept = 0,
                         yintercept = 0,
                         size = 6/.pt,
                         family = 'Arial') +
    labs(x = expression(paste(log[2], "(fold ", change[TE_region], ")")),
         y = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          strip.background = element_rect(fill = tissue.color['blood']),
          strip.text = element_text(size = 10, colour = 'white'))


corr_pan <- ggarrange(corr_brain, corr_skin, corr_blood, 
                      nrow = 1,
                      align = c("hv"))

ggsave(filename = paste0(figure_dir, '10_te_island_host_correlation.pdf'),
       device = cairo_pdf,
       plot = corr_pan,
       width = 20,
       height = 7,
       units = 'cm',
       dpi = 300)





