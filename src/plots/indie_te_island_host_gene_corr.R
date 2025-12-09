# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_annotations()

load_gene_ranges()
load_te_island_annotation()

gene_te_island_association_merged = data.table::fread(file = paste0(table_dir, '02_deseq_results_te_island_extended.csv'))


# ------------------------------ Correlation Plots -----------------------------
#
# The tables may contain multiple entries per TE region, as there is a 
# possibility that more than one CAGE-peak may exists per TE region. 
# Consequently, a filtered table is used for the log2FC correlation plots.
# 
# ============================== Brain =========================================

brain_data <- gene_te_island_association_merged %>% 
    filter(tissue == 'brain', !duplicated(te_island_id), individually, gene_association == 'intragenic')

faded <- brain_data %>% filter(abs(log2FoldChange) < 0.1)
focus <- brain_data %>% filter(abs(log2FoldChange) > 0.1)


ggplot(brain_data, aes(log2FoldChange, gene_log2FoldChange)) +
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

skin_data <- gene_te_island_association_merged %>% 
    filter(tissue == 'skin', !duplicated(te_island_id))

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

