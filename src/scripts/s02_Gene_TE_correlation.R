# The purpose of that script is to make a log2FC correlation plot between genes and their associated TEs.
# A distinction is made between inter- and intragenic TEs.

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}

getIntergenic = FALSE

# Get all intergenic TEs to store them in a bed to make a bedtools closest 
# to get the neighbored genes.
te.annotation <- te.annotation %>% 
    filter(chromosome %in% chr_of_interest) #%>% 
   # mutate(distance = 0)
    

if (getIntergenic) {
    
    intergenic_TEs_bed <- te.annotation %>% 
        filter(position == 'intergenic') %>% 
        dplyr::select(chromosome, start, end, te_id, Kimura, strand) %>% 
        mutate(Kimura = 1)
    
    write.table(intergenic_TEs_bed, 
                file = paste0(table_dir, "intergeneic_tes.bed"),
                row.names = F, 
                col.names = F, 
                quote = F, 
                sep = '\t')
    
    gene_bed <- data.frame(geneRanges) %>% 
        dplyr::select(seqnames, start, end, ensembl_gene_id, width, strand)
    
    write.table(gene_bed,
                file = paste0(table_dir, "mm10_gene_annotation_v102.bed"),
                row.names = F, 
                col.names = F, 
                quote = F, 
                sep = '\t')
    
    # The gene and intergenic TE bed was used to get the closest genes of intergenic TEs. Therefore,
    # bedtools closest -a x -b y -s -d was used
    
}


#########################################################################
##### 
### Correlation between intragenic TEs and their host genes
##  
#
#########################################################################

intragenic_combination <- te.annotation %>% 
    filter(position == 'intronic') %>% 
    dplyr::select(te_id, ensembl_gene_id)

gene_deseq <- fread('results/tables/02_deseq_results_gene.csv')
te_instances_deseq <- fread('results/tables/02_deseq_results_te_instances.csv')

intraList <- list()

##### brain

genes_brain <- gene_deseq %>% filter(tissue == 'brain', 
                               ensembl_gene_id %in% intragenic_combination$ensembl_gene_id,
                               # baseMean > 5,
                               # abs(log2FoldChange) > 0.1
                               ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange)

tes_brain <- te_instances_deseq %>% filter(tissue == 'brain',
                                     te_id %in% intragenic_combination$te_id,
                                     # baseMean > 5,
                                     # abs(log2FoldChange) > 0.1
                                     ) %>% 
    dplyr::select(te_id, log2FoldChange, ensembl_gene_id) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange)


combination_brain <- merge(genes_brain, tes_brain, by = 'ensembl_gene_id')

 # ### Blood

genes_blood <- gene_deseq %>% filter(tissue == 'blood', 
                               ensembl_gene_id %in% intragenic_combination$ensembl_gene_id,
                               #baseMean > 5,
                               #abs(log2FoldChange) > 0.1
                               ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange)

tes_blood <- te_instances_deseq %>% filter(tissue == 'blood',
                                     te_id %in% intragenic_combination$te_id,
                                     #baseMean > 5,
                                     #abs(log2FoldChange) > 0.1
                                     ) %>% 
    dplyr::select(te_id, log2FoldChange, ensembl_gene_id) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange)


combination_blood <- merge(genes_blood, tes_blood, by = 'ensembl_gene_id')

###### skin

genes_skin <- gene_deseq %>% filter(tissue == 'skin', 
                               ensembl_gene_id %in% intragenic_combination$ensembl_gene_id,
                               #baseMean > 5,
                               #abs(log2FoldChange) > 0.1
                               ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange)

tes_skin <- te_instances_deseq %>% filter(tissue == 'skin',
                                     te_id %in% intragenic_combination$te_id,
                                     #baseMean > 5,
                                     #abs(log2FoldChange) > 0.1
                                     ) %>% 
    dplyr::select(te_id, log2FoldChange, ensembl_gene_id) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange)


combination_skin <- merge(genes_skin, tes_skin, by = 'ensembl_gene_id')

intra_combination <- rbind(combination_brain, rbind(combination_blood, combination_skin))

intra_combination$position <- 'intra'
intra_combination$distance <- 0

ggplot(intra_combination, aes(gene_log2FoldChange, te_log2FoldChange)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(method = 'lm') +
    coord_cartesian() +
    facet_grid(.~tissue)


intra_combination$position <- 'intra'
intra_combination$distance <- 0

intra_combination$tissue <- factor(intra_combination$tissue, level = tissues)

pl <- ggplot(intra_combination, aes(gene_log2FoldChange, te_log2FoldChange)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(method = 'lm', color = '#028763') +
    stat_quadrant_counts(xintercept = 0, yintercept = 0) +
    coord_cartesian() +
    facet_grid(.~tissue) +
    theme_rob() +
    theme(text = element_text(family = 'arial'))


pl <- color_strips(pl, 
                   bg_cols = c(tissue.color[2:4]), 
                   text_cols = c( "#ffffff", "#000000","#ffffff"))

pl <- grid.draw(pl)

show(pl)

#########################################################################
##### 
### Correlation between intergenic TEs and their associated (closest)
##  genes
#
#########################################################################

te_gene_distance <- read.csv(paste0(table_dir,
                                    "intragenic_te_gene_distance.csv"))


#### brain

inter_te_deseq_brain <-
    te_instances_deseq %>% filter(tissue == 'brain',
                                  te_id %in% te_gene_distance$te_id,
                                 # baseMean > 5,
                                 # abs(log2FoldChange) > 0.1
                                  ) %>%
    dplyr::select(te_id, log2FoldChange, tissue) %>%
    dplyr::rename(te_log2FoldChange = log2FoldChange)

inter_te_deseq_brain <-
    merge(inter_te_deseq_brain, te_gene_distance, by = 'te_id')

inter_gene_deseq_brain <- gene_deseq %>% filter(tissue == 'brain',
                                                ensembl_gene_id %in% te_gene_distance$ensembl_gene_id,
                                               # baseMean > 5,
                                               # abs(log2FoldChange) > 0.1
                                                ) %>%
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name) %>%
    dplyr::rename(gene_log2FoldChange = log2FoldChange)

inter_combination_brain <-
    merge(inter_te_deseq_brain, inter_gene_deseq_brain, by = 'ensembl_gene_id')


##### blood ####

inter_te_deseq_blood <-
    te_instances_deseq %>% filter(tissue == 'blood',
                                  te_id %in% te_gene_distance$te_id,
                                  #baseMean > 5,
                                  #abs(log2FoldChange) > 0.1
                                  ) %>%
    dplyr::select(te_id, log2FoldChange, tissue) %>%
    dplyr::rename(te_log2FoldChange = log2FoldChange)

inter_te_deseq_blood <-
    merge(inter_te_deseq_blood, te_gene_distance, by = 'te_id')

inter_gene_deseq_blood <- gene_deseq %>% filter(tissue == 'blood',
                                                ensembl_gene_id %in% te_gene_distance$ensembl_gene_id,
                                               # baseMean > 5,
                                                #abs(log2FoldChange) > 0.1
                                                ) %>%
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name) %>%
    dplyr::rename(gene_log2FoldChange = log2FoldChange)

inter_combination_blood <-
    merge(inter_te_deseq_blood, inter_gene_deseq_blood, by = 'ensembl_gene_id')

##### skin ####

inter_te_deseq_skin <-
    te_instances_deseq %>% filter(tissue == 'skin',
                                  te_id %in% te_gene_distance$te_id,
                                  #baseMean > 5,
                                  #abs(log2FoldChange) > 0.1
                                  ) %>%
    dplyr::select(te_id, log2FoldChange, tissue) %>%
    dplyr::rename(te_log2FoldChange = log2FoldChange)

inter_te_deseq_skin <-
    merge(inter_te_deseq_skin, te_gene_distance, by = 'te_id')

inter_gene_deseq_skin <- gene_deseq %>% filter(tissue == 'skin',
                                                ensembl_gene_id %in% te_gene_distance$ensembl_gene_id,
                                               #baseMean > 5,
                                               #abs(log2FoldChange) > 0.1
                                               ) %>%
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name) %>%
    dplyr::rename(gene_log2FoldChange = log2FoldChange)

inter_combination_skin <-
    merge(inter_te_deseq_skin, inter_gene_deseq_skin, by = 'ensembl_gene_id')


inter_combination <- rbind(inter_combination_brain, rbind(inter_combination_skin, inter_combination_blood))


inter_combination$position <- 'inter'

new_order <- names(intra_combination)

inter_combination <- inter_combination[,..new_order]

data <- rbind(inter_combination, intra_combination)

data$tissue <- factor(data$tissue, levels = c('brain', 'skin', 'blood'))
data$position <- factor(data$position, levels = c('intra', 'inter'))

pl <- ggplot(data, aes( gene_log2FoldChange, te_log2FoldChange)) +
    #geom_point(alpha = 0.05, color = '#192E37') +
    geom_point_rast(alpha = 0.05, color = '#192E37') +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    #geom_smooth(method = 'lm', color = '#028763') +
    xlim(c(-10,10)) +
    ylim(c(-10,10)) +
    stat_quadrant_counts(xintercept = 0, yintercept = 0) +
    labs(y = expression(paste(log[2], "(fold ", change[TE], ")")),
         x = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    coord_cartesian() +
    stat_cor(method = "spearman", label.y = 7, label.x = -9) +
    #stat_poly_eq(label.y = 0.7, label.x = 0.1) +
    facet_grid(position~tissue) +
    theme_rob() +
    theme()

pl <- color_strips(pl,
                   bg_cols = c(tissue.color[2:4], "#ffffff","#ffffff"), 
                   text_cols = c( "#ffffff", "#000000","#ffffff", "#000000", "#000000"))

pl <- grid.draw(pl)

show(pl)


### the stuff which is relevant for a submission.

ggsave(
    filename = paste0(figure_dir, 's02_TE_gene_corr.pdf'),
    device = cairo_pdf(),
    plot = last_plot(),
    width = 20,
    height = 15,
    units = "cm",
    dpi = 300)

for (tissue in c('skin', 'blood', 'brain')) {
    
    summary(
        lm(gene_log2FoldChange~te_log2FoldChange, data %>% filter(tissue == tissue, position == 'inter'))
    )
    
    #cor.test(brain_data$gene_L2FC, brain_data$rna_L2FC, method = 'spearman')
    
}

summary(
    lm(gene_log2FoldChange~te_log2FoldChange, data %>% filter(tissue == 'skin', position == 'intra'))
)


#########
# RESIS #
#########

ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_2/S2_TE_gene_correlation.pdf',
    plot = last_plot(),
    device = cairo_pdf,
    width = 20,
    height = 15,
    units = "cm",
    dpi = 300)

write.csv(data,
          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_2/S2_TE_gene_correlation.csv')

#########################################################################
##### 
### Correlation between exonic TEs and their host genes
##
#
#########################################################################


exonic_combination <- te.annotation %>% 
    filter(position == 'exonic') %>% dplyr::select(te_id, ensembl_gene_id)

gene_deseq <- fread('results/tables/02_deseq_results_gene.csv')
te_instances_deseq <- fread('results/tables/02_deseq_results_te_instances.csv')

exonicList <- list()

##### brain

genes_brain <- gene_deseq %>% filter(tissue == 'brain', 
                               ensembl_gene_id %in% exonic_combination$ensembl_gene_id,
                               # baseMean > 5,
                               # abs(log2FoldChange) > 0.1
                               ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, baseMean, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange,
                  gene_baseMean = baseMean)

exonic_tes_brain <- te_instances_deseq %>% filter(tissue == 'brain',
                                     te_id %in% exonic_combination$te_id,
                                      #baseMean > 10,
                                     # abs(log2FoldChange) > 0.1
                                     ) %>% 
    dplyr::select(te_id, log2FoldChange, ensembl_gene_id, baseMean) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange,
                  te_baseMean = baseMean)


exonic_combination_brain <- merge(genes_brain, exonic_tes_brain, by = 'ensembl_gene_id')

#########
# Blood #
#########

genes_blood <- gene_deseq %>% filter(tissue == 'blood', 
                               ensembl_gene_id %in% exonic_combination$ensembl_gene_id,
                               # baseMean > 5,
                               # abs(log2FoldChange) > 0.1
                               ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, baseMean, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange,
                  gene_baseMean = baseMean)

exonic_tes_blood <- te_instances_deseq %>% filter(tissue == 'blood',
                                     te_id %in% exonic_combination$te_id,
                                     # baseMean > 10,
                                     # abs(log2FoldChange) > 0.1
                                     ) %>% 
    dplyr::select(te_id, log2FoldChange, ensembl_gene_id, baseMean) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange,
                  te_baseMean = baseMean)


exonic_combination_blood <- merge(genes_blood, exonic_tes_blood, by = 'ensembl_gene_id')

#########
# skin #
#########

genes_skin <- gene_deseq %>% filter(tissue == 'skin', 
                               ensembl_gene_id %in% exonic_combination$ensembl_gene_id,
                               # baseMean > 5,
                               # abs(log2FoldChange) > 0.1
                               ) %>% 
    dplyr::select(ensembl_gene_id, log2FoldChange, external_gene_name, baseMean, tissue) %>% 
    dplyr::rename(gene_log2FoldChange = log2FoldChange,
                  gene_baseMean = baseMean)


exonic_tes_skin <- te_instances_deseq %>% filter(tissue == 'skin',
                                     te_id %in% exonic_combination$te_id,
                                      #baseMean > 10,
                                     # abs(log2FoldChange) > 0.1
                                     ) %>% 
    dplyr::select(te_id, log2FoldChange, ensembl_gene_id, baseMean) %>% 
    dplyr::rename(te_log2FoldChange = log2FoldChange,
                  te_baseMean = baseMean)


exonic_combination_skin <- merge(genes_skin, exonic_tes_skin, by = 'ensembl_gene_id')

exonic_combination <- rbind(exonic_combination_brain, rbind(exonic_combination_blood, exonic_combination_skin))

exonic_combination$position <- 'exonic'
exonic_combination$distance <- 0
exonic_combination$tissue <- factor(exonic_combination$tissue, levels = tissues)

pl <- ggplot(exonic_combination, aes(gene_log2FoldChange, te_log2FoldChange)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(method = 'lm', color = '#028763') +
    stat_quadrant_counts(xintercept = 0, yintercept = 0) +
    stat_poly_eq(label.y = 0.8, label.x = 0.1) +
    geom_text_repel(aes(label = external_gene_name), 
                    size = 10*0.36,
                    box.padding = 0.7,
                    #label.padding = 0.1,
                    max.overlaps = 50) +
    coord_cartesian() +
    facet_grid(.~tissue) +
    theme_rob() +
    theme(text = element_text(family = 'arial'))


pl <- color_strips(pl, 
                   bg_cols = c(tissue.color), 
                   text_cols = c( "#ffffff", "#000000","#ffffff"))

pl <- grid.draw(pl)

show(pl)


ggplot(exonic_combination, aes(gene_baseMean, te_baseMean, fill = tissue)) +
    geom_point()
