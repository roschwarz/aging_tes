# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}


te.family.enrichment <- loadRdata(file = paste0(data_dir, "28_te_family_quant_enrichment.Rdata"))

#### Kimura distance stuff


b1_family <- te.annotation %>% 
    filter(super_family == "Alu") %>% 
    mutate(Kimura = as.numeric(Kimura))


sine_as_bg <- te.annotation %>% 
    filter(order == 'SINE') %>% 
    mutate(Kimura = as.numeric(Kimura),
           family = 'background')


kim_df <- rbind(b1_family, sine_as_bg) 

family_order <- kim_df %>% 
    group_by(family) %>% 
    dplyr::summarise(kimura_mean = mean(as.numeric(Kimura), na.rm = TRUE), .groups = 'drop') %>% 
    arrange(kimura_mean) %>% 
    pull(family)

kim_df$family <- factor(kim_df$family, levels = family_order)

kimura_pl <- ggplot(kim_df, aes(family, as.numeric(Kimura), fill = family)) + 
    geom_violin() +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6) +
    labs(y = "Kimura distance") +
    coord_flip() +
    theme_rob() +
    theme(legend.position = 'None')

#### correlation

kim_mean <- kim_df %>% 
    group_by(family) %>% 
    summarise(kim_age = mean(Kimura))

enriched_b1s <- te.family.enrichment %>% filter(category %in% family_order)

enriched_b1s <- merge(enriched_b1s, kim_mean, by.x = 'category', by.y = 'family')

ggplot(enriched_b1s, aes(kim_age, ratio)) +
    geom_point(aes(color = tissue)) +
    geom_smooth(method = 'lm', , color = '#028763', linewidth = 0.3) +
    scale_color_manual(values = tissue.color[2:4]) + 
    stat_cor(method = 'spearman',
             cor.coef.name = 'Rho',
             p.accuracy = .Machine$double.xmin,
             size = 8/.pt,
             label.sep = '\n',
             label.y = 0.5, 
             label.x = 8,
             family = 'Arial',)  +
    labs(x = "Kimura distance",
         y = "log(odds ratio)") +
        theme(legend.position = c(0.82,0.75),
              legend.background = element_rect(color = 'black', linewidth = 0.3),
              axis.title = element_text(size = 8),
              axis.text = element_text(size = 8))


ggsave(filename = paste0(figure_dir, 's09_Kimura_corr.pdf'),
       device = cairo_pdf,
       plot = last_plot(),
       width = 7,
       height = 7,
       units = 'cm',
       dpi = 300)


#########
# RESIS #
#########

write.csv(enriched_b1s,
          './manuscripts/nature_aging/resis/figures/supplemental_figure_9/S9A_Kimura_corr.csv')

ggsave(filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_9/S9A_Kimura_corr.pdf',
       device = cairo_pdf,
       plot = last_plot(),
       width = 7,
       height = 7,
       units = 'cm',
       dpi = 300)




#############
# Pie-Chart #
#############
# 1. Get all DETEs
# 2. add column where it is noted to which category the dete contains:
#    a - DETE contains a CAGE signal
#    b - DETE is part of a canoncial TE island
#    c - DETE is co-regulated with their host gene
#    d - all DETEs that are not assigned are categorized as other
#
# Keep the order as named above so when a DETE is intersected with a TE island and a CAGE-signal than
# it is assigned to the CAE signal category.

if (!exists("deseq.te.merged")) {
    
    deseq.te.merged <- read.csv(paste0(table_dir, '02_deseq_results_te_instances.csv'))
    
}

cage.tes <- read.csv(paste0(table_dir, '07_deseq_cage_te.csv')) %>% 
    filter(tissue == 'brain') %>% 
    dplyr::pull(te.te.id)

island.tes <- blackRcloud::loadRdata(paste0(data_dir, "30_TE_island_composition.Rdata")) %>% 
    filter(tissue == 'brain') %>% 
    pull(te_id)


dete.gene <- read.csv(paste0(table_dir, '02_deseq_results_gene.csv')) %>% filter(tissue == 'brain', 
                                                                                 padj <= 0.1 ) %>% 
    dplyr::pull(ensembl_gene_id)


gene.tes <- deseq.te.merged %>% 
    filter(ensembl_gene_id %in% dete.gene) %>% 
    dplyr::pull(te_id)


detes <- deseq.te.merged %>%
    filter(tissue == 'brain', padj <= 0.05)

detes$category <- 'other' 

detes <- detes %>%  
    mutate(category = case_when(te_id %in% gene.tes ~ 'host gene', .default = category))

detes <- detes %>% 
    mutate(category = case_when(te_id %in% island.tes ~ 'TE island', .default = category))

detes <- detes %>% 
    mutate(category = case_when(te_id %in% cage.tes ~ 'CTSS', .default = category))

# foo <- foo %>% 
#     mutate(category = case_when(te_id %in% cage.tes ~ 'CTSS',
#                                 te_id %in% island.tes ~ 'TE island',
#                                 te_id %in% gene.tes ~ 'host Gene',
#                                 .default = 'other'))

detes$category <- factor(detes$category, levels = c("CTSS", "TE island", "host gene", "other"))

ggplot(detes, aes(category)) +
    geom_bar(position = 'stack') +
    labs(y = 'counts of DETEs') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme(legend.position = "None",
          legend.background = element_rect(color = 'black', linewidth = 0.3),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          axis.title.x = element_blank())


ggsave(filename = paste0(figure_dir, 's09_DETE_regulation.pdf'),
       device = cairo_pdf,
       plot = last_plot(),
       width = 7,
       height = 7,
       units = 'cm',
       dpi = 300)


#########
# RESIS #
#########

ggsave(filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_9/S9_DETE_regulation.pdf',
       device = cairo_pdf,
       plot = last_plot(),
       width = 7,
       height = 7,
       units = 'cm',
       dpi = 300)

write.csv(detes, './manuscripts/nature_aging/resis/figures/supplemental_figure_9/S9_DETE_regulation.csv')

