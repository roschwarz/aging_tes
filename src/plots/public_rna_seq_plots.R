
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_public_data_env()
aging_tes::load_plotting_env()
aging_tes::load_annotations()
load_te_annotation()

color_set <- unlist(tissue.color[c('Gastrocnemius muscle', "liver", 'White adipose tissue')])
strip_text_color <- c( "#ffffff","#000000", "#000000", "#ffffff")

# ------------------------------------------------------------------------------
# Volcano plots for public data
# ------------------------------------------------------------------------------

# all transposable elements

deseq_results_public <- fread(paste0(rna_seq_deseq_dir, "deseq_results_te_instances_public.csv"))

deseq_results_public_female <- deseq_results_public %>% filter(sex == 'female')

volcano_female <- volcanoPlot(deseq_results_public_female, FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    facet_grid(sex~tissue) +
    theme(legend.position = 'None')

volcano_female <- color_strips(volcano_female, 
                             bg_cols = c(color_set,"#136f63"), 
                             text_cols = strip_text_color)

deseq_results_public_male <- deseq_results_public %>% filter(sex == 'male')

volcano_male <- volcanoPlot(deseq_results_public_male, FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    facet_grid(sex~tissue) +
    theme(legend.position = 'None') 

volcano_male <- color_strips(volcano_male, 
                           bg_cols = c(color_set, "#ffba08"), 
                           text_cols = strip_text_color)


pl_all_TEs <- gridExtra::grid.arrange(volcano_female, volcano_male, nrow = 2)


# young transposable elements (Kimura distance <= 5)

young_TEs <- deseq_results_public %>% 
    filter(Kimura <= 5)

deseq_results_public_female_young_TEs <- young_TEs %>% filter(sex == 'female')

volcano_female_young_TEs <- volcanoPlot(deseq_results_public_female_young_TEs, FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    facet_grid(sex~tissue) +
    theme(legend.position = 'None')

volcano_female_young_TEs <- color_strips(volcano_female_young_TEs, 
                               bg_cols = c(color_set,"#136f63"), 
                               text_cols = strip_text_color)

deseq_results_public_male_young_TEs <- young_TEs %>% filter(sex == 'male')

volcano_male_young_TEs <- volcanoPlot(deseq_results_public_male_young_TEs, FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    facet_grid(sex~tissue) +
    theme(legend.position = 'None') 

volcano_male_young_TEs <- color_strips(volcano_male_young_TEs, 
                             bg_cols = c(color_set, "#ffba08"), 
                             text_cols = strip_text_color)

pl_young_TEs <- gridExtra::grid.arrange(volcano_female_young_TEs, volcano_male_young_TEs, nrow = 2)


public_volcano_plots <- gridExtra::grid.arrange(pl_all_TEs, pl_young_TEs, nrow = 2)

ggsave(plot = public_volcano_plots,
       filename = paste0(figure_dir, 'public_te_instances_volcano_13x15_300.pdf'),
       device = cairo_pdf,
       width = 13,
       height = 15,
       units = "cm",
       dpi = 300
)


# ------------------------------------------------------------------------------
# Overlap of differentially expressed TEs in female and male
# ------------------------------------------------------------------------------

male <- deseq_results_public_male %>% 
    filter(padj <= 0.1, baseMean > 5) %>% 
    dplyr::select(te_id, log2FoldChange, padj, tissue) %>% 
    dplyr::rename(male.log2FoldChange = log2FoldChange,
                  male.padj = padj)
    
female <- deseq_results_public_female %>% 
    filter(padj <= 0.1, baseMean > 5) %>% 
    dplyr::select(te_id, log2FoldChange, padj, tissue) %>% 
    dplyr::rename(female.log2FoldChange = log2FoldChange,
                  female.padj = padj)

expressed_tes <- merge(male, female, by = c('te_id', 'tissue'))

ggplot(expressed_tes, aes(female.log2FoldChange, male.log2FoldChange, colour = tissue)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_bw()

# ------------------------------------------------------------------------------
# Overlap of differentially expressed TEs in female and male in liver
# ------------------------------------------------------------------------------

intersect(deseq_results_public_male$te_id, deseq_results_public_female$te_id)

deg_male_liver <- deseq_results_public_male %>% 
    filter(padj <= 0.05, tissue == "Liver") %>% 
    dplyr::select(te_id, log2FoldChange) %>% 
    dplyr::rename(male.log2FoldChange = log2FoldChange)

deg_female_liver <- deseq_results_public_female %>% 
    filter(padj <= 0.05, tissue == "Liver") %>% 
    dplyr::select(te_id, log2FoldChange) %>% 
    dplyr::rename(female.log2FoldChange = log2FoldChange)


deg <- merge(deg_male_liver, deg_female_liver, by = 'te_id')

deg <- merge(deg, te.annotation, by = 'te_id')

ggplot(deg, aes(male.log2FoldChange, female.log2FoldChange)) +
    geom_point(aes(shape = position, color = order)) +
    geom_smooth(method = 'lm') +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    theme_bw()


