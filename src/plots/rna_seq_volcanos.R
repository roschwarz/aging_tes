
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_env()
aging_tes::load_plotting_env()


# ------------------------------------------------------------------------------
# Volcano plots for all TEs and tissues in male
# ------------------------------------------------------------------------------

deseq.te.merged <- fread(paste0(table_dir, deseq_results_te_csv))

# change order of tissues
deseq.te.merged$tissue <- factor(deseq.te.merged$tissue, 
                                 levels = c('brain', 'skin', 'blood'))

volcano.te <- volcanoPlot(cutPvalue(deseq.te.merged), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')


volcano.te <- color_strips(volcano.te, 
                           bg_cols = tissue.color[2:4], 
                           text_cols = c( "#ffffff", "#000000","#ffffff"))

ggsave(plot = volcano.te,
       filename = paste0(figure_dir, 'te_instances_volcano_9.5x5.5_300.pdf'),
       device = cairo_pdf,
       width = 9.5,
       height = 5.5,
       units = "cm",
       dpi = 300
)


if (requireNamespace("grid", quietly = TRUE)) {
    library(grid)
}else{
    message("grid is not installed, please install for printing the merged plot")
}

pl <- grid.draw(volcano.te)

show(pl)


# ------------------------------------------------------------------------------
# Volcano plots separated by age of TEs using Kimura distances
# ------------------------------------------------------------------------------
# There are 3 groups of TEs based on Kimura distances:
# young TEs (kimura <= 5)
# mid-age TEs (5 < kimura <= 25)
# old TEs (kimura > 25)  

young_TEs <- deseq.te.merged %>% 
    filter(Kimura <= 5)

mid_age <- deseq.te.merged %>% 
    filter(Kimura > 5 & Kimura <= 25)

old_TEs <- deseq.te.merged %>%
    filter(Kimura > 25)

# Change order of tissues

young_TEs$tissue <- factor(young_TEs$tissue, 
                                 levels = c('brain', 'skin', 'blood'))

mid_age$tissue <- factor(mid_age$tissue, 
                                 levels = c('brain', 'skin', 'blood'))

old_TEs$tissue <- factor(old_TEs$tissue, 
                                 levels = c('brain', 'skin', 'blood'))

# Volcano plots for each group of TEs

volcano.young <- volcanoPlot(cutPvalue(young_TEs), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.young <- color_strips(volcano.young, 
                           bg_cols = tissue.color[2:4], 
                           text_cols = c( "#ffffff", "#000000","#ffffff"))

ggsave(plot = volcano.young,
       filename = paste0(figure_dir, 'te_young_volcano_9.5x5.5_300.pdf'),
       device = cairo_pdf,
       width = 9.5,
       height = 5.5,
       units = "cm",
       dpi = 300
)

volcano.mid <- volcanoPlot(cutPvalue(mid_age), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.mid <- color_strips(volcano.mid, 
                                  bg_cols = tissue.color[2:4], 
                                  text_cols = c( "#ffffff", "#000000","#ffffff"))

ggsave(plot = volcano.mid,
       filename = paste0(figure_dir, 'te_midage_volcano_9.5x5.5_300.pdf'),
       device = cairo_pdf,
       width = 9.5,
       height = 5.5,
       units = "cm",
       dpi = 300
)

volcano.old <- volcanoPlot(cutPvalue(old_TEs), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.old <- color_strips(volcano.old, 
                                  bg_cols = tissue.color[2:4], 
                                  text_cols = c( "#ffffff", "#000000","#ffffff"))

ggsave(plot = volcano.old,
       filename = paste0(figure_dir, 'te_old_volcano_9.5x5.5_300.pdf'),
       device = cairo_pdf,
       width = 9.5,
       height = 5.5,
       units = "cm",
       dpi = 300
)


if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_env()
aging_tes::load_plotting_env()


# ------------------------------------------------------------------------------
# Volcano plots for all TEs and tissues in male and female
# ------------------------------------------------------------------------------

deseq_te_merged_male <- fread(paste0(table_dir, deseq_results_te_csv_male)) %>% 
    mutate(sex = 'male',
           sex_tissue = paste0(sex, "_",tissue))

deseq_te_merged_female <- fread(paste0(table_dir, deseq_results_te_csv_female)) %>% 
    mutate(sex = 'female',
           sex_tissue = paste0(sex, "_", tissue))

deseq_tes <- rbind(deseq_te_merged_male, deseq_te_merged_female)


# === Collect some Numbers ===

deseq_tes %>% filter(!is.na(padj)) %>% group_by(tissue, sex) %>% summarize(count = n())
deseq_tes %>% filter(!is.na(padj)) %>% group_by(sex) %>% summarize(count = n())
deseq_tes %>% filter(padj <= FDR) %>% group_by(sex, tissue) %>% summarize(count = n())


# change order of tissues
deseq_tes$sex_tissue <- factor(deseq_tes$sex_tissue, 
                                 levels = c('male_brain', 'male_skin', 'male_blood', 'female_brain', 'female_skin'))

volcano_te <- volcanoPlot(cutPvalue(deseq_tes), FDR = 0.05, "sex_tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')

tissue_sex_colors = c(male_brain = "#58B2AA",
                      male_skin = "#EFA081",
                      male_blood = "#685299",
                      female_brain = "#58B2AA",
                      female_skin = "#EFA081")

volcano_te <- color_strips(volcano_te, 
                           bg_cols = tissue_sex_colors, 
                           text_cols = c( "#ffffff", "#000000","#ffffff", "#ffffff", "#000000"))


# Save figure with metadata to index
meta <- list(name = 'te_instances_volcano_both_sexes',
             description = 'Volcano plot of all expressed transposable elements in male and female brain, skin and blood',
             tags = c('expression', 'rna-seq', 'TE'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin', 'blood'), sex = c('male', 'female')),
             script = 'rna_seq_volcanoR'
)

fig_index(plot = volcano_te,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 17.5,
          height = 5.5,
          dpi = 300,
          format = 'pdf')
