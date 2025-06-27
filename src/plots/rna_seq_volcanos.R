
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_env()
aging_tes::load_plotting_env()


# ------------------------------------------------------------------------------
# Volcano plots for all TEs and tissues
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

# young TEs (kimura <= 5)

young_TEs <- deseq.te.merged %>% 
    filter(Kimura <= 5)


# change order of tissues
young_TEs$tissue <- factor(young_TEs$tissue, 
                                 levels = c('brain', 'skin', 'blood'))

volcanoPlot(cutPvalue(young_TEs), FDR = 0.05, facet = 'tissue')


