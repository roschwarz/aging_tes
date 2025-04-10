
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_env()
aging_tes::load_plotting_env()


# ------------------------------------------------------------------------------
# Volcano plots for all TEs and tissues
# ------------------------------------------------------------------------------

# change order of tissues
deseq.mixed.merged$tissue <- factor(deseq.mixed.merged$tissue, 
                                 levels = c('brain', 'skin', 'blood'))


deseq_mixed_tes <- deseq.mixed.merged %>% filter(grepl("^chr", te_id))

volcano_mixed_tes <- volcanoPlot(cutPvalue(deseq_mixed_tes), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')


volcano_mixed_tes <- color_strips(volcano_mixed_tes, 
                           bg_cols = tissue.color[2:4], 
                           text_cols = c( "#ffffff", "#000000","#ffffff"))

ggsave(plot = volcano_mixed_tes,
       filename = paste0(figure_dir, 'te_instances_mixed_volcano_9.5x5.5_300.pdf'),
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

pl <- grid.draw(volcano_mixed_tes)

show(pl)