
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_female_env()
aging_tes::load_plotting_env()

# ------------------------------------------------------------------------------
# Volcano plots for all TEs and tissues
# ------------------------------------------------------------------------------

deseq.te.merged <- fread(paste0(table_dir, deseq_results_te_csv))

# change order of tissues
deseq.te.merged$tissue <- factor(deseq.te.merged$tissue, 
                                 levels = c('brain', 'skin'))

volcano.te <- volcanoPlot(cutPvalue(deseq.te.merged), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')


volcano.te <- color_strips(volcano.te, 
                           bg_cols = tissue.color[2:3], 
                           text_cols = c( "#ffffff", "#000000"))


# Save figure with metadata to index
meta <- list(name = 'te_instances_female_volcano',
             description = 'Volcano plot of all expressed transposable elements in female brain and skin',
             tags = c('expression', 'rna-seq', 'female'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin')),
             script = 'rna_seq_volcano_female.R'
)

fig_index(plot = volcano.te,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 9.5,
          height = 5.5,
          dpi = 300,
          format = 'pdf')

# Display plot

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
                           levels = c('brain', 'skin'))

mid_age$tissue <- factor(mid_age$tissue, 
                         levels = c('brain', 'skin'))

old_TEs$tissue <- factor(old_TEs$tissue, 
                         levels = c('brain', 'skin'))

# Volcano plots for each group of TEs

volcano.young <- volcanoPlot(cutPvalue(young_TEs), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.young <- color_strips(volcano.young, 
                              bg_cols = tissue.color[2:3], 
                              text_cols = c( "#ffffff", "#000000"))

# Save figure with metadata to index
meta <- list(name = 'te_young_instances_female_volcano',
             description = 'Volcano plot of young transposable elements (kimura <= 5) in female brain and skin',
             tags = c('expression', 'rna-seq', 'female'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin'), kimura = '<=5'),
             script = 'rna_seq_volcano_female.R'
)

fig_index(plot = volcano.young,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 9.5,
          height = 5.5,
          dpi = 300,
          format = 'pdf')

volcano.mid <- volcanoPlot(cutPvalue(mid_age), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.mid <- color_strips(volcano.mid, 
                            bg_cols = tissue.color[2:3], 
                            text_cols = c( "#ffffff", "#000000"))

# Save figure with metadata to index
meta <- list(name = 'te_midage_instances_female_volcano',
             description = 'Volcano plot of mid aged expressed transposable elements (kimura between 5 and 25) in female brain and skin',
             tags = c('expression', 'rna-seq', 'female'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin'), kimura = '5-25'),
             script = 'rna_seq_volcano_female.R'
)

fig_index(plot = volcano.mid,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 9.5,
          height = 5.5,
          dpi = 300,
          format = 'pdf')

volcano.old <- volcanoPlot(cutPvalue(old_TEs), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.old <- color_strips(volcano.old, 
                            bg_cols = tissue.color[2:3], 
                            text_cols = c( "#ffffff", "#000000"))

# Save figure with metadata to index
meta <- list(name = 'te_old_instances_female_volcano',
             description = 'Volcano plot of old expressed transposable elements (Kimura > 25) in female brain and skin',
             tags = c('expression', 'rna-seq', 'female'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin'), kimura = '>25'),
             script = 'rna_seq_volcano_female.R'
)

fig_index(plot = volcano.old,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 9.5,
          height = 5.5,
          dpi = 300,
          format = 'pdf')
