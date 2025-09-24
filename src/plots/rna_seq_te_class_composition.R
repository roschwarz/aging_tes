# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_rna_seq_env()
aging_tes::load_annotations()

load_te_ranges()


# ------------------------------------------------------------------------------
# Step 1: Calculate the proportion of TE classes for all TEs (background)
# ------------------------------------------------------------------------------
deseq.te.merged <- fread(paste0(table_dir, deseq_results_te_csv))

expressed.tes <- deseq.te.merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'expressed')

background <- transGrange(teRanges) %>% 
    splitTEID('te.id') %>% 
    filter(order %in% orders.of.interest) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           tissue = 'background',
           group = 'expressed')

# ------------------------------------------------------------------------------
# Step 2: Calculate the proportion of TE classes for expressed and differentially expressed TEs
# ------------------------------------------------------------------------------

dete.up <- deseq.te.merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= FDR, log2FoldChange > 0) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'up')

dete.down <- deseq.te.merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= FDR, log2FoldChange < 0) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'down')


# ------------------------------------------------------------------------------
# Step 3: Merge the data frames 
# ------------------------------------------------------------------------------

df.classes <- rbind(background, expressed.tes, dete.up, dete.down)

df.classes$group = factor(df.classes$group, levels = c('expressed','down','up'))
df.classes$tissue = factor(df.classes$tissue, levels = c('background', 'brain', 'skin', 'blood'))


# ------------------------------------------------------------------------------
# Step 4: Create the plot
# ------------------------------------------------------------------------------

df.class.updated <- df.classes %>% 
    mutate(x = paste0(tissue, '_', group))

order.of.x = c('background_expressed',
               'brain_expressed',
               'skin_expressed',
               'blood_expressed',
               'brain_down',
               'skin_down',
               'blood_down',
               'brain_up',
               'skin_up',
               'blood_up')

df.class.updated$x <- factor(df.class.updated$x, levels = order.of.x)


pl.te.percent.updated <- ggplot(df.class.updated, aes(x, percent, fill = order)) +
    geom_flow(aes(alluvium = order), 
              alpha = .3, 
              color = 'white',
              curve_type = 'linear', 
              width = .3) +
    geom_col(width = .7, color = "white") +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0,0))) +
    scale_x_discrete(expand = c(0,0),labels = c('background_expressed' = 'background',
                                               'brain_expressed' = 'brain',
                                               'skin_expressed' = 'skin',
                                               'blood_expressed' = 'blood',
                                               'brain_down' = 'brain',
                                               'skin_down' = 'skin',
                                               'blood_down' = 'blood',
                                               'brain_up' = 'brain',
                                               'skin_up' = 'skin',
                                               'blood_up' = 'blood')) +
    scale_fill_manual(values = order.color,
                      name = 'TE class:') +
    labs(y = 'Percent of\n(differentially) expressed TEs') +
    geom_vline(xintercept = c(1.5, 4.5, 7.5), color = "black", linewidth = .5) +
    theme_rob(base_family = 'arial',
              base_size = 10) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          legend.position = 'bottom',
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1
          ),
          legend.background = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-15,-0, 0,-10),
          legend.key.size = unit(0.5,"line"),
    )

# Save figure with metadata to index
meta <- list(name = 'percentage_expr_te_class',
             description = 'Proportion of TE classes among expressed or differentially expressed TE instances compared to the background',
             tags = c('expression', 'rna-seq'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin', 'blood')),
             script = 'rna_seq_te_class_composition.R'
)

fig_index(plot = pl.te.percent.updated,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 9.2,
          height = 7.5,
          dpi = 300,
          format = 'pdf')

#ggsave(pl.te.percent.updated,
#       filename = paste0(figure_dir, 'panel_1_expr_te_class_updated_4x3_300.pdf'),
#       device = cairo_pdf,
#       width = 9.2,
#       height = 7.5,
#       units = "cm",
#       dpi = 300
#)

