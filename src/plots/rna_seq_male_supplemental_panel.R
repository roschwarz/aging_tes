
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_env()
aging_tes::load_plotting_env()
aging_tes::load_annotations()

load_te_ranges()


deseq_te_merged <- fread(paste0(table_dir, deseq_results_te_csv))

# change order of tissues
deseq_te_merged$tissue <- factor(deseq_te_merged$tissue, 
                                 levels = c('brain', 'skin', 'blood'))

# ------------------------------------------------------------------------------------------------------------
# Volcano plots for different aged TEs
# ------------------------------------------------------------------------------------------------------------
brain_detes_pl <- create_te_analysis_panel(deseq_data = deseq_te_merged %>% filter(tissue == 'brain'), 
                                           tissue_name = 'brain', 
                                           order_colors = order.color,
                                           fdr_threshold = 0.05,
                                           base_font_size = 10,
                                           min_font_size = 8,
                                           font_family = 'sans',
                                           section = 'bottom',
                                           common.legend = FALSE)

skin_detes_pl <- create_te_analysis_panel(deseq_data = deseq_te_merged %>% filter(tissue == 'skin'), 
                                          tissue_name = 'skin', 
                                          order_colors = order.color,
                                          fdr_threshold = 0.05,
                                          base_font_size = 10,
                                          min_font_size = 8,
                                          font_family = 'sans',
                                          section = 'bottom',
                                          common.legend = FALSE)

blood_detes_pl <- create_te_analysis_panel(deseq_data = deseq_te_merged %>% filter(tissue == 'blood'), 
                                           tissue_name = 'blood', 
                                           order_colors = order.color,
                                           fdr_threshold = 0.05,
                                           base_font_size = 10,
                                           min_font_size = 8,
                                           font_family = 'sans',
                                           section = 'bottom',
                                           common.legend = FALSE)

combined_detes_pl <- ggarrange(brain_detes_pl, 
                               skin_detes_pl,
                               blood_detes_pl,
                               ncol = 1,
                               nrow = 3,
                               common.legend = TRUE,
                               legend = 'bottom')

# ------------------------------------------------------------------------------------------------------------
# Volcano plots for TEs when TEs and genes are in the count matrix used by DESeq2
# ------------------------------------------------------------------------------------------------------------

deseq_mixed_merged <- data.table::fread(paste0(table_dir, deseq_results_mixed_csv))

deseq_mixed_merged$tissue <- factor(deseq_mixed_merged$tissue, 
                                    levels = c('brain', 'skin', 'blood'))


deseq_mixed_tes <- deseq_mixed_merged %>% filter(grepl("^chr", te_id))

volcano_mixed_tes <- volcanoPlot(cutPvalue(deseq_mixed_tes), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')


volcano_mixed_tes <- color_strips(volcano_mixed_tes, 
                                  bg_cols = tissue.color[2:4], 
                                  text_cols = c( "#ffffff", "#000000","#ffffff"))

# ------------------------------------------------------------------------------------------------------------
# Overlap of TEs
# ------------------------------------------------------------------------------------------------------------

expressed_TEs <- sapply(c('brain', 'skin', 'blood'), function(x){
    
    deseq_te_merged %>% 
        filter(tissue == x, !is.na(padj)) %>% 
        pull(te_id) %>% 
        unique()
})


venn <- venn.diagram(
    x = expressed_TEs,
    category.names = names(expressed_TEs),
    
    
    # Circles
    lwd = 2,  
    fill = tissue.color[2:4], #c('#264653', '#2A9D8F',  '#E9C46A'),
    alpha = c(0.7, 0.7, 0.7),
    
    
    # Number
    cex = 1, # font size
    fontface = "bold",
    fontfamily = "arial",
    # 
    # # Set names
    cat.cex = 1.5,
    cat.default.pos = "outer",
    cat.fontface = "bold",
    cat.fontfamily = "arial",
    # cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055, 0.055),
    # main = header,
    scaled = F,
    print.mode = c("raw", "percent"),
    
    # Output
    filename = NULL, #paste0(figures, 'Panel_2C_VennDiagram.svg'),
    imagetype = "svg",
    output = FALSE,
    width = 200,
    height = 500,
    resolution = 300,
    disable.logging = TRUE
    
)

venn <- grid.grabExpr(grid.draw(venn))

venn <- ggplotify::as.ggplot(venn)


# ------------------------------------------------------------------------------------------------------------
# Categorize expressed TEs
# ------------------------------------------------------------------------------------------------------------

# Step 1: Calculate the proportion of TE classes for all TEs (background)

expressed_tes <- deseq_te_merged %>% 
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

dete_up <- deseq_te_merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= 0.05, log2FoldChange > 0) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'up')

dete_down <- deseq_te_merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= 0.05, log2FoldChange < 0) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'down')

df_classes <- rbind(background, expressed_tes, dete_up, dete_down)

df_classes$group = factor(df_classes$group, levels = c('expressed','down','up'))
df_classes$tissue = factor(df_classes$tissue, levels = c('background', 'brain', 'skin', 'blood'))
df_classes <- df_classes %>% 
    mutate(x_axis = paste0(tissue, '_', group))

order_of_x_axis = c('background_expressed',
               'brain_expressed',
               'skin_expressed',
               'blood_expressed',
               'brain_down',
               'skin_down',
               'blood_down',
               'brain_up',
               'skin_up',
               'blood_up')

df_classes$x <- factor(df_classes$x, levels = order_of_x)


pl_te_percent_updated <- ggplot(df_classes, aes(x, percent, fill = order)) +
    geom_flow(aes(alluvium = order), 
              alpha = .3, 
              color = 'white',
              curve_type = 'linear', 
              width = .3) +
    geom_col(width = .7, color = "white") +
    #facet_grid(.~group, scales = 'free_x',space='free') +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0,0))) +
    scale_x_discrete(expand = c(0,0),labels= c('background_expressed' = 'background',
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
    # geom_text(aes(y = 1, x = 2, label = c('expressed'))) +
    #theme(plot.margin = margin(20, 2, 2, 2, unit = "pt"))
    geom_vline(xintercept = c(1.5, 4.5, 7.5), color = "black", linewidth = .5) +
    #geom_hline(yintercept = 1, color = "black", linewidth = .5) +
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
          #plot.margin = margin(20, 2, 2, 2, unit = "pt")
    )

# ------------------------------------------------------------------------------------------------------------
# Create the Panel
# ------------------------------------------------------------------------------------------------------------

last_row <- ggarrange(pl_te_percent_updated, venn, ncol = 2, nrow = 1, widths = c(1.2, 0.8), heights = c(1,0.8))

supplemental_pannel_male <- ggarrange(brain_detes_pl, 
          skin_detes_pl,
          blood_detes_pl,
          volcano_mixed_tes,
          last_row,
          nrow = 5, ncol = 1)

meta <- list(name = 'male _supplemental_panel',
             description = 'Volcano differentially expressed transposable elements separated by their age. Volcanos with gnes considered and overview about the expressed TEs.',
             tags = c('expression', 'rna-seq', 'TE', 'kimura'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin', 'blood'), sex = c('male')),
             script = 'rna_seq_male_supplemental_panel.R'
)

fig_index(plot = supplemental_pannel_male,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 17,
          height = 22,
          dpi = 300,
          format = 'pdf')

