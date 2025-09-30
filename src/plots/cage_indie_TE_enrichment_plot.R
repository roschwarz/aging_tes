# ------------------------------------------------------------------------------
# Enrichment analysis for transposable elements with their own TSS
#
# Prerequisite is that cage_store_annotation was applied. 
#
# ------------------------------------------------------------------------------


# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_cage_seq_env()
aging_tes::load_annotations()
aging_tes::load_plotting_env()

load_cage_peak_annotation()
load_te_ranges()
load_te_annotation()

# ------------------------------------------------------------------------------
# Identify transposable elements with their own TSS
# ------------------------------------------------------------------------------

te_cages <- sapply(names(cageRanges), simplify = F, function(x) {
    
    df <- process_overlapping(
        query = teRanges,
        subject = cageRanges[[x]],
        col_query_to_subj = "te.id"
    )
    
    return(df)
    
})

# ------------------------------------------------------------------------------
# Calculate enrichment of TE instances with a CTSS within their super families 
# ------------------------------------------------------------------------------


te_super_enrichment <- do.call('rbind', sapply(names(te_cages), simplify = F, function(x) {
    df <- data.frame(te_id = names(unique(te_cages[[x]]$query_with_hit)),
                     row.names = NULL)
    df <- splitTEID(df, 'te_id')
    df <- blackRcloud::binoRich(
        df,
        te_annotation,
        c(target = 'super_family', background = 'super_family'),
        FDR.cap = 1e-10
    )
    df <- df %>%
        filter(n.target >= 10) %>%
        filter(!grepl("[?]", category))
    
    df$tissue <- x
    
    return(df)
    
}))

# ------------------------------------------------------------------------------
# Clean up the data frame 
#
# - rename of Alu to B1
# - sort the data frame for the plot order
# ------------------------------------------------------------------------------

te_super_enrichment <- te_super_enrichment %>% 
    filter(category != 'NA') %>% 
    mutate(category = case_when(category == 'Alu' ~ 'B1',
                                TRUE ~ category))

te_super_enrichment$tissue <- 
    factor(te_super_enrichment$tissue, levels = c('brain', 'skin', 'blood'))

cat_Sort <- te_super_enrichment %>% 
    group_by(category) %>% 
    summarise(meanOdds = sum(ratio)) %>% 
    arrange(desc(meanOdds))

te_super_enrichment$category <-
    factor(te_super_enrichment$category, levels = cat_Sort$category)

# ------------------------------------------------------------------------------
# Create Plot
# ------------------------------------------------------------------------------

x_max <- max(te_super_enrichment$ratio) + 0.2

label_pos <- te_super_enrichment %>% 
    group_by(tissue) %>% 
    dplyr::count() %>% pull(n) %>% max() + 2

cage_enrichment_pl <- ggplot(te_super_enrichment, aes(ratio, category, size = log10.padj)) +
    geom_point(aes(fill = tissue), shape = 21, color = 'black') +
    scale_fill_manual(values = tissue.color) +
    scale_size(range = c(1, 5), name = expression(paste(log[10], "(FDR)"))) +
    annotate("text", x = -1, y = 0, label = "depleted", vjust = -1.5, color = 'white') + # trick to have a gap at the bottom
    geom_hline(yintercept = label_pos - 1, linewidth = 0.4) + 
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    annotate("text", x = x_max/2, y = label_pos, label = "enriched", vjust = 1.2, size = 8/.pt) +
    annotate("text", x = -x_max/2, y = label_pos, label = "depleted", vjust = 1.2, size = 8/.pt) +
    theme_rob(base_size = 12, base_family = 'arial') +
    xlim(c(-x_max, x_max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    theme_rob(base_size = 8, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.position = c(0.84, 0.65),
          legend.key.size = unit(0, "lines"),
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, size = 0.4),
          legend.box = "vertical",
          legend.box.just = 'right',) +
    guides(size = guide_legend(direction = 'vertical'))

cage_enrichment_pl

# Save figure with metadata to index
meta <- list(name = 'cage_TE_enrichment',
             description = 'Enrichment of CAGE peaks within transposable element superfamilies in brain, skin and blood',
             tags = c('CAGE-Seq', 'TEs', 'Enrichment'),
             parameters = list(n_fam_count= ">10", tissues = c('brain', 'skin', 'blood')),
             script = 'cage_indie_TE_enrichment_plot.R'
)

fig_index(plot = cage_enrichment_pl,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 7,
          height = 8,
          dpi = 300,
          format = 'pdf')
