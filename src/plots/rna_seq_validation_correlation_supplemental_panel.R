if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_public_data_env()
aging_tes::load_plotting_env()
aging_tes::load_annotations()
load_te_annotation()

color_set <- unlist(tissue.color[c('Gastrocnemius muscle', "liver", 'White adipose tissue')])
strip_text_color <- c( "#ffffff","#000000", "#000000", "#ffffff", "#ffffff")

# ------------------------------------------------------------------------------
# Volcano plots for public data
# ------------------------------------------------------------------------------

# all transposable elements

deseq_results_public <- fread(paste0(rna_seq_deseq_dir, "deseq_results_te_instances_public.csv"))


public_volcano_all_TEs <- volcanoPlot(deseq_results_public, FDR = 0.05, c("sex", "tissue")) +
    theme_rob(base_size = 8) +
    theme(legend.position = 'None')

public_volcano_all_TEs <- color_strips(public_volcano_all_TEs, 
                               bg_cols = c(color_set, "#136f63", "#ffba08"), 
                               text_cols = strip_text_color)



# young transposable elements (Kimura distance <= 5)

young_TEs <- deseq_results_public %>% 
    filter(Kimura <= 5)

public_volcano_young_TEs <- volcanoPlot(young_TEs, FDR = 0.05, c("sex", "tissue")) +
    theme_rob(base_size = 8) +
    theme(legend.position = 'None')

public_volcano_young_TEs <- color_strips(public_volcano_young_TEs, 
                               bg_cols = c(color_set, "#136f63", "#ffba08"), 
                               text_cols = strip_text_color)



public_volcano_plots <- gridExtra::grid.arrange(public_volcano_all_TEs, public_volcano_young_TEs, nrow = 2)

###### Correlation #####

correlation_df <- data.table::fread(paste0(table_dir,
                                             "te_host_gene_correlation_data.csv"))

correlation_df$tissue <- factor(correlation_df$tissue, levels = tissues)

# Calculate some value beforhand to inprove the perfomance during the creation of the plot

# Visible area of the final plot
xlim <- c(-10, 10)
ylim <- c(-10, 10)

# filter for visible are of the final plot
df_plot <- correlation_df[
    gene_log2FoldChange >= xlim[1] & gene_log2FoldChange <= xlim[2] &
        te_log2FoldChange   >= ylim[1] & te_log2FoldChange   <= ylim[2]
]

# Spearman-correlation for each facet
cor_dt <- df_plot[, {
    x <- gene_log2FoldChange
    y <- te_log2FoldChange
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 3) {
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        list(rho = unname(ct$estimate), p = ct$p.value)
    } else list(rho = NA_real_, p = NA_real_)
}, by = .(position, tissue)]


cor_dt[, label := sprintf("rho = %.2f, p =%.2g", rho, p)]
cor_dt[, `:=`(label_x = -9, label_y = 7)]

# Count per facet

df_plot[, quadrant :=
            fifelse(gene_log2FoldChange >= 0 & te_log2FoldChange >= 0, "Q1",
            fifelse(gene_log2FoldChange < 0 & te_log2FoldChange >= 0, "Q2",
            fifelse(gene_log2FoldChange < 0 & te_log2FoldChange < 0, "Q3", "Q4")))]

q_counts <- df_plot[, .N, by = .(position, tissue, quadrant)]

# fix position into the cornors of each facet
quad_pos <- data.table(
    quadrant = c("Q1","Q2","Q3","Q4"),
    x = c(7, -9, -9, 7),
    y = c(7, 7, -9, -9)
)

q_ann <- q_counts[quad_pos, on = 'quadrant']
q_ann[, label := sprintf("n = %s", format(N, big.mark = ",", decimal.mark = "."))]


pl <- ggplot(df_plot, aes(gene_log2FoldChange, te_log2FoldChange)) +
    geom_point_rast(alpha = 0.05, color = '#192E37', raster.dpi = 150) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    facet_grid(position ~ tissue) +
    geom_text(data = cor_dt,
              aes(x = label_x, y = label_y, label = label),
              inherit.aes = FALSE, size = 6/.pt) +
    geom_text(data = q_ann,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 6/.pt) +
    labs(y = expression(paste(log[2], "(fold ", change[TE], ")")),
         x = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    theme_rob(base_size = 8)

pl <- color_strips(pl,
                   bg_cols = c(tissue.color[2:4], "#ffffff","#ffffff"), 
                   text_cols = c( "#ffffff", "#000000","#ffffff", "#000000", "#000000"))

pl <- grid.grabExpr(grid.draw(pl)) %>% 
    ggplotify::as.ggplot()

suplemental_plot <- gridExtra::grid.arrange(public_volcano_plots, pl, nrow = 2, heights = c(2,1))

meta <- list(name = 'suplemental_validation_corr',
             description = 'Supplemental figure: Public RNA-seq volcano plots and TE-host gene correlation',
             tags = c('rna-seq', 'volcano plot', 'correlation'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin'), sex = c('female')),
             script = 'rna_seq_validation_correlation_supplemental_panel.R'
)

fig_index(plot = suplemental_plot,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 18,
          height = 25,
          dpi = 300,
          format = 'pdf')

