# The script loads the .csv-files created by calc_z-scores_pcdhb15_locus.R and makes a heatmap out of it.
#
# INPUT:
# ./tables/pcdhb15_locus_heatmap_zscore_vst.csv
# ./tabels/pcdhb15_locus_heatmap_log2fc_fdr.csv
#
# OUTPUT:
# ./figures/heatmap_pcdhb15_locus_vst.pdf

# Author: Robert Schwarz
# Created: 2025-11-07


# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

vst_counts_pro_locus <- read.csv(paste0(table_dir, "pcdhb15_locus_heatmap_zscore_vst.csv"))
df <- read.csv(paste0(table_dir, "pcdhb15_locus_heatmap_log2fc_fdr.csv"))
# ------------------------------------------------------------------------------------------------------------
# Create heat map
# ------------------------------------------------------------------------------------------------------------

col_fun_l2fc = colorRamp2(c(-2,-1,0,1,2), c('#264653','#2a9d8f', 'white','#f4a261', '#e76f51'))
col_fun_fdr = colorRamp2(c(0.000001,0.00001,0.0001, 0.001, 0.01, 0.1, 1), 
                         c('#AC0A0A','#B82C2C','#C44F4F','#D17272','#DD9595', '#E9B8B8', '#F6DBDB'))

col_fun_fdr = colorRamp2(c(0.00000001,0.0000001,0.000001,
                           0.00001,0.0001,0.001, 
                           0.05, 0.1, 0.5,
                           1), 
                         c("#ea698b","#d55d92","#c05299","#ac46a1","#973aa8","#822faf","#6d23b6","#6411ad","#571089","#47126b"))

row_legend = list(direction = "horizontal",
                  legend_width = unit(2, 'cm'),
                  grid_height = unit(0.1, 'cm'),
                  title_gp = gpar(fontsize = 6),
                  labels_gp = gpar(fontsize = 6),
                  title_position = "topcenter")

row_ha = rowAnnotation(L2FC = df$log2FoldChange,
                       FDR = df$padj,
                       col = list(L2FC = col_fun_l2fc,
                                  FDR = col_fun_fdr),
                       simple_anno_size = unit(0.3, 'cm'),
                       annotation_name_gp = gpar(fontsize = c(6,6)),
                       annotation_name_side = c('top', 'top'),
                       annotation_legend_param = list(L2FC = row_legend, 
                                                      FDR = row_legend))



ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points") 
col_fun = colorRamp2(c(-1.5, 0, 1.5), c('#457b9d','white','#e63946'))

row_text_size = 8

heatmap_pcdhb15_locus <- Heatmap(vst_counts_pro_locus,
                                 col = col_fun ,
                                 column_title = "Protocadherin beta cluster",
                                 column_title_gp = gpar(color = 'white', border = 'black', fontsize = 8),
                                 border = TRUE,
                                 column_gap = unit(0, 'mm'),
                                 row_names_gp = gpar(fontsize = row_text_size),
                                 column_names_gp = gpar(fontsize = 8),
                                 show_heatmap_legend = TRUE,
                                 row_names_side = "left",
                                 cluster_columns = FALSE,
                                 column_names_rot = 0,
                                 column_names_centered = TRUE,
                                 row_title_gp = gpar(fontsize = 8),
                                 top_annotation = HeatmapAnnotation(age = anno_block(gp = gpar(fill = c('white', 'white')),
                                                                                     labels = c('young', 'old'),
                                                                                     labels_gp = gpar(col = 'black', fontsize = 8))),
                                 right_annotation = row_ha,
                                 column_km = 2,
                                 cluster_rows = FALSE,
                                 column_labels = rep('', 10),
                                 heatmap_legend_param = list(title = "z-score of TPM", 
                                                             direction = 'horizontal',
                                                             legend_width = unit(2, 'cm'),
                                                             grid_height = unit(0.1, 'cm'),
                                                             title_gp = gpar(fontsize = 6),
                                                             labels_gp = gpar(fontsize = 6),
                                                             title_position = "topcenter"))

meta <- list(name = 'heatmap_pcdhb15_locus_vst',
             description = 'Heatmap of the z-score of the vst of a part of the protocadherin beta cluster in brain',
             tags = c('expression', 'rna-seq', 'TE island', 'z-score', 'vst'),
             parameters = list(tissue = 'brain', metric = 'z-score of vst'),
             script = 'heatmap_pcdhb15_locus.R')

fig_index(plot = draw(heatmap_pcdhb15_locus, heatmap_legend_side = 'bottom'),
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 8,
          height = 11)