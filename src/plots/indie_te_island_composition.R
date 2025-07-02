# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()


indie_te_island_categorized <- blackRcloud::loadRdata(paste0(tables_and_co, "indie_te_island_categorized.Rdata"))

# ==============================================================================
# Proportion of TE islands
# ==============================================================================


indie_te_islands_composition <- indie_te_island_categorized$instance

indie_te_island_types_pls <- sapply(names(indie_te_islands_composition), simplify = F, function(x){
    
    df <- indie_te_islands_composition[[x]]
    
    df <- df %>% filter(!duplicated(te_island_id))
    
    df$island_type <- factor(df$island_type, levels = c('single', 'double', 'multiple'))
    
    print(paste(x, ":", nrow(df)))
    
    pl <- ggplot(df, aes(island_type)) +
        geom_bar() +
        labs(title = x) +
        theme(axis.text.x = element_text(angle = 90))
    
    return(pl)
    
})


gridExtra::grid.arrange(grobs = indie_te_island_types_pls, nrow = 1)

indie_single_te_island_pls <- sapply(names(indie_te_islands_composition), simplify = F, function(x){
    
    df <- indie_te_islands_composition[[x]]
    
    df <- df %>%
        filter(member == 1) %>%
        order.TEs('super_family', decreasing = F)
    
    pl <- ggplot(df, aes(super_family)) +
        geom_bar() +
        coord_flip()
    
    return(pl)
})

gridExtra::grid.arrange(grobs = indie_single_te_island_pls, nrow = 1)

# ==============================================================================
# Categorization of TE islands
# ==============================================================================

indie_te_island_super_family_composition <- indie_te_island_categorized$super_fam

ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points")    # setting to get the titels in the heatmaps centered from a horizontal perspective see: https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html?q=color%20title#heatmap-titles
col_fun = colorRamp2(c(0, 0.4), c('white','#e63946'))

pl_brain <-  Heatmap(indie_te_island_super_family_composition$brain,
                     col = col_fun,
                     border = TRUE,
                     column_gap = unit(0, 'mm'),
                     column_title_gp = gpar(fill = tissue.color['brain'], color = 'black', border = 'black', fontsize = 8), 
                     column_title = 'brain',
                     row_names_side = 'left',
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     show_heatmap_legend = FALSE,
                     cluster_columns = FALSE,
                     column_names_rot = 0,
                     column_names_centered = TRUE,
                     cluster_rows = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.2f", indie_te_island_super_family_composition$brain[i, j]), x, y, gp = gpar(fontsize = 6))
                     },
                     heatmap_legend_param = list(title = "Proportion of TE-families", 
                                                 direction = 'horizontal',
                                                 legend_width = unit(3, 'cm'))
)

pl_skin <-  Heatmap(indie_te_island_super_family_composition$skin,
                    col = col_fun,
                    border = TRUE,
                    column_title = "skin",
                    column_title_gp = gpar(fill = tissue.color['skin'], color = 'black', border = 'black', fontsize = 8), 
                    column_gap = unit(0, 'mm'),
                    row_names_side = 'left',
                    row_names_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 8),
                    show_heatmap_legend = TRUE,
                    cluster_columns = FALSE,
                    column_names_rot = 0,
                    column_names_centered = TRUE,
                    cluster_rows = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", indie_te_island_super_family_composition$skin[i, j]), x, y, gp = gpar(fontsize = 6))
                    },
                    heatmap_legend_param = list(title = "Proportion of TE-families", 
                                                direction = 'horizontal',
                                                legend_width = unit(3, 'cm'),
                                                grid_height = unit(0.2, 'cm'),
                                                title_gp = gpar(fontsize = 6),
                                                labels_gp = gpar(fontsize = 6),
                                                title_position = "topcenter"))

pl_blood <-  Heatmap(indie_te_island_super_family_composition$blood,
                     col = col_fun,
                     column_title = "blood",
                     column_title_gp = gpar(fill = tissue.color['blood'], color = 'white', border = 'black', fontsize = 8), 
                     border = TRUE,
                     column_gap = unit(0, 'mm'),
                     row_names_side = 'left',
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     show_heatmap_legend = FALSE,
                     cluster_columns = FALSE,
                     column_names_rot = 0,
                     column_names_centered = TRUE,
                     cluster_rows = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.2f", indie_te_island_super_family_composition$blood[i, j]), x, y, gp = gpar(fontsize = 6))
                     },
                     heatmap_legend_param = list(title = "Proportion of TE-families", 
                                                 direction = 'horizontal',
                                                 legend_width = unit(3, 'cm'),
                                                 legend_heigt = unit(0.25, 'cm'))
)


ht_list <- pl_brain + pl_skin + pl_blood

pdf(file = paste0(figure_dir, 'panel_3_indie_te_island_categorized_3.7x2.5_300.pdf'), width = 3.78, height = 2.8)
draw(ht_list, heatmap_legend_side = "bottom")
dev.off()

