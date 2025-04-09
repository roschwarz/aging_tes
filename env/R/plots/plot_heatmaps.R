#' ------------------------------------------------------------------------------
#' Heatmaps, z-score, tpms
#' ------------------------------------------------------------------------------
#' @param results A DESeq2 result object
#' @return A ggplot2 object
#' @export
heatmap_top50 <- function(df, tissue){
    
    ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points") 
    col_fun = colorRamp2(c(-2, 0, 2), c('#457b9d','white','#e63946'))
    
    row_text_size = 6
    
    row_split_vector <- data.frame(df) %>% 
        rownames_to_column("te_id") %>% 
        tidyr::separate(col = "te_id",
                        into = c("order", "super"),
                        sep = "[|]") %>%
        pull(order)
    
    row_annotation_text <-  data.frame(df) %>% 
        rownames_to_column("te_id") %>% 
        tidyr::separate(col = "te_id", into = c("order", "super", "fam"), sep = "[|]") %>% pull(fam) 
    
    fam_anno <- rowAnnotation(labels = anno_text(row_annotation_text, which = "row"))
    
    
    row_annotation <- rowAnnotation(labels = anno_text(row_annotation_text, 
                                                       which = "row",
                                                       gp = gpar(fontsize = rep(6, 50)))
    )
    
    Heatmap(
        df,
        col = col_fun ,
        column_title = tissue,
        column_title_gp = gpar(
            fill = tissue.color[tissue],
            color = 'white',
            border = 'black',
            fontsize = 8
        ),
        border = TRUE,
        column_gap = unit(0, 'mm'),
        row_names_gp = gpar(fontsize = row_text_size),
        column_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_title_gp = gpar(fontsize = 8),
        top_annotation = HeatmapAnnotation(age = anno_block(
            gp = gpar(fill = c('white', 'white')),
            labels = c('young', 'old'),
            labels_gp = gpar(col = 'black', fontsize = 8)
        )),
        right_annotation = row_annotation,
        column_km = 2,
        cluster_rows = FALSE,
        column_labels = rep('', ncol(df)),
        row_labels = rep('', 50),
        row_split = row_split_vector,
        heatmap_legend_param = list(
            title = "z-score of TPM",
            direction = 'horizontal',
            legend_width = unit(3, 'cm'),
            grid_height = unit(0.2, 'cm'),
            title_gp = gpar(fontsize = 6),
            labels_gp = gpar(fontsize = 6),
            title_position = "topcenter"
        )
    )
    
}