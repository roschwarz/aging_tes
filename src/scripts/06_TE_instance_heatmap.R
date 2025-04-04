# --------------------------------- Notes--- -----------------------------------
# 
# This script is used for the RNA-Seq quantification
#
# Output data (../data/rna; ../tables):
#

# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    source('./01_load_environment.R')
} else if (!ENVIRONMENT_LOADED) {
    
    source('./01_load_environment.R')
    
}

if (!file.exists(paste0(data_dir, "rna_RS_tpms_TE_instances_SalmonTE.Rdata"))) {
    
    dds.own.rna <- loadRdata(paste0(data_dir, "rna_RS_dds_TE_instances_SalmonTE.Rdata"))
    
    tpm.own.rna.exp <- sapply(tissues, simplify = F, function(x){
        
        counts <- counts(dds.own.rna[[x]])
        
        te.length <- data.frame(te_id = row.names(counts))
        
        te.length <- splitTEID(te.length, "te_id")
        
        te.length$length <- as.numeric(te.length$end) - as.numeric(te.length$start)
        
        te.length <- te.length %>% 
            dplyr::select(te_id, length)
        
        norm <- normalizeCountMatrix(counts, te.length)
        
        names(norm) <- str_replace(names(norm), "TPM." , "age_")
        
        norm <- rownames_to_column(norm, var = "te_id")
        
        norm <- norm %>% 
            gather(key = 'sample', value = "TPM", names(norm)[2:length(norm)]) %>% 
            splitTEID("te_id")
        
        
        
    })
    
    save(tpm.own.rna.exp, file = paste0(data_dir, "rna_RS_tpms_TE_instances_SalmonTE.Rdata"))
    
}else{
    
    tpm.own.rna.exp <- loadRdata(paste0(data_dir, "rna_RS_tpms_TE_instances_SalmonTE.Rdata"))
    
}

deseq.own.rna.merged <- read.csv(paste0(table_dir,
                                        "02_deseq_results_te_instances.csv"))

top_50 <- sapply(tissues, simplify = F, function(x){
    
    norm <-  tpm.own.rna.exp[[x]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
                                             !grepl("[?]", super_family),
                                             super_family != "NA") %>% 
        
        filter(te_id %in% (deseq.own.rna.merged %>% 
                                             filter(tissue == x) %>%
                                             dplyr::slice_min(order_by = padj, n = 50) %>% 
                                             dplyr::slice_max(order_by = baseMean, n = 50) %>% 
                                             dplyr::pull(te_id))) #%>%  
       
    
    norm <- norm %>% 
        separate('sample', c("age.name", "age", "sample.id", "tissue"), sep = "[_]", remove = F) %>% 
        dplyr::select(-age.name)
    
    norm <- norm %>% 
        mutate(te.annotation = paste0(order, '|', super_family, '|', family, '|', start)) %>% 
        dplyr::select(te.annotation, sample, TPM) %>% 
        spread(sample, TPM) %>%
        column_to_rownames(var = 'te.annotation')
        
    # the order is different in skin because of the naming
    if (x == 'skin') {
        norm <- norm[,rev(names(norm))]
    }
    
    norm <- as.matrix(norm)
    
    norm <- t(apply(norm, 1, function(x) calc_z_score(x, capped = 2)))
    
    # removes rows full of NA --> Why such rows are there?
    norm <- norm[rowSums(is.na(norm)) != ncol(norm), ]
    
    
})





################### Next to doooo ################

ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points") 
col_fun = colorRamp2(c(-2, 0, 2), c('#457b9d','white','#e63946'))

row_text_size = 6

brain_row_split <- data.frame(top_50$brain) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super"), sep = "[|]") %>% pull(order)

brain_row_annotation_text <-  data.frame(top_50$brain) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super", "fam"), sep = "[|]") %>% pull(fam)

brain_fam_anno <- rowAnnotation(labels = anno_text(brain_row_annotation_text, which = "row"))


brain_row_annotation <- rowAnnotation(labels = anno_text(brain_row_annotation_text, 
                                                      which = "row",
                                                      gp = gpar(fontsize = rep(6, 50)))
                                      )

brain_heatMap <- Heatmap(top_50$brain,
        col = col_fun ,
        column_title = 'brain',
        column_title_gp = gpar(fill = tissue.color['brain'], color = 'white', border = 'black', fontsize = 8),
        border = TRUE,
        column_gap = unit(0, 'mm'),
        row_names_gp = gpar(fontsize = row_text_size),
        column_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_title_gp = gpar(fontsize = 8),
        top_annotation = HeatmapAnnotation(age = anno_block(gp = gpar(fill = c('white', 'white')),
                                                            labels = c('young', 'old'),
                                                            labels_gp = gpar(col = 'black', fontsize = 8))),
        right_annotation = brain_row_annotation,
        column_km = 2,
        cluster_rows = FALSE,
        column_labels = rep('', 10),
        row_labels = rep('', 50),
        row_split = brain_row_split,
        heatmap_legend_param = list(title = "z-score of TPM", 
                                    direction = 'horizontal',
                                    legend_width = unit(3, 'cm'),
                                    grid_height = unit(0.2, 'cm'),
                                    title_gp = gpar(fontsize = 6),
                                    labels_gp = gpar(fontsize = 6),
                                    title_position = "topcenter"))

brain_heatMap

pdf(file = paste0(figure_dir, '06_te_instance_heat_brain_2.6x5.12.pdf'), width = 2.6, height = 5.12)
draw(brain_heatMap)
dev.off()


skin_row_split <- data.frame(top_50$skin) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super"), sep = "[|]") %>% pull(order)

skin_row_annotation_text <-  data.frame(top_50$skin) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super", "fam"), sep = "[|]") %>% pull(fam)

skin_fam_anno <- rowAnnotation(labels = anno_text(skin_row_annotation_text, which = "row"))


skin_row_annotation <- rowAnnotation(labels = anno_text(skin_row_annotation_text, 
                                                         which = "row",
                                                         gp = gpar(fontsize = rep(6, 50)))
)

skin_HeatMap <- Heatmap(top_50$skin,
                        col = col_fun ,
                        column_title = 'skin',
                        column_title_gp = gpar(fill = tissue.color['skin'], color = 'white', border = 'black', fontsize = 8),
                        border = TRUE,
                        column_gap = unit(0, 'mm'),
                        row_names_gp = gpar(fontsize = row_text_size),
                        column_names_gp = gpar(fontsize = 8),
                        show_heatmap_legend = FALSE,
                        cluster_columns = FALSE,
                        column_names_rot = 0,
                        column_names_centered = TRUE,
                        row_title_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(age = anno_block(gp = gpar(fill = c('white', 'white')),
                                                                            labels = c('young', 'old'),
                                                                            labels_gp = gpar(col = 'black', fontsize = 8))),
                        right_annotation = skin_row_annotation,
                        column_km = 2,
                        cluster_rows = FALSE,
                        column_labels = rep('', 10),
                        row_labels = rep('', 50),
                        row_split = skin_row_split,
                        heatmap_legend_param = list(title = "z-score of TPM", 
                                                    direction = 'horizontal',
                                                    legend_width = unit(3, 'cm'),
                                                    grid_height = unit(0.2, 'cm'),
                                                    title_gp = gpar(fontsize = 6),
                                                    labels_gp = gpar(fontsize = 6),
                                                    title_position = "topcenter"))


pdf(file = paste0(figure_dir, '06_te_instance_heat_skin_2.6x5.12.pdf'), width = 2.6, height = 5.12)
draw(skin_HeatMap)
dev.off()



blood_row_split <- data.frame(top_50$blood) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super"), sep = "[|]") %>% pull(order)

blood_row_split <- data.frame(top_50$blood) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super"), sep = "[|]") %>% pull(order)

blood_row_annotation_text <-  data.frame(top_50$blood) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super", "fam"), sep = "[|]") %>% pull(fam)

blood_fam_anno <- rowAnnotation(labels = anno_text(blood_row_annotation_text, which = "row"))

blood_row_annotation <- rowAnnotation(labels = anno_text(blood_row_annotation_text, 
                                                        which = "row",
                                                        gp = gpar(fontsize = rep(6, 50)))
)

blood_HeatMap <- Heatmap(top_50$blood,
                        col = col_fun ,
                        column_title = 'blood',
                        column_title_gp = gpar(fill = tissue.color['blood'], color = 'white', border = 'black', fontsize = 8),
                        border = TRUE,
                        column_gap = unit(0, 'mm'),
                        row_names_gp = gpar(fontsize = row_text_size),
                        column_names_gp = gpar(fontsize = 8),
                        show_heatmap_legend = FALSE,
                        cluster_columns = FALSE,
                        column_names_rot = 0,
                        column_names_centered = TRUE,
                        row_title_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(age = anno_block(gp = gpar(fill = c('white', 'white')),
                                                                            labels = c('young', 'old'),
                                                                            labels_gp = gpar(col = 'black', fontsize = 8))),
                        right_annotation = blood_row_annotation,
                        column_km = 2,
                        cluster_rows = FALSE,
                        column_labels = rep('', 9),
                        row_labels = rep('', 50),
                        row_split = blood_row_split,
                        heatmap_legend_param = list(title = "z-score of TPM", 
                                                    direction = 'horizontal',
                                                    legend_width = unit(3, 'cm'),
                                                    grid_height = unit(0.2, 'cm'),
                                                    title_gp = gpar(fontsize = 6),
                                                    labels_gp = gpar(fontsize = 6),
                                                    title_position = "topcenter"))




pdf(file = paste0(figure_dir, '06_te_instance_heat_blood_2.6x5.12.pdf'), width = 2.6, height = 5.12)
draw(blood_HeatMap)
dev.off()

# with legend to get one legend

blood_HeatMap <- Heatmap(top_50$blood,
                        col = col_fun ,
                        column_title = 'blood',
                        column_title_gp = gpar(fill = tissue.color['blood'], color = 'white', border = 'black', fontsize = 8),
                        border = TRUE,
                        column_gap = unit(0, 'mm'),
                        row_names_gp = gpar(fontsize = row_text_size),
                        column_names_gp = gpar(fontsize = 8),
                        #show_heatmap_legend = FALSE,
                        cluster_columns = FALSE,
                        column_names_rot = 0,
                        column_names_centered = TRUE,
                        row_title_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(age = anno_block(gp = gpar(fill = c('white', 'white')),
                                                                            labels = c('young', 'old'),
                                                                            labels_gp = gpar(col = 'black', fontsize = 8))),
                        right_annotation = blood_row_annotation,
                        column_km = 2,
                        cluster_rows = FALSE,
                        column_labels = rep('', 9),
                        row_labels = rep('', 50),
                        row_split = blood_row_split,
                        heatmap_legend_param = list(title = "z-score of TPM", 
                                                    direction = 'horizontal',
                                                    legend_width = unit(3, 'cm'),
                                                    grid_height = unit(0.2, 'cm'),
                                                    title_gp = gpar(fontsize = 6),
                                                    labels_gp = gpar(fontsize = 6),
                                                    title_position = "topcenter"))


pdf(file = paste0(figure_dir, '06_te_instance_heat_blood_LEGEND_2.6x5.12.pdf'), width = 2.6, height = 5.12)
draw(blood_HeatMap, heatmap_legend_side = "bottom")
dev.off()


#########
# RESIS #
#########


pdf(file = './manuscripts/nature_aging/resis/Figure1/1D_heatmap_brain.pdf', width = 2.6, height = 5.12)
draw(brain_heatMap)
dev.off()

write.table(as.data.frame(top_50$brain) %>% rownames_to_column(var = "te_id"), 
            file = './manuscripts/nature_aging/resis/Figure1/1D_heatmap_brain.csv', 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

pdf(file = './manuscripts/nature_aging/resis/Figure1/1D_heatmap_skin.pdf', width = 2.6, height = 5.12)
draw(skin_HeatMap)
dev.off()

write.table(as.data.frame(top_50$skin) %>% rownames_to_column(var = "te_id"), 
            file = './manuscripts/nature_aging/resis/Figure1/1D_heatmap_skin.csv', 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

pdf(file = './manuscripts/nature_aging/resis/Figure1/1D_heatmap_blood.pdf', width = 2.6, height = 5.12)
draw(blood_HeatMap)
dev.off()

write.table(as.data.frame(top_50$blood) %>% rownames_to_column(var = "te_id"), 
            file = './manuscripts/nature_aging/resis/Figure1/1D_heatmap_blood.csv', 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


##### Next merge the plots


ggarrange(brain_heatMap, skin_HeatMap, brain_heatMap)


blood_HeatMap <- Heatmap(top_50$skin,
                        col = col_fun ,
                        column_title = 'skin',
                        column_title_gp = gpar(fill = tissue.color['skin'], color = 'white', border = 'black', fontsize = 8),
                        border = TRUE,
                        column_gap = unit(0, 'mm'),
                        row_names_gp = gpar(fontsize = row_text_size),
                        column_names_gp = gpar(fontsize = 8),
                        show_heatmap_legend = FALSE,
                        cluster_columns = FALSE,
                        column_names_rot = 0,
                        column_names_centered = TRUE,
                        row_title_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(age = anno_block(gp = gpar(fill = c('white', 'white')),
                                                                            labels = c('young', 'old'),
                                                                            labels_gp = gpar(col = 'black', fontsize = 8))),
                        right_annotation = skin_row_annotation,
                        column_km = 2,
                        cluster_rows = FALSE,
                        column_labels = rep('', 10),
                        row_labels = rep('', 50),
                        row_split = skin_row_split,
                        heatmap_legend_param = list(title = "z-score of TPM", 
                                                    direction = 'horizontal',
                                                    legend_width = unit(3, 'cm'),
                                                    grid_height = unit(0.2, 'cm'),
                                                    title_gp = gpar(fontsize = 6),
                                                    labels_gp = gpar(fontsize = 6),
                                                    title_position = "topcenter"))


#ht_list <- brain_heatMap + skin_HeatMap + blood_HeatMap
draw(brain_heatMap, heatmap_legend_side = "bottom")
draw(skin_HeatMap, heatmap_legend_side = "bottom")
draw(blood_HeatMap, heatmap_legend_side = "bottom")


#pdf(file = paste0(figure_dir, '09_superFam_composition_te_island_3.7x2.5_300.pdf'), width = 3.78, height = 2.8)
draw(ht_list, heatmap_legend_side = "bottom")
dev.off()



click <- data.frame(top_50$skin) %>% 
    rownames_to_column("te_id") %>% 
    tidyr::separate(col = "te_id", into = c("order", "super"), sep = "[|]") %>% pull(order)

Heatmap(auto.region.super.composition$blood,
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
            grid.text(sprintf("%.2f", auto.region.super.composition$blood[i, j]), x, y, gp = gpar(fontsize = 6))
        },
        heatmap_legend_param = list(title = "Proportion of TE-families", 
                                    direction = 'horizontal',
                                    legend_width = unit(3, 'cm'),
                                    legend_heigt = unit(0.25, 'cm'))
)


colnames(top_50$brain) <- NULL

pheatmap(top_50$brain, 
         border_color = NA,
         clustering_method = 'average',
         fontsize = 8,
         cluster_rows = F,
         cluster_cols = F,
         
         
         )


         clustering_method = 'average',
         #annotation_col = col_annotation,
         fontsize = 8,
         fontfamily = 'Arial',
         cluster_rows = F,
         cluster_cols = F,
         # breaks = seq(-2, 2, length.out = 100),
         # color = colorRampPalette(c('#457b9d', 'white','#e63946'))(101),
         cutree_rows = 2,
         cutree_cols = 2,
         treeheight_row = 0, # removes dendrogram
         treeheight_col = 0 
)




write.table(as.data.frame(top_50$skin) %>% rownames_to_column(var = "te_id"), 
            file = '../../submission/Resis/figure1/figure_1d_skin.csv', 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


colnames(top_50$skin) <- NULL

pheatmap(top_50$skin, 
         border_color = NA,
         clustering_method = 'average',
         annotation_names_col = F,
         fontsize = 8,
         fontfamily = 'Arial',
         cluster_rows = F,
         cluster_cols = F,
         breaks = seq(-2, 2, length.out = 100),
         color = colorRampPalette(c('#457b9d', 'white','#e63946'))(101),
         cutree_rows = 2,
         cutree_cols = 2,
         treeheight_row = 0, # removes dendrogram
         treeheight_col = 0 
)

write.table(as.data.frame(top_50$blood) %>% rownames_to_column(var = "te_id"), 
            file = '../../submission/Resis/figure1/figure_1d_blood.csv', 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')



colnames(top_50$blood) <- NULL

pheatmap(top_50$blood, 
         border_color = NA,
         clustering_method = 'average',
         #annotation_col = col_annotation,
         fontsize = 8,
         fontfamily = 'Arial',
         cluster_rows = F,
         cluster_cols = F,
         breaks = seq(-2, 2, length.out = 100),
         color = colorRampPalette(c('#457b9d', 'white','#e63946'))(101),
         cutree_rows = 2,
         cutree_cols = 2,
         treeheight_row = 0, # removes dendrogram
         treeheight_col = 0 
)
