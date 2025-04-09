# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_rna_seq_env()


# Using ComplexHeatmap package
te.mean.log2 <- deseq.te.merged %>% 
    filter(!is.na(padj), 
           order %in% orders.of.interest, 
           !grepl("[?]", super_family),
           super_family != "NA") %>% 
    mutate(super_family = case_when(super_family == "Alu" ~ "B1",
                                    TRUE ~ super_family)) %>% 
    group_by(tissue, 
             super_family) %>% 
    summarise(mean.log2FC = round(mean(log2FoldChange),2),
              super_family.members = n()) %>%
    filter(super_family.members > 10) %>% 
    dplyr::select(-super_family.members) %>% 
    spread(key = tissue, 
           value = mean.log2FC) %>%
    column_to_rownames('super_family') %>% 
    as.matrix()

te.mean.log2 <- te.mean.log2[,c('brain', 'skin', 'blood')]

col_fun = colorRamp2(c(-0.1, 0, 0.1), c(direction.color[['down']], 'white', direction.color[['up']]))

lgd = Legend(col_fun = col_fun, title = "log2(mean(fold change))", direction = "horizontal")


hm <- Heatmap(te.mean.log2,
              column_split = c(1,2,3),
              col = col_fun,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(tissue.color[['brain']],
                                                                                     tissue.color[['skin']],
                                                                                     tissue.color[['blood']])),
                                                                  labels = c("brain", "skin", "blood"), 
                                                                  labels_gp = gpar(col = "white", fontsize = 10))),
              border = TRUE,
              column_title = NULL, # c('','',''),
              column_gap = unit(0, 'mm'),
              row_names_side = 'left',
              row_names_gp = gpar(fontsize = 10),
              column_labels = c('','',''),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              heatmap_legend_param = list(title = "log2(mean(fold change))", 
                                          direction = 'horizontal',
                                          legend_width = unit(5, 'cm'))
)


pl <- draw(hm, heatmap_legend_side = "bottom")

pdf(file = paste0(figure_dir, '04_log2fc_heat_superFam_9.6x13_300.pdf'), width = 3.78, height = 5.12)
    draw(hm, heatmap_legend_side = "bottom")
dev.off()
