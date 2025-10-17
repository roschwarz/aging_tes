# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_rna_seq_env()


deseq.te.merged <- fread(paste0(table_dir, deseq_results_te_csv))


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

# === horizontal ===

hm_horizontal <-  Heatmap(t(te.mean.log2),
                          width = ncol(t(te.mean.log2)) * unit(3, "mm"),
                          height = nrow(t(te.mean.log2)) * unit(3, "mm"),
                          row_names_rot = 0,
                          col = col_fun,
                          column_names_rot = 90,
                          column_names_side = "top",
                          row_names_side = "left",
                          cluster_columns = FALSE,
                          cluster_rows = FALSE,
                          #row_split = sex_split,
                          # left_annotation = rowAnnotation(
                          #     block = anno_block(
                          #         gp = gpar(fill = sex_colors),
                          #         labels = c("female", "male"),
                          #         labels_rot = 0,
                          #         labels_gp = gpar(col = c("black", "white"), fontsize = 8)
                          #     )
                          # ),
                          heatmap_legend_param = list(title = "log2(mean(fold change))", 
                                                      direction = 'vertical',
                                                      grid_width = unit(0.2, "cm"),
                                                      legend_width = unit(3.5, 'cm'),
                                                      labels_gp = gpar(fontsize = 8),
                                                      title_gp = gpar(fontsize = 8)),
                          border = TRUE,
                          row_title = NULL,
                          column_names_gp = gpar(fontsize = 8),
                          row_names_gp = gpar(fontsize = 8)
                          
)

draw(hm_horizontal, heatmap_legend_side = "right")

meta <- list(name = 'te_super_families_log2fc_heatmap_flipped',
             description = 'Heatmap of log2 mean fold change of transposable element super families in male in a horizontal version',
             tags = c('expression', 'rna-seq', 'super_family'),
             parameters = list(tissues = c('brain', 'skin', 'blood')),
             script = 'rna_seq_family_heatmap.R'
)

fig_index(plot = hm_horizontal,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 15,
          height = 3,
          dpi = 300,
          format = 'pdf')



#### Female ####

deseq.te.merged.female <- fread(paste0(table_dir, "02_deseq_results_te_female_instances.csv"))


# Using ComplexHeatmap package
te.mean.log2_female <- deseq.te.merged.female %>% 
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

te.mean.log2_female <- te.mean.log2_female[,c('brain', 'skin')]

col_fun = colorRamp2(c(-0.1, 0, 0.1), c(direction.color[['down']], 'white', direction.color[['up']]))

lgd = Legend(col_fun = col_fun, title = "log2(mean(fold change))", direction = "horizontal")


hm_female <- Heatmap(te.mean.log2_female,
              column_split = c(1,2),
              col = col_fun,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(tissue.color[['brain']],
                                                                                     tissue.color[['skin']])),
                                                                  labels = c("brain", "skin"), 
                                                                  labels_gp = gpar(col = "white", fontsize = 10))),
              border = TRUE,
              column_title = NULL, # c('','',''),
              column_gap = unit(0, 'mm'),
              row_names_side = 'left',
              row_names_gp = gpar(fontsize = 10),
              column_labels = c('',''),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              heatmap_legend_param = list(title = "log2(mean(fold change))", 
                                          direction = 'horizontal',
                                          legend_width = unit(5, 'cm'))
)


pl_female <- draw(hm_female, heatmap_legend_side = "bottom")


# === Combine male and female in one plot ===

female_te_mean_log2 <- data.frame(te.mean.log2_female) %>% 
    dplyr::rename('f_brain' = brain,
                  'f_skin' = skin) %>% 
    rownames_to_column('super_family')

male_te_mean_log2 <- data.frame(te.mean.log2) %>% 
    dplyr::rename('m_brain' = brain,
                  'm_skin' = skin,
                  'm_blood' = blood) %>% 
    rownames_to_column('super_family')

te_mean_log2 <- merge(male_te_mean_log2, female_te_mean_log2, by = 'super_family', all = TRUE) %>% 
    column_to_rownames('super_family') %>% 
    as.matrix()


colMeans(te_mean_log2, na.rm = TRUE)

col_fun = colorRamp2(c(-0.1, 0, 0.1), c(direction.color[['down']], 'white', direction.color[['up']]))
lgd = Legend(col_fun = col_fun, title = "log2(mean(fold change))", direction = "horizontal")

# === Vertical ===

hm_merged <- Heatmap(te_mean_log2,
              column_split = c(1,2,3,4,5),
              col = col_fun,
              top_annotation = (
                  HeatmapAnnotation(
                      Sex = anno_block(
                          gp = gpar(fill = c("#cad2c5", "#cad2c5", "#cad2c5", "#", "#0f4c5c")),
                          labels = c("male", "male", "male", "female", "female"),
                          labels_gp = gpar(col = "white", fontsize = 12)
                      ),
                      tissue = anno_block(
                              gp = gpar(fill = c(
                                  tissue.color[['brain']],
                                  tissue.color[['skin']],
                                  tissue.color[['blood']],
                                  tissue.color[['brain']],
                                  tissue.color[['skin']]
                              )),
                              labels = c("brain", "skin", "blood", "brain", "skin"),
                              labels_gp = gpar(col = "white", fontsize = 10)
                          )
                      )),
              border = TRUE,
              column_title = NULL, # c('','',''),
              column_gap = unit(0, 'mm'),
              row_names_side = 'left',
              row_names_gp = gpar(fontsize = 10),
              column_labels = c('','','','',''),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              heatmap_legend_param = list(title = "log2(mean(fold change))", 
                                          direction = 'horizontal',
                                          legend_width = unit(5, 'cm'))
)

pl_merged <- draw(hm_merged, heatmap_legend_side = "bottom")




# Save figure with metadata to index
meta <- list(name = 'te_super_families_log2fc_heatmap',
             description = 'Heatmap of log2 mean fold change of transposable element super families in male and female',
             tags = c('expression', 'rna-seq', 'super_family'),
             parameters = list(tissues = c('brain', 'skin', 'blood'), sex = c('male', 'female')),
             script = 'rna_seq_family_heatmap.R'
)

fig_index(plot = pl_merged,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 10,
          height = 10,
          dpi = 300,
          format = 'pdf')



# === Horizontal ===

te_mean_log2_t <- t(te_mean_log2)

# Split by sex: first 3 male, last 2 female
sex_split <- c(rep("male", 3), rep("female", 2))
sex_colors <- c(female = "#efd6ac", male = "#52796f")


hm_merged_flipped <- Heatmap(te_mean_log2_t,
        width = ncol(te_mean_log2_t) * unit(3, "mm"),
        height = nrow(te_mean_log2_t) * unit(3, "mm"),
        row_names_rot = 0,
        col = col_fun,
        column_names_rot = 90,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_split = sex_split,
        left_annotation = rowAnnotation(
            block = anno_block(
                gp = gpar(fill = sex_colors),
                labels = c("female", "male"),
                labels_rot = 0,
                labels_gp = gpar(col = c("black", "white"), fontsize = 8)
            )
            ),
        heatmap_legend_param = list(title = "log2(mean(fold change))", 
                                    direction = 'horizontal',
                                    grid_width = unit(0.5, "cm"),
                                    legend_width = unit(3, 'cm'),
                                    labels_gp = gpar(fontsize = 8),
                                    title_gp = gpar(fontsize = 8)),
        border = TRUE,
        row_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8)
            
)


hm_merged_flipped <- draw(hm_merged_flipped, heatmap_legend_side = "bottom")

# Save figure with metadata to index
meta <- list(name = 'te_super_families_log2fc_heatmap_flipped',
             description = 'Heatmap of log2 mean fold change of transposable element super families in male and female',
             tags = c('expression', 'rna-seq', 'super_family'),
             parameters = list(tissues = c('brain', 'skin', 'blood'), sex = c('male', 'female')),
             script = 'rna_seq_family_heatmap.R'
)

fig_index(plot = hm_merged_flipped,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 10,
          height = 6,
          dpi = 300,
          format = 'pdf')

