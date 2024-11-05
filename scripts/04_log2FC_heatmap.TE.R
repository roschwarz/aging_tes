# --------------------------------- Notes--- -----------------------------------
# This script uses the log2foldchanges, calculated by DESeq2, and calculates a mean
# log2FoldChange of members of a super family. This shows that the analysis of TEs at a family level
# is not appropriated as we get only on value per family, which is also not that strong.
#
# Input:
#
# - deseq result table
# 
# Output:
#
# - heat map of log2(fold change) of super families

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}

if (!exists("deseq.te.merged")) {
    
    deseq.te.merged <- read.csv(paste0(table_dir, '02_deseq_results_te_instances.csv'))
    
}

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
    #mutate(mean.log2FC = case_when(mean.log2FC == 'NA' ~ 0, .default = mean.log2FC)) %>% 
    spread(key = tissue, 
           value = mean.log2FC) %>%
    # mutate(brain = replace_na(brain, 0),skin = replace_na(skin, 0),blood = replace_na(blood, 0)
    # ) %>% 
    column_to_rownames('super_family') %>% 
    as.matrix()

te.mean.log2 <- te.mean.log2[,c('brain', 'skin', 'blood')]
colorRampPalette(c('#457b9d', 'white','#e63946'))
seq(-0.1, 0.1, length.out = 40)

col_fun = colorRamp2(c(-0.1, 0, 0.1), c('#457b9d', 'white','#e63946'))
col_fun(seq(-0.1, 0.1))

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

#########
# Resis #
#########

pdf(file = './manuscripts/nature_aging/resis/Figure1/1A_heat.pdf', width = 3.78, height = 5.12)
draw(hm, heatmap_legend_side = "bottom")
dev.off()

write.csv(te.mean.log2, file = 'manuscripts/nature_aging/Figure1/resis/1A_heat.csv')


# -------------------------- means and max in tissues --------------------------

colMeans(te.mean.log2, na.rm = T)
colMaxs(abs(te.mean.log2), na.rm = T)


#### Old Version

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

pheatmap(te.mean.log2, 
         border_color = NA,
         breaks = seq(-0.1, 0.1, length.out = 40),
         color = colorRampPalette(c('#457b9d', 'white','#e63946'))(41),
         fontsize = 12,
         fontfamily = 'arial',
         cluster_rows = F,
         cluster_cols = F)

write.table(te.mean.log2, 
            file = '../../submission/Resis/figure1/figure_1a.csv', 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')



