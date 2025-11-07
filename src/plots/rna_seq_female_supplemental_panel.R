if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_female_env()
aging_tes::load_plotting_env()

n_top = 50

# ---------------------------------- Volcanos for TEs -------------------------------------------------------- 

deseq_te_merged_female <- fread(paste0(table_dir, deseq_results_te_csv_female)) %>% 
    mutate(sex = 'female',
           sex_tissue = paste0(sex, "_", tissue))

female_brain_full <- create_te_analysis_panel(deseq_data = deseq_te_merged_female %>% filter(tissue == 'brain'), 
                                              tissue_name = 'brain', 
                                              order_colors = order.color,
                                              fdr_threshold = 0.05,
                                              base_font_size = 10,
                                              min_font_size = 8,
                                              font_family = 'sans',
                                              section = 'full',
                                              common.legend = FALSE)

female_skin_full <- create_te_analysis_panel(deseq_data = deseq_te_merged_female %>% filter(tissue == 'skin'), 
                                             tissue_name = 'skin', 
                                             order_colors = order.color,
                                             fdr_threshold = 0.05,
                                             base_font_size = 10,
                                             min_font_size = 8,
                                             font_family = 'sans',
                                             section = 'full',
                                             common.legend = FALSE)



combined_female_volcano_pl <- ggarrange(female_brain_full,
                                female_skin_full, 
                                ncol = 1,
                                nrow = 2)

# ---------------------------------- Heatmap for Top50 TEs --------------------------------------------------- 
# Load vst counts and calculate z-scores out of it


if (!exists("deseq_te_merged_female")) {
    deseq_te_merged_female <- fread(paste0(table_dir, deseq_results_te_csv_female)) %>% 
        mutate(sex = 'female',
           sex_tissue = paste0(sex, "_", tissue))
    
}


female_vst_counts_te_instances <- blackRcloud::loadRdata(paste0(rna_seq_deseq_dir, "vst_TE_instances_SalmonTE.Rdata"))

female_top_detes <- sapply(names(female_vst_counts_te_instances), simplify = F, function(t){
    
    z_scores <-  female_vst_counts_te_instances[[t]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
               !grepl("[?]", super_family),
               super_family != "NA") %>% 
        filter(te_id %in% (deseq_te_merged_female %>% 
                               filter(tissue == t) %>%
                               dplyr::slice_min(order_by = padj, n = n_top) %>% 
                               dplyr::slice_max(order_by = baseMean, n = n_top) %>% 
                               dplyr::pull(te_id))) 
    
    z_scores <- z_scores %>% 
        separate('sample', c("age.name", "id", "tissue", "age"), sep = "[_]", remove = F) %>% 
        dplyr::select(-age.name, -id) %>% 
        dplyr::mutate(age = case_when(age == '124w' ~ 'old',
                                      age == '18w' ~ 'young'))
    
    z_scores <- z_scores %>% 
        mutate(te.annotation = paste0(order, '|', super_family, '|', family, '|', start)) %>% 
        dplyr::select(te.annotation, sample, vst) %>% 
        spread(sample, vst) %>%
        column_to_rownames(var = 'te.annotation')
    
    #the order is different in skin because of the naming
    if (t == 'skin') {
        logmsg(paste("Reorder the columns for skin", t))
        
        order_of_samples = c(
            "no036_skin_18w_1_w2_R1",
            "no037_skin_18w_2_w2_R1",
            "no038_skin_18w_3_w2_R1",
            "no039_skin_18w_4_w2_R1",
            "no040_skin_18w_5_w2_R1",
            "no016_skin_124w_1_R1",
            #"no032_skin_124w_2w1_R1",
            "no042_skin_124w_3_w2_R1",
            "no043_skin_124w_4_w2_R1",
            "no045_skin_124w_5_w3_R1"
        )
        
        z_scores <- z_scores[,order_of_samples]
    }
    # 
    z_scores <- as.matrix(z_scores)
    
    z_scores <- t(apply(z_scores, 1, function(x) calc_z_score(x, capped = 2)))
    
    # removes rows full of NA --> Why such rows are there?
    z_scores <- z_scores[rowSums(is.na(z_scores)) != ncol(z_scores), ]
    
})

female_brain_heatMap <- heatmap_top50(female_top_detes$brain, 'brain', sample_size = c(5,4))
female_skin_heatMap <- heatmap_top50(female_top_detes$skin, 'skin', sample_size = c(5,4))

female_brain_heatMap_gg <- grid.grabExpr(draw(female_brain_heatMap)) %>% 
    ggplotify::as.ggplot()

female_skin_heatMap_gg <- grid.grabExpr(draw(female_skin_heatMap)) %>% 
    ggplotify::as.ggplot()

combined_female_heatmap_pl <- ggarrange(female_brain_heatMap_gg, female_skin_heatMap_gg)


pl <- ggarrange(combined_female_volcano_pl, 
          combined_female_heatmap_pl,
          ncol = 1,
          heights = c(1.6,1.3))


meta <- list(name = 'female_supplemental_panel',
             description = 'Volcano and Kimura distance plots of differentially expressed transposable elements in female brain and skin. Heatmaps of z-scores of vst counts of top50 DETEs',
             tags = c('expression', 'rna-seq', 'TE', 'kimura'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin'), sex = c('female')),
             script = 'rna_seq_female_supplemental_panel.R'
)

fig_index(plot = pl,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 18,
          height = 30,
          dpi = 300,
          format = 'pdf')
