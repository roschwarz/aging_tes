if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_rna_seq_env()

deseq_te_merged_male <- fread(paste0(table_dir, deseq_results_te_csv_male)) %>% 
    mutate(sex = 'male',
           sex_tissue = paste0(sex, "_",tissue))

# ------------------------------------------------------------------------------
# Step 1: Take vst count from the dds objects got from DESeq2
# ------------------------------------------------------------------------------

if (!file.exists(paste0(rna_seq_results_dir, "vst_rna_TE_instances_SalmonTE.Rdata"))) {
    
    dds_list <- loadRdata(paste0(rna_seq_deseq_dir, "dds_TE_instances_salmonTE.Rdata"))
    
    vst <- sapply(tissues, simplify = F, function(x){
        
        vst_tmp <- data.frame(getVarianceStabilizedData(dds_list[[x]]))
        
        vst_tmp <- rownames_to_column(vst_tmp, var = "te_id")
        
        vst_tmp <- vst_tmp %>% 
            gather(key = 'sample', value = "vst", names(vst_tmp)[2:length(vst_tmp)]) %>% 
            splitTEID("te_id")
        
        return(vst_tmp)
        
    })
    
    save(vst, file = paste0(rna_seq_results_dir, "vst_rna_TE_instances_SalmonTE.Rdata"))
    
}else{
    
    logmsg("Load vst values")
    vst <- loadRdata(paste0(rna_seq_results_dir, "vst_rna_TE_instances_SalmonTE.Rdata"))
    
}

# ------------------------------------------------------------------------------
# Step 2: Filter for the top 50 of DETEs and calculate the z-score
# ------------------------------------------------------------------------------

n_top = 50

male_top_detes <- sapply(tissues, simplify = F, function(t){
    
    z_score <-  vst[[t]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
               !grepl("[?]", super_family),
               super_family != "NA") %>% 
        filter(te_id %in% (deseq_te_merged_male %>% 
                               filter(tissue == t) %>%
                               dplyr::slice_min(order_by = padj, n = n_top) %>% 
                               dplyr::slice_max(order_by = baseMean, n = n_top) %>% 
                               dplyr::pull(te_id))) 
    
    z_score <- z_score %>% 
        separate('sample', c("age.name", "age", "sample.id", "tissue"), sep = "[_]", remove = F) %>% 
        dplyr::select(-age.name)
    
    z_score <- z_score %>% 
        mutate(te.annotation = paste0(order, '|', super_family, '|', family, '|', start)) %>% 
        dplyr::select(te.annotation, sample, vst) %>% 
        spread(sample, vst) %>%
        column_to_rownames(var = 'te.annotation')
    
    # #the order is different in skin because of the naming
    if (t == 'skin') {
        logmsg(paste("Flip the order of the columns, to have young and old in the right order for", t))
        z_score <- z_score[,rev(names(z_score))]
    }

    z_score <- as.matrix(z_score)
    
    z_score <- t(apply(z_score, 1, function(x) calc_z_score(x, capped = 2)))
    
    # removes rows full of NA --> Why such rows are there?
    z_score <- z_score[rowSums(is.na(z_score)) != ncol(z_score), ]
    
})


# ------------------------------------------------------------------------------
# Step 3: Create the heat maps for each tissue separately
# ------------------------------------------------------------------------------

brain_heatMap <- heatmap_top50(male_top_detes$brain, 'brain')
skin_heatMap <- heatmap_top50(male_top_detes$skin, 'skin')
blood_heatmap <- heatmap_top50(male_top_detes$blood, 'blood', sample_size = c(4,5))

heat_maps <- list(male_brain = heatmap_top50(male_top_detes$brain, 'brain'),
                  male_skin = heatmap_top50(male_top_detes$skin, 'skin'),
                  male_blood = heatmap_top50(male_top_detes$blood, 'blood', sample_size = c(4,5)),
                  )

for (hm in names(heat_maps)) {
    
    meta <- list(name = paste0(hm, '_te_instance_heatmap_top50_vst'),
                 description = paste0('Heatmap of the z-score of the vst of the top 50 differentially expressed transposable element instances in ', hm, ' tissue'),
                 tags = c('expression', 'rna-seq', 'TE', 'z-score', 'vst'),
                 parameters = list(tissue = hm, n = 50, metric = 'z-score of vst'),
                 script = 'rna_seq_heatmap_top_detes.R')
    
    fig_index(plot = heat_maps[[hm]],
              outdir = figure_dir,
              meta = meta,
              index_file = 'figure_index.tsv',
              width = 6.6,
              height = 13)
    
}

# ------------------------------------------------------------------------------
# Step 4: Collect some numbers
# ------------------------------------------------------------------------------

df <- data.frame(id = c(rownames(male_top_detes$brain), rownames(male_top_detes$skin), rownames(male_top_detes$blood)),
                 sex = c(rep('male', nrow(male_top_detes$brain)), rep('male', nrow(male_top_detes$skin)), rep('male', nrow(male_top_detes$blood))),
                 tissue = c(rep('brain', nrow(male_top_detes$brain)), 
                            rep('skin', nrow(male_top_detes$skin)), 
                            rep('blood', nrow(male_top_detes$blood)))
                 )


df %>% separate(col = 'id', 
                      into = c('order', 'super', 'fam', 'start'), 
                      sep = "[|]", remove = F) %>% 
    group_by(tissue) %>%
    summarise(n = n(),
              LINE = sum(order == 'LINE'),
              SINE = sum(order == 'SINE'),
              LTR = sum(order == 'LTR'),
              DNA = sum(order == 'DNA')) %>%
    ungroup()

df %>% separate(col = 'id', 
                      into = c('order', 'super', 'fam', 'start'), 
                      sep = "[|]", remove = F) %>% 
    group_by(sex, tissue, super) %>%
    summarise(n = n()) %>% view()
    print(n = 60)