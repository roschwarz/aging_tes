if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_rna_seq_env()

deseq_te_merged_male <- fread(paste0(table_dir, deseq_results_te_csv_male)) %>% 
    mutate(sex = 'male',
           sex_tissue = paste0(sex, "_",tissue))

deseq_te_merged_female <- fread(paste0(table_dir, deseq_results_te_csv_female)) %>%
    mutate(sex = 'female',
           sex_tissue = paste0(sex, "_", tissue))

# ------------------------------------------------------------------------------
# Step 1: Calculate TPM values for each tissue if results do not exist
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

# ================= Female ===============================

if (!file.exists(paste0(rna_seq_results_dir, "female_vst_rna_TE_instances_SalmonTE.Rdata"))) {
    
    dds_list <- loadRdata(paste0(rna_seq_deseq_dir_female, "dds_TE_instances_salmonTE.Rdata"))
    
    female_vst <- sapply(names(dds_list), simplify = F, function(x){
        
        vst_tmp <- data.frame(getVarianceStabilizedData(dds_list[[x]]))
        
        vst_tmp <- rownames_to_column(vst_tmp, var = "te_id")
        
        vst_tmp <- vst_tmp %>% 
            gather(key = 'sample', value = "vst", names(vst_tmp)[2:length(vst_tmp)]) %>% 
            splitTEID("te_id")
        
        return(vst_tmp)
        
    })
    
    save(female_vst, file = paste0(rna_seq_results_dir, "female_vst_rna_TE_instances_SalmonTE.Rdata"))
    
}else{
    logmsg("Load vst values")
    female_vst <- loadRdata(paste0(rna_seq_results_dir, "female_vst_rna_TE_instances_SalmonTE.Rdata"))
    
}



# ------------------------------------------------------------------------------
# Step 2: Filter for the top 50 of DETEs and calculate the z-score
# ------------------------------------------------------------------------------
n_top = 50

male_top_detes <- sapply(tissues, simplify = F, function(t){
    
    norm <-  vst[[t]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
               !grepl("[?]", super_family),
               super_family != "NA") %>% 
        filter(te_id %in% (deseq_te_merged_male %>% 
                               filter(tissue == t) %>%
                               dplyr::slice_min(order_by = padj, n = n_top) %>% 
                               dplyr::slice_max(order_by = baseMean, n = n_top) %>% 
                               dplyr::pull(te_id))) 
    
    norm <- norm %>% 
        separate('sample', c("age.name", "age", "sample.id", "tissue"), sep = "[_]", remove = F) %>% 
        dplyr::select(-age.name)
    
    norm <- norm %>% 
        mutate(te.annotation = paste0(order, '|', super_family, '|', family, '|', start)) %>% 
        dplyr::select(te.annotation, sample, vst) %>% 
        spread(sample, vst) %>%
        column_to_rownames(var = 'te.annotation')
    
    # #the order is different in skin because of the naming
    if (t == 'skin') {
        logmsg(paste("Flip the order of the columns, to have young and old in the right order for", t))
        norm <- norm[,rev(names(norm))]
    }

    norm <- as.matrix(norm)
    
    norm <- t(apply(norm, 1, function(x) calc_z_score(x, capped = 2)))
    
    # removes rows full of NA --> Why such rows are there?
    norm <- norm[rowSums(is.na(norm)) != ncol(norm), ]
    
})


female_top_detes <- sapply(names(female_vst), simplify = F, function(t){
    
    norm <-  female_vst[[t]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
               !grepl("[?]", super_family),
               super_family != "NA") %>% 
        filter(te_id %in% (deseq_te_merged_female %>% 
                               filter(tissue == t) %>%
                               dplyr::slice_min(order_by = padj, n = 50) %>% 
                               dplyr::slice_max(order_by = baseMean, n = 50) %>% 
                               dplyr::pull(te_id))) 
    
    norm <- norm %>% 
        separate('sample', c("age.name", "id", "tissue", "age"), sep = "[_]", remove = F) %>% 
        dplyr::select(-age.name, -id) %>% 
        dplyr::mutate(age = case_when(age == '124w' ~ 'old',
                                 age == '18w' ~ 'young'))
    
    norm <- norm %>% 
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
            "no032_skin_124w_2w1_R1",
            "no042_skin_124w_3_w2_R1",
            "no043_skin_124w_4_w2_R1",
            "no045_skin_124w_5_w3_R1"
        )
        
        norm <- norm[,order_of_samples]
    }
    # 
    norm <- as.matrix(norm)
    
    norm <- t(apply(norm, 1, function(x) calc_z_score(x, capped = 2)))
    
    # removes rows full of NA --> Why such rows are there?
    norm <- norm[rowSums(is.na(norm)) != ncol(norm), ]
    
})


# ------------------------------------------------------------------------------
# Step 3: Create the heat maps for each tissue separately
# ------------------------------------------------------------------------------

brain_heatMap <- heatmap_top50(male_top_detes$brain, 'brain')
skin_heatMap <- heatmap_top50(male_top_detes$skin, 'skin')
blood_heatmap <- heatmap_top50(male_top_detes$blood, 'blood', sample_size = c(4,5))
female_brain_heatMap <- heatmap_topdetes(female_top_detes$brain, 'brain', sample_size = c(5,4))
female_skin_heatMap <- heatmap_topdetes(female_top_detes$skin, 'skin')

heat_maps <- list(male_brain = heatmap_top50(male_top_detes$brain, 'brain'),
                  male_skin = heatmap_top50(male_top_detes$skin, 'skin'),
                  male_blood = heatmap_top50(male_top_detes$blood, 'blood', sample_size = c(4,5)),
                  female_brain = heatmap_top50(female_top_detes$brain, 'brain', sample_size = c(5,4)),
                  female_skin = heatmap_top50(female_top_detes$skin, 'skin')
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

df <- data.frame(id = c(rownames(female_top_50$brain), rownames(female_top_50$skin),
                        rownames(top_50$brain), rownames(top_50$skin), rownames(top_50$blood)),
                 sex = c(rep('female', nrow(female_top_50$brain)), rep('female', nrow(female_top_50$skin)),
                         rep('male', nrow(top_50$brain)), rep('male', nrow(top_50$skin)), rep('male', nrow(top_50$blood))),
                 tissue = c(rep('brain', nrow(female_top_50$brain)), 
                            rep('skin', nrow(female_top_50$skin)), 
                            rep('brain', nrow(top_50$brain)), 
                            rep('skin', nrow(top_50$skin)), 
                            rep('blood', nrow(top_50$blood)))
                 )


df %>% separate(col = 'id', 
                      into = c('order', 'super', 'fam', 'start'), 
                      sep = "[|]", remove = F) %>% 
    group_by(sex, tissue) %>%
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
    
# ------------------------------------------------------------------------------
# Overlap male female - Draft
# ------------------------------------------------------------------------------    
 
        
t <- 'skin'

#### Top 50 male 
    
top50_male <-  tpm[[t]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
               !grepl("[?]", super_family),
               super_family != "NA") %>% 
        filter(te_id %in% (deseq_te_merged_male %>% 
                               filter(tissue == t) %>%
                               dplyr::slice_min(order_by = padj, n = 50) %>% 
                               dplyr::slice_max(order_by = baseMean, n = 50) %>% 
                               dplyr::pull(te_id))) %>% 
    dplyr::pull(te_id) %>% unique()
    

norm <-  female_tpm[[t]] %>%  
    filter(te_id %in% top50_male,
           order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
           !grepl("[?]", super_family),
           super_family != "NA") %>%
    filter(te_id %in% (deseq_te_merged_female %>% 
                           filter(tissue == t) %>%
                           #dplyr::slice_min(order_by = padj, n = 50) %>% 
                           #dplyr::slice_max(order_by = baseMean, n = 50) %>% 
                           dplyr::pull(te_id))) 

norm <- norm %>% 
    separate('sample', c("age.name", "id", "tissue", "age"), sep = "[_]", remove = F) %>% 
    dplyr::select(-age.name, -id) %>% 
    dplyr::mutate(age = case_when(age == '124w' ~ 'old',
                                  age == '18w' ~ 'young'))

norm <- norm %>% 
    mutate(te.annotation = paste0(order, '|', super_family, '|', family, '|', start)) %>% 
    dplyr::select(te.annotation, sample, TPM) %>% 
    spread(sample, TPM) %>%
    column_to_rownames(var = 'te.annotation')

norm <- as.matrix(norm)

norm <- t(apply(norm, 1, function(x) calc_z_score(x, capped = 2)))

# removes rows full of NA --> Why such rows are there?
norm <- norm[rowSums(is.na(norm)) != ncol(norm), ]

heatmap_top50(norm, 'skin')
