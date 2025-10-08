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

if (!file.exists(paste0(rna_seq_results_dir, "tpms_rna_TE_instances_SalmonTE.Rdata"))) {
    
    dds_list <- loadRdata(paste0(rna_seq_deseq_dir, "dds_TE_instances_salmonTE.Rdata"))
    
    tpm <- sapply(tissues, simplify = F, function(x){
        
        counts <- counts(dds_list[[x]])
        
        te_length <- data.frame(te_id = row.names(counts))
        
        te_length <- splitTEID(te_length, "te_id")
        
        te_length$length <- as.numeric(te_length$end) - as.numeric(te_length$start)
        
        te_length <- te_length %>% 
            dplyr::select(te_id, length)
        
        norm <- normalizeCountMatrix(counts, te_length)
        
        names(norm) <- str_replace(names(norm), "TPM." , "age_")
        
        norm <- rownames_to_column(norm, var = "te_id")
        
        norm <- norm %>% 
            gather(key = 'sample', value = "TPM", names(norm)[2:length(norm)]) %>% 
            splitTEID("te_id")
        
        
        
    })
    
    save(tpm, file = paste0(rna_seq_results_dir, "tpms_rna_TE_instances_SalmonTE.Rdata"))
    
}else{
    logmsg("Load tpm values")
    tpm <- loadRdata(paste0(rna_seq_results_dir, "tpms_rna_TE_instances_SalmonTE.Rdata"))
    
}

# ================= Female ===============================

if (!file.exists(paste0(rna_seq_results_dir, "female_tpms_rna_TE_instances_SalmonTE.Rdata"))) {
    
    dds_list <- loadRdata(paste0(rna_seq_deseq_dir_female, "dds_TE_instances_salmonTE.Rdata"))
    
    female_tpm <- sapply(names(dds_list), simplify = F, function(x){
        
        counts <- counts(dds_list[[x]])
        
        te_length <- data.frame(te_id = row.names(counts))
        
        te_length <- splitTEID(te_length, "te_id")
        
        te_length$length <- as.numeric(te_length$end) - as.numeric(te_length$start)
        
        te_length <- te_length %>% 
            dplyr::select(te_id, length)
        
        norm <- normalizeCountMatrix(counts, te_length)
        
        names(norm) <- str_replace(names(norm), "TPM." , "age_")
        
        norm <- rownames_to_column(norm, var = "te_id")
        
        norm <- norm %>% 
            gather(key = 'sample', value = "TPM", names(norm)[2:length(norm)]) %>% 
            splitTEID("te_id")
        
        
        
    })
    
    save(female_tpm, file = paste0(rna_seq_results_dir, "female_tpms_rna_TE_instances_SalmonTE.Rdata"))
    
}else{
    logmsg("Load tpm values")
    female_tpm <- loadRdata(paste0(rna_seq_results_dir, "female_tpms_rna_TE_instances_SalmonTE.Rdata"))
    
}



# ------------------------------------------------------------------------------
# Step 2: Filter for the top 50 of DETEs and calculate the z-score
# ------------------------------------------------------------------------------

top_50 <- sapply(tissues, simplify = F, function(t){
    
    norm <-  tpm[[t]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
               !grepl("[?]", super_family),
               super_family != "NA") %>% 
        filter(te_id %in% (deseq_te_merged_male %>% 
                               filter(tissue == t) %>%
                               dplyr::slice_min(order_by = padj, n = 50) %>% 
                               dplyr::slice_max(order_by = baseMean, n = 50) %>% 
                               dplyr::pull(te_id))) 
    
    norm <- norm %>% 
        separate('sample', c("age.name", "age", "sample.id", "tissue"), sep = "[_]", remove = F) %>% 
        dplyr::select(-age.name)
    
    norm <- norm %>% 
        mutate(te.annotation = paste0(order, '|', super_family, '|', family, '|', start)) %>% 
        dplyr::select(te.annotation, sample, TPM) %>% 
        spread(sample, TPM) %>%
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


female_top_50 <- sapply(names(female_tpm), simplify = F, function(t){
    
    norm <-  female_tpm[[t]] %>%  
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
        dplyr::select(te.annotation, sample, TPM) %>% 
        spread(sample, TPM) %>%
        column_to_rownames(var = 'te.annotation')
    
    #the order is different in skin because of the naming
    if (t == 'skin') {
        logmsg(paste("Reorder the columns for skin", t))
        
        order_of_samples = c(
            "age_no036_skin_18w_1_w2_R1",
            "age_no037_skin_18w_2_w2_R1",
            "age_no038_skin_18w_3_w2_R1",
            "age_no039_skin_18w_4_w2_R1",
            "age_no040_skin_18w_5_w2_R1",
            "age_no016_skin_124w_1_R1",
            "age_no032_skin_124w_2w1_R1",
            "age_no042_skin_124w_3_w2_R1",
            "age_no043_skin_124w_4_w2_R1",
            "age_no045_skin_124w_5_w3_R1"
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

brain_heatMap <- heatmap_top50(top_50$brain, 'brain')
skin_heatMap <- heatmap_top50(top_50$skin, 'skin')
blood_heatmap <- heatmap_top50(top_50$blood, 'blood', sample_size = c(4,5))
female_brain_heatMap <- heatmap_top50(female_top_50$brain, 'brain')
female_skin_heatMap <- heatmap_top50(female_top_50$skin, 'skin')

heat_maps <- list(male_brain = heatmap_top50(top_50$brain, 'brain'),
                  male_skin = heatmap_top50(top_50$skin, 'skin'),
                  male_blood = heatmap_top50(top_50$blood, 'blood', sample_size = c(4,5)),
                  female_brain = heatmap_top50(female_top_50$brain, 'brain'),
                  female_skin = heatmap_top50(female_top_50$skin, 'skin')
                  )

for (hm in names(heat_maps)) {
    
    meta <- list(name = paste0(hm, '_te_instance_heatmap_top50_'),
                 description = paste0('Heatmap of the z-score of the tpm of the top 50 differentially expressed transposable element instances in ', hm, ' tissue'),
                 tags = c('expression', 'rna-seq', 'TE', 'z-score', 'TPM'),
                 parameters = list(tissue = hm, n = 50, metric = 'z-score of TPM'),
                 script = 'rna_seq_heatmap_top_detes.R')
    
    fig_index(plot = heat_maps[[hm]],
              outdir = figure_dir,
              meta = meta,
              index_file = 'figure_index.tsv',
              width = 6.6,
              height = 13)
    
}
