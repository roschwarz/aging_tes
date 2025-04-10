
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_rna_seq_env()


# ------------------------------------------------------------------------------
# Step 1: Calculate tpm values for each tissue if results do not exist
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

# ------------------------------------------------------------------------------
# Step 2: Filter for the top 50 of DETEs and calculate the z-score
# ------------------------------------------------------------------------------

top_50 <- sapply(tissues, simplify = F, function(t){
    
    norm <-  tpm[[t]] %>%  
        filter(order %in% c('LINE', 'SINE', 'LTR', 'DNA'),
               !grepl("[?]", super_family),
               super_family != "NA") %>% 
        filter(te_id %in% (deseq.te.merged %>% 
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
    
    #the order is different in skin because of the naming
    if (t == 'skin') {
        logmsg(paste("Flip the order of the columns, to have young and old in the right order for", t))
        norm <- norm[,rev(names(norm))]
    }
    
    norm <- as.matrix(norm)
    
    norm <- t(apply(norm, 1, function(x) calc_z_score(x, capped = 2)))
    
    # removes rows full of NA --> Why such rows are there?
    norm <- norm[rowSums(is.na(norm)) != ncol(norm), ]
    
})

# ------------------------------------------------------------------------------
# Step 3: Create the heatmaps for each tissue separatley
# ------------------------------------------------------------------------------
# NOTE: The columns of young and old for skin need to be flipped afterwards in 
#       a vector graphic program. The heatmap functions doesn't keep the order
#       of the data frame that is submitted. The order of the columns are adapted
#       as well later in the vector graphic program.

brain_heatMap <- heatmap_top50(top_50$brain, 'brain')
skin_heatMap <- heatmap_top50(top_50$skin, 'skin')
blood_heatmap <- heatmap_top50(top_50$blood, 'blood')

heat_maps <- list(brain = heatmap_top50(top_50$brain, 'brain'),
                  skin = heatmap_top50(top_50$skin, 'skin'),
                  blood = heatmap_top50(top_50$blood, 'blood'))

for (hm in names(heat_maps)) {
    pdf(file = paste0(figure_dir, 'panel1_te_instance_heat_', hm, '_2.6x5.12.pdf'), width = 2.6, height = 5.12)
    draw(heat_maps[[hm]])
    dev.off()
    
}
