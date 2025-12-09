# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()
aging_tes::load_rna_seq_env()
aging_tes::load_annotations()
aging_tes::load_analysis_env()


# Load annotations
load_te_ranges()
load_te_island_annotation()
load_te_island_instance_annotation()
load_gene_ranges()

#needs te_region_id as header for blackRcloud functions
names(te_island_instances) <- c("te_region_id", "te_id", "te_position")   

# ----- Quantification of summed TE region expression during aging in different tissues ---
# te_island_instances <- read.csv()
logmsg('Load count tables.')   
count_tables_rna <-
    sapply(names(counts_rna), simplify = FALSE,
           function(x) {
               load_salmonTE_counts(counts_rna[[x]],
                                    "region",
                                    te_island_instances)[["region"]]
           })
    
# assign the sample names to a condition (age)
conditions <- sapply(names(counts_rna), simplify = FALSE,
                         function(x) {
                             
                             getConditions(names(count_tables_rna[[x]]))
                             
                         })
    
# run DESeq
logmsg('Run DESeq2.') 
dds_te_island <- sapply(names(count_tables_rna), simplify = F,  function(x){
        
        doDEseq(count.matrix = round(count_tables_rna[[x]]),
                col_data = conditions[[x]],
                reference = "young",
                target = "all"
        )
})
    
# collect DESeq results
deseq_te_island <- sapply(names(dds_te_island), simplify = F, function(x){
        
        getDEseqResults(dds_te_island[[x]],
                        coefficient = "condition_old_vs_young",
                        FDR.filter = F)
        
    })
    


logmsg('Save the Results.') 
save(dds_te_island, file = paste0(rna_seq_deseq_dir,
                                      "te_island_rna_seq_dds.Rdata"))
    
save(deseq_te_island, file = paste0(rna_seq_deseq_dir,
                                        "te_island_rna_seq_deseq_results.Rdata"))
    

te_island_coordinates <- data.frame(te_island_5primeRanges)

deseq_te_island_merged <- 
    do.call("rbind",
            sapply(names(deseq_te_island),
                   simplify = F,
                   function(x) {
                       
                       tmp_df <- deseq_te_island[[x]] %>%
                           rownames_to_column(var = "te_island_id") %>%
                           mutate(tissue = x,
                                  diff = case_when(padj <= FDR ~ TRUE, .default = FALSE),
                                  exp_dir = case_when(log2FoldChange < 0 ~ "down", 
                                                      log2FoldChange > 0 ~ "up",
                                                      .default = "None"))
                       
                       tmp_df <- merge(tmp_df,
                                       te_island_coordinates, 
                                       by.x = "te_island_id", by.y = 'names')
                       
                       tmp_df <- tmp_df %>% 
                           mutate(te_island_coordinates = paste0(seqnames, ":", start, "-", end))
                       
                       return(tmp_df)
                       
                       
                   }))

write.table(deseq_te_island_merged, 
            file = paste0(table_dir, "02_deseq_results_te_island.csv"), 
            col.names = TRUE, 
            row.names = FALSE,
            quote = FALSE,
            sep = ',')
