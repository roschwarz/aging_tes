# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}


df <- data.table::fread('./results/te_island/orfs_brain/brain_te_islands_orfs_longest.bed',
                        header = FALSE) %>% 
    mutate(V3 = as.numeric(str_replace(V3, "ORF_len=", "")))

orf_dirs = list.files('./results/te_island/',
                      pattern = "orfs_",
                      include.dirs = TRUE,
                      full.names = TRUE)

data_df <- do.call('rbind', sapply(orf_dirs, simplify = FALSE, function(orf_dir){
    
    tissue = str_replace(basename(orf_dir), "orfs_", "")
    tissue = str_replace(tissue, "_downsampled", "")
    
    
    
    bed_file = list.files(orf_dir,
              pattern = "longest")
    
    bed_file = paste0(orf_dir, "/", bed_file)
    df <- data.table::fread(bed_file, header = FALSE) %>% 
    mutate(V3 = as.numeric(str_replace(V3, "ORF_len=", "")),
           tissue = tissue)
    
    return(df)
})
)


stats_orfs <- data_df %>% 
    group_by(tissue) %>% 
    dplyr::summarize(mean_length = mean(V3),
                     median_length = median(V3),
                     total = dplyr::n())

data_df$tissue <- factor(data_df$tissue, levels = tissues)

ggplot(data_df, aes(tissue, V3)) +
    geom_boxplot() +
    labs(y = 'orf length [bp]') +
    geom_text(data = stats_orfs, aes(label = median_length, y = median_length - 10),
              vjust = -8.5, angle = -90, color = 'gray70', size = 8/.pt) +
    geom_text(data = stats_orfs, aes(label = paste0("n = ", total), y = 85)) +
    scale_y_log10() +
    theme_classic()



