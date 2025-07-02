# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}


####
# Update the directory when you defined the TEITx #
                                            #######

# orf_dirs = list.files('./results/te_island/',
#                       pattern = "orfs_",
#                       include.dirs = TRUE,
#                       full.names = TRUE)

orf_dirs = list.files('./foo/orfify/',
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

orf_pl <- ggplot(data_df, aes(tissue, V3, color = tissue)) +
    geom_boxplot() +
    labs(y = 'ORF length [bp]') +
    geom_text(data = stats_orfs, aes(label = median_length, y = median_length),
              vjust = -.5, color = 'gray70', size = 8/.pt) +
    geom_text(data = stats_orfs, 
              aes(label = paste0("n = ", format(total, big.mark = ",")), y = 85), 
              size = 8/.pt) +
    scale_color_manual(values = tissue.color) + 
    scale_y_log10() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "None",
          panel.grid = element_blank(),
          axis.title.y = element_text(size = 8, face = 'bold'))

ggsave(orf_pl,
       filename = paste0(figure_dir, 'supplemental_teitx_orf_7x7_300.pdf'),
       device = cairo_pdf,
       width = 7,
       height = 7,
       units = "cm",
       dpi = 300
)


