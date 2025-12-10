if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()

# Load canoncial TE island transcript annotations

cannonical_teIslands <- fread(te_island_transcripts)



rpl_density(cannonical_teIslands %>% 
                dplyr::rename(width = te_island.width)) +
    theme(legend.position = c(0.1,0.8))


n_func <- function(df, y = 1){
    return(data.frame(y = y, 
                      label = paste0("n = ", length(df))))
}

max_func <- function(df, y = 4.3){
    
    return(data.frame(y,
                      label = paste0("max = ", 10^max(df), " bp")))
    
}

min_func <- function(df, y = 1.35){
    
    return(data.frame(y,
                      label = paste0("min = ", 10^min(df), " bp")))
    
}

mean_func <- function(df){
    
    return(data.frame(y = mean(df),
                      label = paste0(round(10^mean(df), 2), " bp")))
    
}

cannonical_teIslands$tissue <- factor(cannonical_teIslands$tissue, levels =  c('brain', 'skin', 'blood') )

cannonical_teIslands %>% 
    filter(te_island.width > 50) %>% 
    group_by(tissue) %>% 
    summarise(max.length = max(te_island.width),
              min.length = min(te_island.width))

# Filtered for TE islands with a length bigger than 50 bp
pl_length_violin <- ggplot(cannonical_teIslands %>% filter(te_island.width > 50), aes(tissue, te_island.width)) +
    geom_violin(trim = FALSE,
                aes(fill = tissue),
                color = NA,
                #scale = "count"
    ) +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6) +
    scale_y_log10(labels = scales::comma,
                  expand = expansion(mult = c(0.05, 0.0))) +
    stat_summary(geom = "text",
                 fun.data = n_func,
                 fun.args = c(y = 1.3),
                 size = 6/.pt) +
    # stat_summary(geom = "text",
    #              fun.data = max_func,
    #              fun.args = c(y = 4.5),
    #              size = 6/.pt) +
    stat_summary(geom = "text",
                 fun.data = mean_func,
                 #fun.args = c(y = 4.5),
                 size = 6/.pt,
                 vjust = 1.5,
                 color = "black",
                 angle = 90) +
    # stat_summary(geom = "text",
    #              fun.data = min_func,
    #              fun.args = c(y = 1.4),
    #              size = 6/.pt) +
    scale_fill_manual(values = tissue.color) +
    labs(y = "Length of TE island in bp") +
    theme(legend.position = "None",
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.title.x = element_blank(),
          panel.grid = element_blank())

pl_length_violin
