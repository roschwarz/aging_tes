getConditions <- function(file_names) {
    # Project specific. It uses the naming convention of the fastq file to
    # determine the tissue and age of each sample.
    
    
    samples <- data.frame(
        SampleID = file_names,
        condition = file_names)
    
    # define the conditions for each sample
    if (grepl("skin", file_names[1], ignore.case = T)) {
        samples <- samples %>%
            separate(condition,
                     c("ID", "age"),
                     sep = "[_]") %>%
            mutate(age = str_replace(age, "[0-9]*Skin", "")) %>%
            mutate(condition = as.factor(case_when(
                age == "MO" ~ "old",
                age == "MY" ~ "young"
            ))) %>%
            dplyr::select(SampleID, condition) %>% 
            column_to_rownames("SampleID")
        
        return(samples)
    }
    
    samples <- samples %>%
        separate(condition,
                 c("ID", "org", "age", "tissue", "number"),
                 sep = "[_]") %>%
        mutate(condition = as.factor(case_when(age == "o" ~ "old",
                                               age == "y" ~ "young"))) %>%
        dplyr::select(SampleID, condition) %>% 
        column_to_rownames("SampleID")
    
    
    return(samples)
    
    
}

cutPvalue <- function(df, FDR_CAP = 0.00001){
    # This function allows to set the p-value to a certain max value if the adjuste
    # p-value is below. 
    
    df <- df %>% 
        mutate(padj = case_when(padj <= FDR_CAP ~ FDR_CAP,
                                TRUE ~ padj)) 
    
    return(df)
}

#================== Coloring  strips =============================
# To-Do, make it possible to submit the colors, then I guess you can
# add it to your blackRCloud
color_strips <- function(pl, 
                         bg_cols = c('#264653', '#2A9D8F', "#E9C46A", "#af7ac5", "#F2CC8F", "#3D405B", '#81B29A'),
                         text_cols = c('#ffffff', '#ffffff', "#000000", "#000000", "#000000", "#ffffff", '#000000')){
    
    g <- ggplot_gtable(ggplot_build(pl))
    
    strip_both <- which(grepl('strip-', g$layout$name))
    
    k <- 1
    
    if (length(strip_both) != length(bg_cols)) {
        print('Sorry the number of delivired colours is different compared to the number of facetts.')
        return(g)
    }
    
    for (i in seq_along(strip_both)) {
        
        j <- which(grepl('rect', g$grobs[[strip_both[i]]]$grobs[[1]]$childrenOrder))
        l <- which(grepl('titleGrob', g$grobs[[strip_both[i]]]$grobs[[1]]$childrenOrder))
        
        g$grobs[[strip_both[i]]]$grobs[[1]]$children[[j]]$gp$fill <- bg_cols[i]
        g$grobs[[strip_both[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- text_cols[i]
        k <- k + 1
    }
    
    return(g)
} 

