# from https://benjaminlouis-stat.fr/en/blog/2020-05-21-astuces-ggplot-rmarkdown/
theme_rob <- function(base_size = 14, base_family = 'sans'){
    
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            
            plot.title = element_text(size = rel(1), margin = margin(0,0,5,0), hjust = 0.5),
            
            panel.grid.minor = element_blank(),
            
            # axes
            axis.title = element_text(size = rel(0.85), face = "bold"),
            axis.text = element_text(size = rel(0.70)),
            
            # Legend
            legend.title = element_text(size = rel(0.85), face = "bold"),
            legend.text = element_text(size = rel(0.70)),
            legend.background = element_rect(colour = 'lightgrey'),
            
            # Facett stuff
            panel.spacing.x = unit(0, "lines"),
            strip.background = element_rect(fill = 'white'),
            strip.text = element_text(size = rel(0.85), face = "bold")
        )
    
    
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

volcanoPlot <- function(df,
                        FDR = 0.05,
                        facet = NULL,
                        rastering = TRUE) {
    # FDR - sets the threshold when a measurement is considered as significant
    # facet - takes a string and makes the respective facet for that group
    # rastering - if TRUE all points that are not significant are combined and cannot be select in a graphic tool anymore
    require(tidyverse)
    require(ggrastr)
    
    df <- df %>% filter(!is.na(padj))
    x_max <- max(abs(df$log2FoldChange)) + 0.5
    
    # #colors <- setNames(c("#d1495b", "gray", "#66a182"), c("Up", "Not Sig", "Down"))
    # colors <-
    #     setNames(c("#e63946", "gray", "#457b9d"), c("Up", "Not Sig", "Down"))
    
    df$Expression <-
        ifelse(
            df$padj <= FDR &  df$log2FoldChange > 0,
            'Up',
            ifelse(df$padj <= FDR &
                       df$log2FoldChange < 0, 'Down', "Not Sig")
        )
    
    
    
    if (rastering) {
        pl <- ggplot(df, aes(log2FoldChange,-log10(padj))) +
            geom_point_rast(
                data = df %>% filter(Expression == 'Not Sig'),
                aes(color = Expression),
                size = 0.7
            ) +
            geom_point(
                data = df %>% filter(Expression %in% c('Up', 'Down')),
                aes(color = Expression),
                size = 0.7
            ) +
            geom_hline(yintercept = -log10(FDR),
                       linetype = 2) +
            geom_vline(xintercept = 0, linetype = 2) +
            xlim(c(-x_max, x_max)) +
            labs(y = bquote(-log[10](FDR)),
                 x = expression(log[2] * "(fold change)")) +
            scale_color_manual(values = volcano_colors) +
            theme_bw() +
            theme(
                legend.title = element_blank(),
                panel.grid = element_blank(),
                panel.grid.major = element_blank()
            )
        
    } else{
        pl <- ggplot(df, aes(log2FoldChange,-log10(padj))) +
            geom_point(aes(color = Expression), size = 0.7) +
            geom_hline(yintercept = -log10(FDR),
                       linetype = 2) +
            xlim(c(-x_max, x_max)) +
            labs(y = "-log10(FDR)",
                 x = "logFC") +
            scale_color_manual(values = volcano_colors) +
            theme_bw() +
            theme(
                legend.title = element_blank(),
                panel.grid = element_blank(),
                panel.grid.major = element_blank()
            )
    }
    
    if (!is.null(facet)) {
        
        y_count_label = (-log10(min(df$padj))) - 0.2
        n_up_regulated = df %>% filter(Expression == 'Up') %>% dplyr::count(!!sym(facet))
        n_down_regulated = df %>% filter(Expression == 'Down') %>% dplyr::count(!!sym(facet))
        
        pl <- pl +
            facet_grid(cols = vars(!!sym(facet))) +
            geom_text(
                data = n_up_regulated,
                color = "#e63946",
                mapping = aes(
                    x = x_max,
                    y = y_count_label,
                    label = paste0("n=", n)
                ),
                hjust = 1,
                size = 8/.pt
            ) +
            geom_text(
                data = n_down_regulated,
                color = "#457b9d",
                mapping = aes(
                    x = -x_max,
                    y = y_count_label,
                    label = paste0("n=", n)
                ),
                hjust = 0,
                size = 8/.pt
            )
    }else{
        
        y_count_label = (-log10(min(df$padj))) - 0.2
        n_up_regulated = df %>% filter(Expression == 'Up') %>% nrow()
        n_down_regulated = df %>% filter(Expression == 'Down') %>% nrow()
        
        
        pl <- pl +
            annotate(geom = 'label',
                     x = x_max, y = y_count_label,
                     hjust = 1,
                     color = "#e63946", 
                     label = paste0("n=", n_up_regulated)) +
            annotate(geom = 'label',
                     x = -x_max, y = y_count_label,
                     hjust = 0,
                     color = "#457b9d",
                     label = paste0("n=", n_down_regulated))
        
    }
    
    return(pl)
}