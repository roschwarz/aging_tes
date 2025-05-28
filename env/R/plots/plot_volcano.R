#' ------------------------------------------------------------------------------
#' Volcano plots, DESeq2 results
#' ------------------------------------------------------------------------------
#' @param results A DESeq2 result object
#' @return A ggplot2 object
#' @export
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
            scale_color_manual(values = direction.color) +
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
                color = direction.color[['up']],
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
                color = direction.color[['down']],
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
                     color = direction.color[['up']], 
                     label = paste0("n=", n_up_regulated)) +
            annotate(geom = 'label',
                     x = -x_max, y = y_count_label,
                     hjust = 0,
                     color = direction.color[['down']],
                     label = paste0("n=", n_down_regulated))
        
    }
    
    return(pl)
}
