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


#' ------------------------------------------------------------------------------
#' Filter a DESeq2 
#' ------------------------------------------------------------------------------
create_te_analysis_panel <- function(deseq_data, 
                                     tissue_name, 
                                     order_colors,
                                     fdr_threshold = 0.05,
                                     base_font_size = 10,
                                     min_font_size = 8,
                                     font_family = 'sans',
                                     section = 'full',
                                     ...) {
    
    # Volcano plot: all TEs
    all_volcano <- volcanoPlot(cutPvalue(deseq_data), FDR = fdr_threshold) +
        labs(title = paste('DESeq of TEs in', tissue_name)) +
        theme_rob(base_size = base_font_size, base_family = font_family) +
        theme(legend.position = 'none')
    
    # Downregulated TEs
    down_detes <- deseq_data %>% 
        filter(log2FoldChange < 0, padj <= fdr_threshold)
    
    kimura_down_pl <- ggplot(down_detes, aes(Kimura, fill = order)) +
        geom_histogram(position = 'stack', color = 'black', binwidth = 1) +
        scale_fill_manual(values = order_colors) +
        scale_x_continuous(breaks = seq(0, 50, 5)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = 'Downregulated TEs') +
        theme_bw(base_size = min_font_size) +
        theme(plot.title = element_text(hjust = 0.5, size = base_font_size),
              legend.position = 'none')
    
    # Upregulated TEs
    up_detes <- deseq_data %>% 
        filter(log2FoldChange > 0, padj <= fdr_threshold)
    
    kimura_up_pl <- ggplot(up_detes, aes(Kimura, fill = order)) +
        geom_histogram(position = 'stack', color = 'black', binwidth = 1) +
        scale_fill_manual(values = order_colors) +
        scale_x_continuous(breaks = seq(0, 50, 5)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = 'Upregulated TEs') +
        theme_bw(base_size = min_font_size) +
        theme(plot.title = element_text(hjust = 0.5, size = base_font_size),
              legend.position = 'none')
    
    # Young TEs (Kimura <= 5)
    deseq_young <- deseq_data %>% filter(Kimura <= 5)
    
    young_volcano <- volcanoPlot(cutPvalue(deseq_young), FDR = fdr_threshold) +
        labs(title = 'Young TEs (Kimura \u2264 5)') +
        theme_rob(base_size = base_font_size, base_family = font_family) +
        theme(legend.position = 'none')
    
    # Middle-age TEs (5 < Kimura <= 25)
    deseq_middle <- deseq_data %>% 
        filter(dplyr::between(Kimura, 5.001, 25))
    
    middle_volcano <- volcanoPlot(cutPvalue(deseq_middle), FDR = fdr_threshold) +
        labs(title = 'Middle old TEs (Kimura 5.001-25)') +
        theme_rob(base_size = base_font_size, base_family = font_family) +
        theme(legend.position = 'none')
    
    # Old TEs (Kimura > 25)
    deseq_old <- deseq_data %>% filter(Kimura > 25)
    
    old_volcano <- volcanoPlot(cutPvalue(deseq_old), FDR = fdr_threshold) +
        labs(title = 'Old TEs (Kimura > 25)') +
        theme_rob(base_size = base_font_size, base_family = font_family) +
        theme(legend.position = 'none')
    
    # Create combined plot with legend below
    if (section == 'top'){
        combined_plot <- ggarrange(
            kimura_down_pl, all_volcano, kimura_up_pl,
            nrow = 1, ncol = 3,
            align = 'h',
            ...)
        
        
    }else if(section == 'bottom') {
        combined_plot <- ggarrange(
            young_volcano, middle_volcano, old_volcano, 
            nrow = 1, ncol = 3,
            align = 'h',
            ...)
        
        
    }else{
        combined_plot <- ggarrange(
            kimura_down_pl, all_volcano, kimura_up_pl, 
            young_volcano, middle_volcano, old_volcano, 
            nrow = 2, ncol = 3,
            align = 'h',
            ...
    )
        
    }
    
    return(combined_plot)
}
