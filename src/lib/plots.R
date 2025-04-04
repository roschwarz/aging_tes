######################
# Plot Manipulations #
######################


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


dropBars <- function(pl){
    
    pl + scale_y_continuous(expand = expansion(mult = c(0, .05)))
}

mplot_drop_bars <- function(bottom = 0, top = 0.1){
    
    scale_y_continuous(expand = expansion(mult = c(bottom, top)))
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


#==================== Enrichment Plots ==============================

# got from https://stackoverflow.com/questions/48000292/center-align-legend-title-and-legend-keys-in-ggplot2-for-long-legend-titles 
align_legend <- function(p, hjust = 0.5){
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]

  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")

  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]

    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1 - hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2

    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  g$grobs[[legend_index]] <- legend
  g
}


enrichmentPlot <- function(te.df.enrichment){
    # To-Do: 
    # - How to handle the labels. Currently, I have to adapt them by hand --> substitute your label vector by abs helps
    # - submit the stuff that should be shown (currently super_family) 
    # - make the tissue as header instead of axis text
    # - add enriched and depleted to the legend (see presentation of SAB)
    
    pl <- ggplot(te.df.enrichment, aes(tissue, super_family, fill = log10.padj)) +
        geom_tile(colour = 'white') +
        #geom_text(aes(label = paste0(round(ratio, 1), "(", n.target, ")" )), size =2) +
        geom_text(aes(label = label), size = 2) +
        scale_fill_gradient2(high = "#d1495b",
                             low = "#66a182",
                             mid = "white",
                             na.value = "grey90",
                             labels = abs) + 
        scale_x_discrete(position = 'top') +
        theme_bw() +
        theme_rob(12) +
        labs( fill= '-log10(adj. p-value)',
              y = 'TE family') +
        theme(axis.text.x = element_text(size = 10),
              axis.title.x = element_blank(),
              rect = element_blank(),
              panel.grid = element_blank(),
              axis.ticks = element_blank())
    
    
    pl <- ggdraw(align_legend(pl)) # center the legend
    
    return(pl)
}




# this function makes a bar plot and drops automatically the bars to the
# ground of the x-axis  
# The fill category does not work for the moment. You will get an error by
# submitting your category that you want to fill. It think it is something
# needed like '!!', however, I don't know for what I have to google.
bar <- function(ggpl, fill_cat = NULL){
    
    ggpl +
        geom_bar(aes(fill = fill_cat)) +
        scale_y_continuous(expand = expansion(mult = c(0, .05)))
}

count_raw_reads_pl <- function(df){
    
    require("wesanderson")
    
    ggplot(df, aes(age, reads/1e+6, fill=seq.approach)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 2, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) +
        facet_grid(tissue~seq.approach, scales = 'free_y') +
        scale_color_manual(values = wes_palette("Royal1")) +
        labs(y='# raw reads [Mio]') +
        theme_rob() +
        theme(panel.spacing.y = unit(0, "lines"),
              legend.position = 'bottom')
        
    
}

count_tissue_age <- function(df){
    
    require(ggplot2)
    require("wesanderson")
    require(scales) # allows to separate digit wise huge numbers by comma  
    
    #change the default theme of ggplot to my theme coded above
    
    theme_set(theme_rob(18))
    
    if(!grepl('^ENS', df$V4[1])){
        df <- df %>% filter(order %in% c('DNA', 'LINE','LTR', 'SINE'))
        
        pl <- ggplot(df, aes(order, fill=age)) +
            geom_bar(position = 'dodge') +
            labs(x='order of TE',
                y = '# affected TEs') +
            scale_fill_manual(values = wes_palette("Royal1")) +
            scale_y_continuous(expand = expansion(mult = c(0, .01)), labels = scales::comma_format()) +
            theme(
                legend.position = c(0.12,0.82),
            )
        
    }else{
        
        pl <- ggplot(df, aes(tissue, fill=age)) +
        geom_bar(position = 'dodge') +
        labs(y = '# affected genes') +
        scale_fill_manual(values = wes_palette("Royal1")) +
        scale_y_continuous(expand = expansion(mult = c(0, .01)), labels = scales::comma_format()) +
        theme(
            
            legend.position = c(0.12,0.82),
        )
    }
    
    if(length(unique(df$tissue)) > 1){
        return(pl + facet_grid(.~tissue, scales = 'free_x'))
    }else{
        return(pl)
    }
    
}


correlationPanel <- function(df, pannel.group = 'sig'){
    
# This Panel shows the correlation of the log2FC between different sequencing 
# strategies. The correlations separated into elements that are significantly 
# differential expressed in both, only one of both or in none of the sequencing 
# strategies.
    categories <- unique(df[[pannel.group]])
    log2FC <- names(df)[grepl('log2FoldChange', names(df))]
    
    maxScale <- max(c(abs(df[[log2FC[1]]]), c(abs(df[[log2FC[2]]])))) # determine the max log2FC to make the plot symmetric
    
    pl <- ggplot(df, aes_string(log2FC[1], log2FC[2], color = pannel.group)) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_point(alpha = 0.5) +
        xlim(c(-maxScale, maxScale)) +
        ylim(c(-maxScale, maxScale)) +
        geom_smooth(method = 'lm') +
        facet_wrap(c(pannel.group), ncol = 2, nrow = 2) +
        labs(color = "sig. expression change in") +
        theme(panel.grid.major = element_blank(),
              legend.position = 'bottom'
        ) +
        scale_color_tableau()
    
    pl <- pl+ggpubr::stat_cor(aes(color = pannel.group))
    
    return(pl)
    
}


###------------------------ General plots generated by my own ------------------------------------------------

# rpl stands for robert plots
rpl_count <- function(df){
    
    ggplot(df, aes(tissues, n, fill = tissues)) +
        geom_col() +
        scale_fill_manual(values = tissue.color) +
        mplot_drop_bars() + 
        geom_text(
            x = tissues,
            y = 1000,
            aes(label = paste0("n=", n)),
            vjust = -2.5,
            hjust = 0,
            angle = 90,
            color = "white"
        ) +
        theme(legend.position = "None")
    
}


rpl_density <- function(df){
    
    
    if (!"width" %in% names(df)) {
        return("Error, column width not contained in your data frame.")
    }
    
    length_mean <- df %>%
        group_by(tissue) %>%
        summarize(median = median(width))
    
    ggplot(df, aes(width, after_stat(scaled) , color = tissue)) +
        geom_density(linewidth = 1.2) +
        scale_x_continuous(trans = 'log10') +
        mplot_drop_bars() +
        geom_vline(
            data = length_mean,
            aes(xintercept = median, color = tissue),
            linewidth = 1.2,
            linetype = 2,
            alpha = 0.5
        ) +
        scale_color_manual(values = tissue.color) +
        labs(x = "length in bp") +
        theme(legend.position = "bottom")
}


rpl_composition <- function(df, tissue = NA){
    
    pl <- ggplot(df, aes(x = "", y = per, fill = order)) +
        geom_col(color = "black") +
        geom_label(
            aes(label = labels),
            color = c("white", "white", "white", "black"),
            position = position_stack(vjust = 0.5),
            show.legend = FALSE
        ) +
        coord_polar("y", start = 0) +
        guides(fill = guide_legend(title = "Order of TEs")) +
        scale_fill_manual(values = order.color) +
        theme_void() +
        theme(legend.position =  "bottom")
    
    if (!is.na(tissue)) {
        pl <- pl + ggtitle(tissue) + theme(plot.title = element_text(hjust = 0.5))
    }
    
    return(pl)
    
}

rpl_super_composition <- function(df){
    
    
    df %>% filter(per >= 0.5) %>%  
        ggplot(aes(super_family, per, fill = order)) +
        geom_col() +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        coord_flip() +
        scale_fill_manual(values = order.color) +
        theme(legend.position = 'None',
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
}


rpl_position_density <- function(df){
    ggplot(df, aes(te_rel_start, color = order)) +
        geom_density() +
        scale_color_manual(values = color_palette) +
        labs(x = "relative position in TE Island") +
        theme_bw()
}


make_composition_panel <- function(expr_te_regions, teRegionInstances){
    
    expr_teRegion_composition <-
        sapply(tissues, simplify = FALSE, function(x) {
            islands_of_interest <- expr_te_regions %>%
                filter(tissue == x) %>%
                dplyr::pull(names)
            
            df <- teRegionInstances %>%
                filter(te_region_id %in% islands_of_interest)
            
            super_per <- prepare_super_composition(df)
            
            order_counts <- prepare_Composition(df)
            
            
            pie <- plot_composition(order_counts, tissue = x)
            bar <- plot_super_composition(super_per)
            
            grid.arrange(
                pie +
                    theme(
                        legend.position = 'None',
                        plot.margin = unit(c(0, -10, -3, -10), "mm")
                    ),
                heights = c(0.6, 3),
                widths = c(3),
                bar
            )
            
        })
    
    pl <- ggpubr::ggarrange(plotlist = expr_teRegion_composition,
                            nrow = 1)
    
    
    pl <- annotate_figure(pl, bottom = "TE island coverage (%)", left = "super family")
    
    
    pl <- grid.arrange(pl, order_legend, ncol = 1, heights = c(9, 0.3))
    
    return(pl)
}

## Function that extract the legend of a pl

legendBuilder <- function(df) {
    
    dummy_order_counts <- prepare_Composition(df)
    
    dummy_pie <- rpl_composition(dummy_order_counts) +
        theme(legend.title = element_text(size = 15))
    
    return(get_legend(dummy_pie))
    
}



countTEsPl <- function(list.of.results, fdr = 0.05){
    # DESeq2 was run for different tissues and the results
    # are stored in a list, which is the input of that 
    # function.
    
    df <- do.call('rbind', lapply(names(list.of.results), function(x) {
        
        sig <- getDiffTes(list.of.results[[x]], fdr)
        sig$condition <- x
        
        return(sig)
    }))
    
    df <- splitTEID(df, 'TE.ID')
    df$direction <- ifelse(df$log2FoldChange > 0, 'up', 'down') 
    
    count.instances <- df %>% 
        group_by(condition, direction) %>% 
        dplyr::count() 
    
    count.subfamilies <- df %>% 
        group_by(condition, direction, family) %>% 
        dplyr::count()

    count.instances$n <- ifelse(count.instances$direction == 'down', -(count.instances$n), count.instances$n)
    count.subfamilies$n <- ifelse(count.subfamilies$direction == 'down', -(count.subfamilies$n), count.subfamilies$n)
    
    pl.instances <- ggplot(count.instances, aes(condition, n, fill = direction)) +
            geom_bar(stat="identity", position = "identity") +
            labs(title = 'Instance level', y = 'No. of TE instances', x = 'tissue', fill = 'Differential expression direction: ') +
            scale_fill_tableau() +
            theme_rob(14) +
            theme(axis.title.x = element_blank())
    
    pl.subfamilies <-  ggplot(count.subfamilies, aes(condition, n, fill = direction)) +
        geom_bar(stat="identity", position = "identity") +
        labs(title = 'Subfamily level', y = 'No. of TE subfamilies', x = 'tissue', fill = 'Differential expression direction: ') +
        scale_fill_tableau() +
        theme_rob(14) +
        theme(axis.title.x = element_blank())
    
    pannel <- ggarrange(pl.instances, pl.subfamilies,
              ncol = 2,
              common.legend = T, legend = 'bottom') 

    return(pannel)
    
}


geneTECountPl <- function(df,  tissue = NULL){
    
    df <- filter(df, !duplicated(ensembl_gene_id))
    
    df$ensembl_gene_id <- factor(df$ensembl_gene_id, levels = df[order(-df$number.of.TE),"ensembl_gene_id"])
    
    ggplot(df, aes(ensembl_gene_id, number.of.TE)) +
        geom_col() +
        theme_rob() +
        scale_y_continuous(expand = expansion(mult = c(0, .02))) +
        labs(y = "# of TEs in genomic env",
             x = "differentially expressed genes",
             title = tissue)+
        theme(axis.text.x = element_blank(),
              panel.grid = element_blank(),
              axis.ticks.x = element_blank(), 
              panel.spacing.y = unit(0, "lines"))
}


geneTElogPl <- function(df, sort.by = 'mean', tissue = NULL){
    
    df <- filter(df, !duplicated(ensembl_gene_id))
    y = paste0('te.', sort.by, '.log2FC')
    order.of.genes <- arrange(df, desc(!!as.symbol(y))) %>% pull(ensembl_gene_id)
    
    df$ensembl_gene_id <- factor(df$ensembl_gene_id, levels = order.of.genes)
    
    pl <- ggplot(df, aes(ensembl_gene_id, !!as.symbol(y), fill = gene.exp.direction)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        labs(y = paste(sort.by, "log2FC of associated TEs"),
             x = 'Ensembl Gene Id',
             fill = 'expr. direction gene',
             title = tissue) +
        theme_rob() +
        theme(axis.text.x = element_text(angle = 90),
              panel.grid = element_blank(),
              legend.position = c(0.12,0.12)) +
        scale_fill_manual(values = c(down = '#66a182',up = '#d1495b'))
    
    if(nrow(df) > 100){
        pl <- pl + 
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
    }
    
    return(pl)


}


localCountPl <- function(df){
    ggplot(df, aes(localization)) +
        geom_bar(stat = "count") +
        scale_y_log10() +
        stat_count(geom = "text", size=3.5, aes(label = ..count..), position = position_stack(vjust = 1.2)) +
        facet_grid(tissue~type)
}


arrangeEnrichmentPlots <- function(list.of.plots){
    
    if(length(list.of.plots) == 0 || is.null(list.of.plots)){
        print('Nothing found!')
        return(NULL)
    }
    
    return(ggarrange(plotlist = list.of.plots, nrow = 1))
}

barplot.mod <- function(go.results, header = NULL){
    
    
   pl <- goBarPlot(go.results, header)
    
   return(pl)
   # if(is.null(go.results)){
   #     return(NULL)
   # }
   #  
   #  require(wesanderson)
   #  pal <- wes_palette("Zissou1", 100, type = "continuous")
   #  
   #  barplot(go.results, showCategory=20) + 
   #      labs(title = header)  + 
   #      theme_rob(12) +
   #      scale_fill_gradientn(colours = pal) + 
   #      theme(legend.position = 'bottom')
}

goBarPlot <- function(go.results, top = 10, header = NULL){
    
   if(is.null(go.results)){
       return(NULL)
   }
    
    if(class(go.results) != 'data.frame'){
        go.results <- as.data.frame(go.results)
    }
    
    go.results$Description <- factor(go.results$Description, 
                                     level=go.results[order(go.results$p.adjust, decreasing = T),'Description'])
    
    require(wesanderson)
    pal <- rev(wes_palette("Zissou1", 100, type = "continuous"))
    
    if(!is.null(top)){
        go.results <- go.results[1:10,] 
    }
    
    ggplot(go.results, aes(Count, Description, fill=p.adjust)) +
        geom_col() +
        theme_rob(15, base_family = 'arial') +
        labs(title = header) +
        scale_fill_gradientn(colours = pal)  
}

goBarFacetUpdate <- function(results, top = 10){
    
    goBarFacet(results, top, adjP = T) +
        theme(axis.title.y = element_blank(),
              legend.position = 'none')
}

goBarPPlot <- function(go.results, top=10, header = NULL){
    
    if(is.null(go.results)){
        return(NULL)
    }
    
    if(class(go.results) != 'data.frame'){
        go.results <- as.data.frame(go.results)
    }
    
    go.results$Description <- factor(go.results$Description, 
                                     level=go.results[order(go.results$p.adjust, decreasing = T),'Description'])
    
    require(wesanderson)
    pal <- rev(wes_palette("Zissou1", 100, type = "continuous"))
    
    
    if(!is.null(top)){
        go.results <- go.results[1:10,] 
    }
    
    go.results <- mutate(go.results,
                         ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                              ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                              TRUE ~ 'Molecular\nFunction'))
    
    x.max = max(-log10(go.results$p.adjust))+1
        
    ggplot(go.results, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
        geom_col() +
        theme_rob(15) +
        geom_text(aes(label = Count, hjust = 1.5), color = 'white') +
        scale_fill_manual(values = wes_palette("FantasticFox1")) +
        scale_x_continuous(expand = expansion(mult = c(0, .1)),
                           limits = c(0, round(x.max))) +
        labs(title = header) +
            theme(axis.ticks = element_blank(),
              panel.grid.major.y = element_blank(), 
              strip.text.y = element_text(angle = 360),
              legend.position = 'None',
              strip.background = element_rect(color = 'white'))
    
}

goBarFacet <- function(results, top=10, adjP = FALSE){
    
    if(nrow(results) == 0){
        return(NULL)
    }
    results <- results %>% 
        group_by(ONTOLOGY) %>% 
        dplyr::slice_min(order_by = p.adjust, n = top) %>% 
        ungroup() %>% 
        as.data.frame()
    
    if(adjP){
        
        
        pl <- goBarPPlot(results, top = top) +
            facet_grid(ONTOLOGY~., scales = 'free_y')+
            theme(panel.spacing.y = unit(0.25, "lines"),
                  panel.border = element_blank(),
                  axis.line = element_line())
        
        
        return(pl)
        
    }
    
    pl <- goBarPlot(results) +
        facet_grid(ONTOLOGY~.,  scales = 'free_y')
    
    return(pl)
    
}

dotplot.mod <- function(go.results, header = NULL){
    
    if(is.null(go.results)){
        return(NULL)
    }
    
    require(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    
    dotplot(go.results) + 
        labs(title = header)  + 
        theme_rob(12) +
        scale_fill_gradientn(colours = pal) + 
        theme(legend.position = 'bottom')
}


#======================= Volcano Plots =================================

volPlot <- function(df, FDR=0.05, log2FC = 1, top_gen = 10, type = 'gene'){
    
    # This function creates a standard volcano plot and requires a data frame
    # with an adjusted p-value (padj), a log2FoldChange and a column with a
    # the external gene name. The plot contains dashed lines by the cutoffs
    # of the FDR (default 0.05) and log2FC (default 1). The top ten of the 
    # significant genes (sorted by log2FoldChange) are named by default. Can
    # be set by the top_gen attribute.
    # The type attribute gives the possibility to label the TEs by the TE ids.
    # When the type attribute is set to gene that the label of the TEs are the
    # associated genes.
    
    if(class(df) != "data.frame"){
        df <- as.data.frame(df)
    }
    
    if(type %in% c('TE', 'te', 'Te')){
        # rename the TE.ID to external_gene_name to avoid a second plot call for TE.ID
        df <- df[c('log2FoldChange', 'padj', 'TE.ID')]
        names(df) <- c('log2FoldChange', 'padj', 'external_gene_name')
    }
    
    max_x = max(abs(df[['log2FoldChange']])) 
    
    if(is.null(log2FC)){
        df <- mutate(df, threshold_OE = padj < FDR)
        
    }else{
        
        print('Hi') 
        df <- mutate(df, threshold_OE = padj < FDR & abs(log2FoldChange) >= log2FC)
        
    }
    
    
    # create labels
    df.sig <- df[df$threshold_OE,]
    
    top_genes <- df.sig[order(abs(df.sig$log2FoldChange), decreasing = T), 'external_gene_name'] %>% 
        head(n=top_gen)
   
     
    ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(data = df[!df$threshold_OE,], color = 'black', alpha = 0.5) +
        geom_point(data = df[df$threshold_OE,], color = 'red') +
        xlab("log2 fold change") + 
        xlim(-(max_x), max_x) +
        geom_hline(yintercept = 0) +
        geom_hline(yintercept = -log10(FDR), linetype = 2) +
        geom_vline(xintercept = -(log2FC), linetype = 2) +
        geom_vline(xintercept = log2FC, linetype = 2) +
        geom_text_repel(data = df.sig %>% filter(external_gene_name %in% top_genes),
                        aes(label = external_gene_name)) +
        ylab("-log10 adjusted p-value") +
        theme_rob() +
        theme(panel.grid.major = element_blank())
    
}


volPlot_associated <- function(df, FDR=0.05, log2FC = NULL, top_gen = 10, type = 'gene', facet = FALSE){
    
    # This function creates a standard volcano plot and requires a data frame
    # with an adjusted p-value (padj), a log2FoldChange and a column with a
    # the external gene name. The plot contains dashed lines by the thresholds
    # of the FDR (default 0.05) and log2FC (default 1). The top ten of the 
    # significant genes (sorted by log2FoldChange) are named by default. Can
    # be set by the top_gen attribute.
    # The type attribute gives the possibility to label the TEs by the TE ids.
    # When the type attribute is set to gene that label of the TEs are the
    # associated genes.
    
    if(class(df) != "data.frame"){
        df <- as.data.frame(df)
    }
    
    if(type %in% c('TE', 'te', 'Te')){
       
        df <- df[c('log2FoldChange', 'padj', 'TE.ID', 'localization', 'direction')]
        names(df) <- c('log2FoldChange', 'padj', 'external_gene_name', 'localization', 'direction')
    }
    
    max_x = max(abs(df[['log2FoldChange']])) 
    
    if(is.null(log2FC)){
        df <- mutate(df, threshold_OE = padj < FDR)
        
    }else{
        
        df <- mutate(df, threshold_OE = padj < FDR & abs(log2FoldChange) >= log2FC)
        
    }
    
    # create labels
    df.sig <- df[df$threshold_OE,]
    
    if(facet){
        top_genes <- do.call('c', sapply(c('intergenic', 'intragenic'), 
                                         
                                         simplify = F, function(x) {
                                             
                                             genes.ordered <- df.sig %>% 
                                                 filter(localization == x) %>%
                                                 arrange(desc(abs(log2FoldChange))) %>% 
                                                 pull(external_gene_name) %>% 
                                                 unique()
                                             
                                             return(discard(genes.ordered[1:top_gen], is.na))
                                             
                                         }))
    }else{
        
        top_genes <- df.sig[order(abs(df.sig$log2FoldChange), decreasing = T), 'external_gene_name'] %>% unique() %>% 
            head(n=top_gen)
    }
    
    pl <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
        # non-significant data points
        geom_point(data = df[!df$threshold_OE,], color = 'black', alpha = 0.5) +
        #significant data points
        geom_point(data = df[df$threshold_OE,], aes(color = direction, 
                                                    shape = localization), 
                   size = 3, 
                   alpha=0.5) +
        scale_color_manual(values=c(up = '#d1495b', down = '#66a182')) +
        scale_y_continuous(expand = expansion(mult = c(0, .05))) +
        geom_hline(yintercept = 0) +
        geom_hline(yintercept = -log10(FDR), linetype = 2, color = 'grey') +
        # put the FDR Cutoff the the dashed line the coordinates are 
        # proportions of the max_x and FDR values
        geom_text(aes(x=-max_x + (5*max_x/100), 
                      label=paste("adjusted p-value of", FDR), 
                      y=-log10(FDR-(20*FDR/100))), 
                  colour="grey") +
        geom_text_repel(data = df.sig %>% filter(external_gene_name %in% top_genes),
                        aes(label = external_gene_name),
                        size = 4, 
                        color = "black", 
                        max.overlaps = 20,
                        direction = "both", 
                        box.padding = unit(0.45, "lines"),
                        point.padding = unit(0.3, "lines")) +
        xlab("log2 fold change of DEG associated TEs") + 
        xlim(-(max_x), max_x) +
        ylab("-log10 adjusted p-value of DEG associated TEs") +
        labs(color = "Direction of DEG expression",
             shape = "Localization of DEG associated TEs") +
        theme_rob() +
        theme(panel.grid.major = element_blank(),
              legend.position = "bottom")
    
    if(facet){
        
        pl <- pl + facet_grid(.~localization)
        
    }
    
    if(!is.null(log2FC)){
        
        pl <- pl + geom_vline(xintercept = -(log2FC), linetype = 2, color = 'grey') +
        geom_vline(xintercept = log2FC, linetype = 2, color = 'grey') 
        
    }
    
    return(pl)
    
}


volcanoPlot <- function(df,
                        FDR = 0.05,
                        facet = NULL,
                        rastering = T) {
    # FDR - sets the threshold when a measurement is considered as significant
    # facet - takes a string and makes the respective facet for that group
    # rastering - if TRUE all points that are not significant are combined and cannot be select in a graphic tool anymore
    require(tidyverse)
    require(ggrastr)
    
    df <- df %>% filter(!is.na(padj))
    x_max <- max(abs(df$log2FoldChange)) + 0.5
    
    #colors <- setNames(c("#d1495b", "gray", "#66a182"), c("Up", "Not Sig", "Down"))
    colors <-
        setNames(c("#e63946", "gray", "#457b9d"), c("Up", "Not Sig", "Down"))
    
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
            scale_color_manual(values = colors) +
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
            scale_color_manual(values = colors) +
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



#=======================   Venn-Diagrams  =================================

vennDia <- function(vector.list = NULL, 
                    categories = NULL, 
                    file_name = NULL,
                    header = NULL, 
                    log = TRUE){
    require(VennDiagram)
    
    diagram <- venn.diagram(
        x = vector.list,
        category.names = categories,
        
        # Output
        filename = file_name,
        imagetype = "svg",
        output=TRUE,
        width = 200,
        height = 500,
        resolution = 300,
        disable.logging = TRUE,
        
        
        # Circles
        lwd = 2,  
        fill=c('#264653', '#2A9D8F',  '#E9C46A'),
        alpha = c(0.5, 0.5, 0.5),
        
        
        # Number
        cex = .6,
        fontface = "bold",
        fontfamily = "arial",
        # 
        # # Set names
        cat.cex = 0.6,
        cat.default.pos = "outer",
        cat.fontface = "bold",
        cat.fontfamily = "arial",
        # cat.pos = c(-27, 27),
        # cat.dist = c(0.055, 0.055),
        main = header,
        scaled=T
    )
    
    if(!log){
        
        # remove the log files that are generated by venn.diagram
        file.remove(list.files('.', pattern = 'VennDiagram*'))
    }
    
    return(diagram)
}
