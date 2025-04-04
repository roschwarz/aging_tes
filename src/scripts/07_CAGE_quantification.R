# --------------------------------- Notes--- -----------------------------------
# 
# This script is used for the CAGE-Seq quantification
#
# Output data (../data/cage; ../tables):
#
# [tissue].peakAnnotation - bedfile of all called peaks
#
# dds.cage.Rdata - dds object for called cage peaks (peaks with less than 10 
#    reads across all samples are removed) 
#
# deseq.cage.Rdata - deseq results for all called peaks 
#
# 07_deseq_cage_all.csv - deseq results of all called peaks
#
# 07_deseq_cage_gene.csv - deseq results of all called peaks that intersect with 
#   a gene. The longest version of the gene was used for that intersection and
#   intronic cage peaks are also contained.
#
#
# 07_deseq_cage_te.csv - deseq results of all called peaks that intersect with a
# TE
#
#


# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}

# --------------------------------- Functions ----------------------------------

getCondition <- function(header){
    
    # Returns a df with sample names (finally the row names of the df) and the
    # condition. The header needs to be separated by an underscore and the
    # condition needs to be at the second position. Usable for DESeq2.
    
    df <- data.frame(SampleID = header)
    df$condition <- sapply(df$SampleID, function(x) strsplit(x, "_")[[1]][2])
    
    rownames(df) <- df$SampleID
    df$SampleID <- NULL
    
    df$condition <- as.factor(df$condition)
    return(df)
}

bino.pValue <- function(x,n,p){
    # two sided is the default
    binom.test(x,n,p)$p.value
}


# ------------------------- Quantification of CAGE-peaks ----------------------- 

if (!file.exists('data/processed/cage_RS_deseq_segemehl.Rdata')) {
    
    data <- sapply(names(cage.files),
                   simplify = F,
                   function(x) getFeatureCountTab(cage.files[[x]], filter = F))
    
    count.tables <- sapply(names(data),
                           function(x) updateHeader(data[[x]][['counts']]))
    
    conditions <- sapply(names(count.tables),
                         simplify = F,
                         function(x) getCondition(names(count.tables[[x]])))
    
    # store peak annotation
    for (tissue in names(data)) {
        
        meta <-  data[[tissue]][['meta']]
        directory <- strsplit(cage.files[[tissue]], "/")[[1]][3]
        
        # set the length column to zero to mimic the 5th column in a bed file 
        meta$Length <- 0 
        
        write.table(meta, file = paste0('results/cage_RS/',
                                        directory,
                                        '/raw_peaks/',
                                        tissue,
                                        '.cage_peaks.bed'),
                    quote = F,
                    row.names = F,
                    col.names = F,
                    sep = '\t')
    }
    
    # runDESeq
    dds.cage <- sapply(names(count.tables), simplify = F,  function(x){
        
        doDEseq(count.matrix = count.tables[[x]],
                col_data = conditions[[x]],
                reference = 'young',
                target = 'all'
        )
    })
    
    # collect DESeq results
    deseq.cage <- sapply(names(dds.cage), simplify = F, function(x){
        
        getDEseqResults(dds.cage[[x]], coefficient = 'condition_old_vs_young', FDR.filter = F)
        
    })
    
    save(dds.cage, file = 'data/processed/cage_RS_dds_segemehl.Rdata')
    save(deseq.cage, file = 'data/processed/cage_RS_deseq_segemehl.Rdata')
    
    rm(dds.cage)
    
}else{
    
    deseq.cage <- loadRdata('data/processed/cage_RS_deseq_segemehl.Rdata')
}

deseq.cage.merged <- do.call('rbind',
                             sapply(names(deseq.cage),
                                    simplify = F,
                                    function(x) {
                                        df <- deseq.cage[[x]] %>%
                                            rownames_to_column(var = "peak_id") %>%
                                            mutate(tissue = x)
                                        
                                        return(df)
                                        
                                    }))

write.table(deseq.cage.merged, 
            file = paste0(table_dir, '07_deseq_cage_all.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


cageRanges <- sapply(names(deseq.cage), simplify = F, function(tissue){
    
    directory <- strsplit(cage.files[[tissue]], "/")[[1]][3]
    file <- paste0('results/cage_RS/',
           directory,
           '/raw_peaks/',
           tissue,
           '.cage_peaks.bed')
    
    #file = paste0('../data/cage/',tissue,'.peakAnnotation.bed')
    return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
    
})


# ------------------------- CAGE-peaks in Genes --------------------------------


gene.Cages <- sapply(names(cageRanges), simplify = F, function(x) {
    
    df <- intersectGranger(cageRanges[[x]], geneRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'cage',
                                              'subject' = 'gene'))
    
    return(df)
    
})


deseq.gene.cage.merged <- do.call('rbind', sapply(names(deseq.cage), simplify = F, function(x){
    
    df <- deseq.cage[[x]] %>%
        rownames_to_column(var = "peak_id") %>% 
        mutate(tissue = x)
    
    df <- merge(df, gene.Cages[[x]], by.x = 'peak_id', by.y = 'cage.names') 
    
    return(df)
    
}))

write.table(deseq.gene.cage.merged, 
            file = paste0(table_dir, '07_deseq_cage_gene.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


# ---------------------------- CAGE-peaks in TEs -------------------------------


te.Cages <- sapply(names(cageRanges), simplify = F, function(x) {
    
    df <- intersectGranger(cageRanges[[x]], teRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'cage',
                                              'subject' = 'te'))
    
    return(df)
    
})

deseq.te.cage.merged <- do.call('rbind',
                                sapply(names(deseq.cage),
                                       simplify = F,
                                       function(x) {
                                           df <- deseq.cage[[x]] %>%
                                               rownames_to_column(var = "peak_id") %>%
                                               mutate(tissue = x)
                                           
                                           df <-
                                               merge(df, te.Cages[[x]], by.x = 'peak_id', by.y = 'cage.names')
                                           
                                       }))

write.table(deseq.te.cage.merged, 
            file = paste0(table_dir, '07_deseq_cage_te.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


cage.TEs <- do.call('rbind', sapply(names(te.Cages), simplify = F, function(x){
    
    te.Cages[[x]] %>% 
        dplyr::select('te.te.id', 'cage.names') %>% 
        splitTEID("te.te.id") %>% 
        #filter(order %in% orders.of.interest) %>% 
        mutate(tissue = x) %>% 
        dplyr::rename(te.id = te.te.id,
                      peak_id = cage.names)
    
}))


independent.TEs <- sapply(c('brain', 'skin', 'blood'), function(x){
    
    cage.TEs %>% 
        filter(tissue == x) %>% 
        pull(te.id) %>% 
        unique()
    
})


te.super.enrichment <- do.call('rbind',
                               
    sapply(names(independent.TEs), simplify = F, function(x){
    
        df <- data.frame(TE.ID = unique(independent.TEs[[x]]), row.names = NULL) 
        df <- splitTEID(df, 'TE.ID')
        df <- binoRich(df, te.annotation, c(target = 'super_family',
                                            background = 'super_family'),
                       FDR.cap = 1e-10)
        df <- df %>% 
            filter(n.target >= 10) %>% 
            filter(!grepl("[?]", category))
        
        df$tissue <- x
        
        return(df)
    
}))

te.super.enrichment <- te.super.enrichment %>% 
    filter(category != 'NA') %>% 
    mutate(category = case_when(category == 'Alu' ~ 'B1',
                                TRUE ~ category))

te.super.enrichment$tissue <- 
    factor(te.super.enrichment$tissue, levels = c('brain', 'skin', 'blood'))

cat.Sort <- te.super.enrichment %>% 
    group_by(category) %>% 
    summarise(meanOdds = sum(ratio)) %>% 
    arrange(desc(meanOdds))


x.max <- max(te.super.enrichment$ratio) + 0.2

te.super.enrichment$category <-
    factor(te.super.enrichment$category, levels = cat.Sort$category)

label_pos <- te.super.enrichment %>% 
    group_by(tissue) %>% 
    dplyr::count() %>% pull(n) %>% max() + 2

cage.enrichment.pl <- ggplot(te.super.enrichment, aes(ratio, category, size = log10.padj)) +
    geom_point(aes(fill = tissue), shape = 21, , color = 'black') +
    scale_fill_manual(values = tissue.color) +
    scale_size(range = c(1, 5), name = expression(paste(log[10], "(FDR)"))) +
    annotate("text", x = -1, y = 0, label = "depleted", vjust = -1.5, color = 'white') + # trick to have a gap at the bottom
    geom_hline(yintercept = label_pos-1, linewidth = 0.4) + 
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    annotate("text", x = x.max/2, y = label_pos, label = "enriched", vjust = 1.2, size = 8/.pt) +
    annotate("text", x = -x.max/2, y = label_pos, label = "depleted", vjust = 1.2, size = 8/.pt) +
    theme_rob(base_size = 12, base_family = 'arial') +
    xlim(c(-x.max, x.max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    #facet_grid(category~., scale = 'free_y') +
    theme_rob(base_size = 8, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          #legend.position = 'bottom',
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.position = c(0.84, 0.65),
          legend.key.size = unit(0, "lines"),
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, size = 0.4),
          legend.box = "vertical",
          legend.box.just = 'right',) +
    guides(size = guide_legend(direction = 'vertical'))

cage.enrichment.pl

# Need the figure in 7x8 cm (wxh)
# 2.8x3.14 in
ggsave(cage.enrichment.pl,
       filename = paste0(figure_dir, '7_TSS_enrichment_bubble_new_style.pdf'),
       device = cairo_pdf,
       width = 7,
       height = 8,
       units = "cm",
       dpi = 300
)



#########
# RESIS #
#########

ggsave(cage.enrichment.pl,
       filename = './manuscripts/nature_aging/resis/figures/figure2/2A_TSS_enrichment_bubbl.pdf',
       device = cairo_pdf,
       width = 7,
       height = 8,
       units = "cm",
       dpi = 300
)


write.csv(te.super.enrichment, file = './manuscripts/nature_aging/resis/figures/figure2/2A_TSS_enrichment_bubbl.csv')

# ------------ Old representation --------


cage.enrichment.pl <- ggplot(te.super.enrichment, aes(ratio, category, color = tissue, size = log10.padj)) +
    geom_point() +
    scale_color_manual(values = tissue.color) +
    scale_size(range = c(3, 12), name = "log10(FDR)") +
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    xlim(c(-x.max, x.max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    #facet_grid(category~., scale = 'free_y') +
    theme_rob(base_size = 12, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          legend.position = 'bottom',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, size = 1.2),
          legend.box = "vertical")

cage.enrichment.pl


ggsave(
    filename = paste0(figure_dir, '07_TSS_enrichment_bubble.svg'),
    plot = last_plot(),
    width = 20,
    height = 30,
    units = "cm",
    dpi = 300)

write.table(te.super.enrichment, 
            file = paste0(table_dir, '07_TSS_enrichment_bubble.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

###### Family Enrichment stuff

te.family.enrichment <- do.call('rbind',
                                
                                sapply(names(independent.TEs), simplify = F, function(x){
                                    
                                    df <- data.frame(te_id = unique(independent.TEs[[x]]), row.names = NULL) 
                                    df <- splitTEID(df, 'te_id')
                                    df <- binoRich(df, te.annotation, c(target = 'family',
                                                                        background = 'family'),
                                                   FDR.cap = 1e-10)
                                    df <- df %>% 
                                        filter(n.target >= 10) %>% 
                                        filter(!grepl("[?]", category))
                                    
                                    df$tissue <- x
                                    
                                    return(df)
                                    
                                }))

te.family.enrichment <- te.family.enrichment %>% 
    filter(category != 'NA') %>% 
    mutate(category = case_when(category == 'Alu' ~ 'B1',
                                TRUE ~ category))

te.family.enrichment$tissue <- 
    factor(te.family.enrichment$tissue, levels = c('brain', 'skin', 'blood'))



cat.Sort <- te.family.enrichment %>% 
    group_by(category) %>% 
    summarise(meanOdds = sum(ratio)) %>% 
    arrange(desc(meanOdds))


x.max <- max(te.family.enrichment$ratio) + 0.2

te.family.enrichment$category <-
    factor(te.family.enrichment$category, levels = cat.Sort$category)


ggplot(te.family.enrichment %>% filter(ratio > 0.5), aes(ratio, category, color = tissue, size = log10.padj)) +
    geom_point() +
    scale_color_manual(values = tissue.color) +
    scale_size(range = c(3, 12), name = "log10(FDR)") +
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    xlim(c(-x.max, x.max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    #facet_grid(category~., scale = 'free_y') +
    theme_rob(base_size = 12, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          legend.position = 'bottom',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, size = 1.2),
          legend.box = "vertical")



#### Alternative presentation

label_pos <- te.super.enrichment %>% 
    group_by(tissue) %>% 
    dplyr::count() %>% pull(n) %>% max() + 2

ggplot(te.super.enrichment, aes(category, ratio, group = tissue)) +
    geom_line(aes(color = tissue),
              linetype = 1,
              linewidth = 1.1,
              alpha = 0.5) +
    geom_point(aes(alpha = log10.padj), size = 3 ) + 
    labs(x = 'TE family', 
         y = 'log(odds ratio)') + 
    coord_flip() +
    scale_color_manual(values = tissue.color) +
    ylim(c(-x.max, x.max)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = label_pos - 1) +
    annotate("text", y = x.max/2, x = label_pos, label = "enriched", vjust = 1.2, size = 12/.pt) +
    annotate("text", y = -x.max/2, x = label_pos, label = "depleted", vjust = 1.2, size = 12/.pt) +
    theme_rob(base_size = 12, base_family = 'arial') +
    scale_alpha(name = expression(paste(log[10], "(FDR)"))) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, size = 0.5),
          legend.position = c(0.86, 0.65),
          legend.direction = 'vertical',
          legend.box = "vertical",
          legend.box.just = 'right',
          #legend.key.size = unit(1, "lines"),
          legend.key = element_blank()
    )



# ---------------------------- CAGE-peaks in TE regions ------------------------


te_region_cages <- sapply(names(cageRanges), 
                          simplify = F, 
                          USE.NAMES = T, 
                          function(x){
                              
                              df <-  intersectGranger(cageRanges[[x]], teregionRanges, tab = 'all')
                              names(df) <- str_replace_all(names(df), c('query' = 'cage', 'subject' = 'te_region'))
                              
                              return(df)
                              
                          })


deseq.te.region.cage.merged <- do.call('rbind', sapply(names(deseq.cage), simplify = F, function(x){
    
    df <- deseq.cage[[x]] %>%
        rownames_to_column(var = "peak_id") %>% 
        mutate(tissue = x)
    
    
    df <- merge(df, te_region_cages[[x]], by.x = 'peak_id', by.y = 'cage.names')
    
    df <- merge(df, te_region_gene_relation, 
                by.x = 'te_region.names', 
                by.y = 'te_region_id',
                all.x = T)
    
    return(df)
    
}))


write.table(deseq.te.region.cage.merged, 
            file = paste0(table_dir, '07_deseq_cage_te_region.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


# ------------ Age correlation analysis ------------
#
# Potential Supplemental figures.
#### Enrichment analysis Family based

te.family.enrichment <- do.call('rbind',
                                
                                sapply(names(independent.TEs), simplify = F, function(x){
                                    
                                    df <- data.frame(TE.ID = unique(independent.TEs[[x]]), row.names = NULL)
                                    df <- splitTEID(df, 'TE.ID')
                                    df <- binoRich(df, te.annotation, c(target = 'family',
                                                                        background = 'family'),
                                                   FDR.cap = 1e-10)
                                    df <- df %>%
                                        filter(n.target >= 10) %>%
                                        filter(!grepl("[?]", category))
                                    
                                    df$tissue <- x
                                    
                                    return(df)
                                    
                                }))


alu_families <- te.annotation %>% 
    filter(super_family == 'Alu') %>% 
    pull(family)

l1_families <- te.annotation %>% 
    filter(super_family == 'L1') %>% 
    pull(family)

erv1_families <- te.annotation %>% 
    filter(super_family == 'ERV1') %>% 
    pull(family)


te.family.enrichment <- te.family.enrichment %>% 
    filter(category != 'NA') %>% 
    mutate(category = case_when(category == 'Alu' ~ 'B1',
                                TRUE ~ category))

te.family.enrichment$tissue <- 
    factor(te.family.enrichment$tissue, levels = c('brain', 'skin', 'blood'))


cat.Sort <- te.family.enrichment %>% 
    group_by(category) %>% 
    summarise(meanOdds = sum(ratio)) %>% 
    arrange(desc(meanOdds))


x.max <- max(abs(te.family.enrichment$ratio)) + 0.2

te.family.enrichment$category <-
    factor(te.family.enrichment$category, levels = cat.Sort$category)

label_pos <- te.family.enrichment %>% filter(category %in% l1_families) %>% 
    group_by(tissue) %>% 
    dplyr::count() %>% pull(n) %>% max() + 2

ggplot(te.family.enrichment %>% filter(category %in% l1_families), aes(ratio, category, size = log10.padj)) +
    geom_point(aes(fill = tissue), shape = 21, , color = 'black') +
    scale_fill_manual(values = tissue.color) +
    scale_size(range = c(1, 3), name = expression(paste(log[10], "(FDR)"))) +
    annotate("text", x = -1, y = 0, label = "depleted", vjust = -1.5, color = 'white') + # trick to have a gap at the bottom
    geom_hline(yintercept = 25, linewidth = 0.4) + 
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    annotate("text", x = x.max/2, y = label_pos, label = "enriched", vjust = 1.2, size = 8/.pt) +
    annotate("text", x = -x.max/2, y = label_pos, label = "depleted", vjust = 1.2, size = 8/.pt) +
    theme_rob(base_size = 12, base_family = 'arial') +
    xlim(c(-x.max, x.max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    #facet_grid(category~., scale = 'free_y') +
    theme_rob(base_size = 8, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          #legend.position = 'bottom',
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.position = c(0.84, 0.7),
          legend.key.size = unit(0, "lines"),
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, size = 0.4),
          legend.box = "vertical",
          legend.box.just = 'right',) +
    guides(size = guide_legend(direction = 'vertical'))


ggplot(te.family.enrichment %>% filter(category %in% alu_families), aes(ratio, category, size = log10.padj)) +
    geom_point(aes(fill = tissue), shape = 21, , color = 'black') +
    scale_fill_manual(values = tissue.color) +
    scale_size(range = c(1, 3), name = expression(paste(log[10], "(FDR)"))) +
    annotate("text", x = -1, y = 0, label = "depleted", vjust = -1.5, color = 'white') + # trick to have a gap at the bottom
    geom_hline(yintercept = 25, linewidth = 0.4) + 
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    # annotate("text", x = x.max/2, y = label_pos, label = "enriched", vjust = 1.2, size = 8/.pt) +
    # annotate("text", x = -x.max/2, y = label_pos, label = "depleted", vjust = 1.2, size = 8/.pt) +
    theme_rob(base_size = 12, base_family = 'arial') +
    xlim(c(-x.max, x.max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    #facet_grid(category~., scale = 'free_y') +
    theme_rob(base_size = 8, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          #legend.position = 'bottom',
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 12),
          legend.position = c(0.84, 0.7),
          legend.key.size = unit(0, "lines"),
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, size = 0.4),
          legend.box = "vertical",
          legend.box.just = 'right',) +
    guides(size = guide_legend(direction = 'vertical'))


l1_family <- te.annotation %>% 
    filter(super_family == "L1") %>% 
    mutate(Kimura = as.numeric(Kimura))


line_as_bg <- te.annotation %>% 
    filter(order == 'LINE') %>% 
    mutate(Kimura = as.numeric(Kimura),
           family = 'background')


kim_df <- rbind(l1_family, line_as_bg) 

family_order <- kim_df %>% 
    group_by(family) %>% 
    dplyr::summarise(kimura_mean = mean(as.numeric(Kimura), na.rm = TRUE), .groups = 'drop') %>% 
    arrange(kimura_mean) %>% 
    pull(family)

kim_df$family <- factor(kim_df$family, levels = family_order)

kimura_pl <- ggplot(kim_df, aes(family, as.numeric(Kimura), fill = family)) + 
    geom_violin() +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6) +
    labs(y = "Kimura distance") +
    coord_flip() +
    theme_rob() +
    theme(legend.position = 'None')



#### correlation

kim_mean <- kim_df %>% 
    group_by(family) %>% 
    summarise(kim_age = mean(Kimura))

enriched_l1s <- te.family.enrichment %>% filter(category %in% kim_mean$family)#l1_family)

#foo <- merge(kim_df, te.family.enrichment %>% filter(category %in% alu_families), by.y = 'category', by.x = 'family')

enriched_l1s <- merge(enriched_l1s, kim_mean, by.x = 'category', by.y = 'family')

ggplot(enriched_l1s, aes(kim_age, ratio)) +
    geom_point(aes(color = tissue)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq() +
    scale_color_manual(values = tissue.color[2:4])

young <- enriched_l1s %>% 
    filter(kim_age < 30) %>% 
    ggplot(aes(category)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90))

old <- enriched_l1s %>% 
    filter(kim_age > 30) %>% 
    ggplot(aes(category)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90))



erv1_family <- te.annotation %>% 
    filter(super_family == "ERV1") %>% 
    mutate(Kimura = as.numeric(Kimura))

LTR_as_bg <- te.annotation %>% 
    filter(order == 'LTR') %>% 
    mutate(Kimura = as.numeric(Kimura),
           family = 'background')


kim_df <- rbind(erv1_family, LTR_as_bg) 

family_order <- kim_df %>% 
    group_by(family) %>% 
    dplyr::summarise(kimura_mean = mean(as.numeric(Kimura), na.rm = TRUE), .groups = 'drop') %>% 
    arrange(kimura_mean) %>% 
    pull(family)

kim_df$family <- factor(kim_df$family, levels = family_order)

kimura_pl <- ggplot(kim_df, aes(family, as.numeric(Kimura), fill = family)) + 
    geom_violin() +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6) +
    labs(y = "Kimura distance") +
    coord_flip() +
    theme_rob() +
    theme(legend.position = 'None')



#### correlation

kim_mean <- kim_df %>% 
    group_by(family) %>% 
    summarise(kim_age = mean(Kimura))

enriched_erv1s <- te.family.enrichment %>% filter(category %in% kim_mean$family)#l1_family)

#foo <- merge(kim_df, te.family.enrichment %>% filter(category %in% alu_families), by.y = 'category', by.x = 'family')

enriched_erv1s <- merge(enriched_erv1s, kim_mean, by.x = 'category', by.y = 'family')

ggplot(enriched_erv1s, aes(kim_age, ratio)) +
    geom_point(aes(color = tissue)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq() +
    scale_color_manual(values = tissue.color[2:4])

