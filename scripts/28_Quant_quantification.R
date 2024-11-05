# --------------------------------- Notes--- -----------------------------------
# 
# This script is used for the Quant-Seq quantification
#
# Output data (../data/quant; ../tables):
#
# [tissue].quant.peakAnnotation - bedfile of all called peaks
#
# dds.quant.Rdata - dds object for called quant peaks (peaks with less than 10 
#    reads across all samples are removed) 
#
# deseq.quant.Rdata - deseq results for all called peaks 
#
# 28_deseq_quant_all.csv - deseq results of all called peaks
#
# 28_deseq_quant_gene.csv - deseq results of all called peaks that intersect with 
#   a gene. The longest version of the gene was used for that intersection and
#   intronic cage peaks are also contained.
#
#
# 28_deseq_quant_te.csv - deseq results of all called peaks that intersect with a
# TE

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

# ------------------------- Quantification of Quant-peaks ----------------------- 

if (!file.exists('data/processed/quant_RS_deseq_segemehl.Rdata')) {
    
    data <- sapply(names(quant.files),
                   simplify = F,
                   function(x) getFeatureCountTab(quant.files[[x]], filter = F))
    
    count.tables <- sapply(names(data),
                           function(x) updateHeader(data[[x]][['counts']]))
    
    conditions <- sapply(names(count.tables),
                         simplify = F,
                         function(x) getCondition(names(count.tables[[x]])))
    
    # store peak annotation
    for (tissue in names(data)) {
        
        meta <-  data[[tissue]][['meta']]
        directory <- strsplit(quant.files[[tissue]], "/")[[1]][3]
        
        # set the length column to zero to mimic the 5th column in a bed file 
        meta$Length <- 0 
        
        write.table(meta, file = paste0('results/quant_RS/',
                                        directory,
                                        '/raw_peaks/',
                                        tissue,
                                        '.quant_peaks.bed'),
                    quote = F,
                    row.names = F,
                    col.names = F,
                    sep = '\t')
    }
    
    # runDESeq
    dds.quant <- sapply(names(count.tables), simplify = F,  function(x){
        
        doDEseq(count.matrix = count.tables[[x]],
                col_data = conditions[[x]],
                reference = 'young',
                target = 'all'
        )
    })
    
    # collect DESeq results
    deseq.quant <- sapply(names(dds.quant), simplify = F, function(x){
        
        getDEseqResults(dds.quant[[x]], coefficient = 'condition_old_vs_young', FDR.filter = F)
        
    })
    
    save(dds.quant, file = 'data/processed/quant_RS_dds_segemehl.Rdata')
    save(deseq.quant, file = 'data/processed/quant_RS_deseq_segemehl.Rdata')
    
    rm(dds.quant)
    
}else{
    
    deseq.quant <- loadRdata('data/processed/quant_RS_deseq_segemehl.Rdata')
}

deseq.quant.merged <- do.call('rbind',
                             sapply(names(deseq.quant),
                                    simplify = F,
                                    function(x) {
                                        df <- deseq.quant[[x]] %>%
                                            rownames_to_column(var = "peak_id") %>%
                                            mutate(tissue = x)
                                        
                                        return(df)
                                        
                                    }))

write.table(deseq.quant.merged, 
            file = paste0(table_dir, '28_deseq_quant_all.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


quantRanges <- sapply(names(deseq.quant), simplify = F, function(tissue){
    
    directory <- strsplit(quant.files[[tissue]], "/")[[1]][3]
    file <- paste0('results/quant_RS/',
                   directory,
                   '/raw_peaks/',
                   tissue,
                   '.quant_peaks.bed')
    
    #file = paste0('../data/cage/',tissue,'.peakAnnotation.bed')
    return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
    
})

# ------------------------- Quant-peaks in Genes --------------------------------


gene.quants <- sapply(names(quantRanges), simplify = F, function(x) {
    
    df <- intersectGranger(quantRanges[[x]], geneRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'quant',
                                              'subject' = 'gene'))
    
    return(df)
    
})


deseq.gene.quant.merged <- do.call('rbind', sapply(names(deseq.quant), simplify = F, function(x){
    
    df <- deseq.quant[[x]] %>%
        rownames_to_column(var = "peak_id") %>% 
        mutate(tissue = x)
    
    df <- merge(df, gene.quants[[x]], by.x = 'peak_id', by.y = 'quant.names') 
    
    return(df)
    
}))

write.table(deseq.gene.quant.merged, 
            file = paste0(table_dir, '28_deseq_quant_gene.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


# ---------------------------- Quant-peaks in TEs -------------------------------


te.quants <- sapply(names(quantRanges), simplify = F, function(x) {
    
    df <- intersectGranger(quantRanges[[x]], teRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'quant',
                                              'subject' = 'te'))
    
    return(df)
    
})

deseq.te.quant.merged <- do.call('rbind',
                                sapply(names(deseq.quant),
                                       simplify = F,
                                       function(x) {
                                           df <- deseq.quant[[x]] %>%
                                               rownames_to_column(var = "peak_id") %>%
                                               mutate(tissue = x)
                                           
                                           df <-
                                               merge(df, te.quants[[x]], by.x = 'peak_id', by.y = 'quant.names')
                                           
                                       }))



write.table(deseq.te.quant.merged, 
            file = paste0(table_dir, '28_deseq_quant_te.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

##### Position Composition ###

composition_df <- deseq.te.quant.merged %>% 
    filter(!is.na(padj)) %>% 
    group_by(tissue) %>% 
    dplyr::count(te.position) %>% 
    mutate(per = n/sum(n)*100,
           label = paste0(round(per,1), " %")) %>% 
    ungroup() 

composition_pl <- ggplot(composition_df, aes(tissue, per, fill = te.position)) +
    geom_col() +
    scale_fill_brewer(palette = "Dark2",
                      name = "position") +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 6/.pt) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    ylab("qTE position in genomic context [%]") +
    theme_rob(base_size = 8) +
    theme(axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.text = element_text(size = 8))


ggsave(composition_pl,
       filename = paste0(figure_dir, '28_TTS_TE_composition.png'),
       device = "png",
       width = 10,
       height = 10,
       units = "cm",
       dpi = 300
)


#### Enrichment (bubble ) plot 

quant.TEs <- do.call('rbind', sapply(names(te.quants), simplify = F, function(x){
    
    te.quants[[x]] %>% 
        dplyr::select('te.te.id', 'quant.names') %>% 
        splitTEID("te.te.id") %>% 
        #filter(order %in% orders.of.interest) %>% 
        mutate(tissue = x) %>% 
        dplyr::rename(te_id = te.te.id,
                      peak_id = quant.names)
    
}))


finisher.TEs <- sapply(c('brain', 'skin', 'blood'), function(x){
    
    quant.TEs %>% 
        filter(tissue == x) %>% 
        pull(te_id) %>% 
        unique()
    
})

te.super.enrichment <- do.call('rbind',
                               
                               sapply(names(finisher.TEs), simplify = F, function(x){
                                   
                                   df <- data.frame(te_id = unique(finisher.TEs[[x]]), row.names = NULL) 
                                   df <- splitTEID(df, 'te_id')
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


save(te.super.enrichment, file = paste0(data_dir, "28_te_super_fam_quant_enrichment.Rdata"))

cat.Sort <- te.super.enrichment %>% 
    group_by(category) %>% 
    summarise(meanOdds = sum(ratio)) %>% 
    arrange(desc(meanOdds))


x.max <- max(abs(te.super.enrichment$ratio)) + 0.2

te.super.enrichment$category <-
    factor(te.super.enrichment$category, levels = cat.Sort$category)

label_pos <- te.super.enrichment %>% 
    group_by(tissue) %>% 
    dplyr::count() %>% pull(n) %>% max() + 2

quant.enrichment.pl <- ggplot(te.super.enrichment, aes(ratio, category, size = log10.padj)) +
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

quant.enrichment.pl

# Need the figure in 7x8 cm (wxh)
# 2.8x3.14 in
ggsave(quant.enrichment.pl,
       filename = paste0(figure_dir, '28_TTS_enrichment_bubble.pdf'),
       device = cairo_pdf,
       width = 7,
       height = 8,
       units = "cm",
       dpi = 300
)


ggsave(quant.enrichment.pl,
       filename = paste0(figure_dir, '28_TTS_enrichment_bubble.png'),
       device = "png",
       width = 8,
       height = 10,
       units = "cm",
       dpi = 300
)

ggsave(
    filename = paste0(figure_dir, '28_TTS_enrichment_bubble.svg'),
    plot = last_plot(),
    width = 20,
    height = 30,
    units = "cm",
    dpi = 300)

write.table(te.super.enrichment, 
            file = paste0(table_dir, '28_TTS_enrichment_bubble.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

####### Alternative Plot

ggplot(te.super.enrichment, aes(category, ratio, group = tissue)) +
    geom_line(aes(color = tissue),
              linewidth = 1.3,
              alpha = 0.5) +
    geom_point(aes(alpha = log10.padj), size = 3 ) + 
    coord_flip() +
    scale_color_manual(values = tissue.color) +
    ylim(c(-x.max, x.max)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 25) +
    annotate("text", y = 1, x = 26, label = "enriched", vjust = 1.2, size = 12/.pt) +
    annotate("text", y = -1, x = 26, label = "depleted", vjust = 1.2, size = 12/.pt) +
    theme_rob(base_size = 12, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          legend.position = c(0.9, 0.65),
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, size = 0.5),
          legend.box = "vertical")


#### Enrichment analysis Family based

te.family.enrichment <- do.call('rbind',

                               sapply(names(finisher.TEs), simplify = F, function(x){

                                   df <- data.frame(TE.ID = unique(finisher.TEs[[x]]), row.names = NULL)
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



te.family.enrichment <- te.family.enrichment %>% 
    filter(category != 'NA') %>% 
    mutate(category = case_when(category == 'Alu' ~ 'B1',
                                TRUE ~ category))

te.family.enrichment$tissue <- 
    factor(te.family.enrichment$tissue, levels = c('brain', 'skin', 'blood'))

save(te.family.enrichment, file = paste0(data_dir, "28_te_family_quant_enrichment.Rdata"))


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

ggplot(te.family.enrichment %>% filter(category %in% l1_families, ratio > 0), aes(ratio, category, size = log10.padj)) +
    geom_point(aes(fill = tissue), shape = 21, , color = 'black') +
    scale_fill_manual(values = tissue.color) +
    scale_size(range = c(1, 3), name = expression(paste(log[10], "(FDR)"))) +
    geom_hline(yintercept = label_pos - 1, linewidth = 0.4) + 
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    theme_rob(base_size = 12, base_family = 'arial') +
    xlim(c(-x.max, x.max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    facet_grid(category~., scale = 'free_y') +
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


id_families <- te.annotation %>% filter(super_family == 'ID') %>% pull(family)

ggplot(te.family.enrichment %>% filter(category %in% id_families), aes(ratio, category, size = log10.padj)) +
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


#### Kimura distance stuff


b1_family <- te.annotation %>% 
    filter(super_family == "Alu") %>% 
    mutate(Kimura = as.numeric(Kimura))


sine_as_bg <- te.annotation %>% 
    filter(order == 'SINE') %>% 
    mutate(Kimura = as.numeric(Kimura),
           family = 'background')


kim_df <- rbind(b1_family, sine_as_bg) 

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

enriched_b1s <- te.family.enrichment %>% filter(category %in% alu_families)

#foo <- merge(kim_df, te.family.enrichment %>% filter(category %in% alu_families), by.y = 'category', by.x = 'family')

enriched_b1s <- merge(enriched_b1s, kim_mean, by.x = 'category', by.y = 'family')

ggplot(enriched_b1s, aes(kim_age, ratio)) +
    geom_point(aes(color = tissue)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq() +
    scale_color_manual(values = tissue.color[2:4])

# ---------------------------- Quant-peaks in TEs Lollipop -------------------------------

############################################ TE dot chart ####################################################

# Counts the entries of a data frame filled with TEs. Can count order, super_family, or order, depending on the
# set te_level. Features with less then 5 counts are removed by default, but can set by min_count.

# data preparations 
prep_te_count_plot <- function(df, te_level = 'super_family', min_count = 5){
    # takes a data frame with a te_id column and counts the TE entries respective of the 
    # selected TE level.
    
    require(tidyverse) 
    
    df <- blackRcloud::splitTEID(df, 'te_id')
    
    te_assignment <- df %>% 
        dplyr::select(order, super_family, family) %>% 
        filter(!duplicated(family))
    
    df_overlapp <- df %>% 
        #filter(distance == 0) %>% 
        group_by(!!sym(te_level)) %>% 
        summarise(count = n())
    
    df_overlapp <- merge(df_overlapp,
                         te_assignment,
                         by = te_level)
    
    df_overlapp <- df_overlapp %>% 
        filter(count >= min_count,
               te_level != 'NA')
    
    return(df_overlapp)
    
}

prep_te_proportion_plot <- function(df,
                                    member_count_file,
                                    te_level = 'super_family',
                                    min_count = 5
                                    ){
    
    # Calculates the proportion of each member to the whole order/super_family/family.
    
    
    require(tidyverse)
    
    df <- prep_te_count_plot(df, te_level, min_count) 
    
    te_member_counts <- member_count_file[[te_level]]
    
    
    df <- merge(df, te_member_counts, by = te_level)
    
    df$proportion <- df$count/df$total_count  
    
    return(df)
    
}

# plot
dotChart <- function(df,
                     te_level = 'super_family',
                     y_value = 'count', # can be count or proportion
                     color_palette = c("#00AFBB", "#e9c46a", "#FC4E07", "#264653", "#2a9d8f")){
    
    # Keep the total counts in the bubbles, even when you plot the proportion
    require(ggpubr)
    
    pl <- ggdotchart(df, x = te_level, y = y_value,
                     color = "order",                                # Color by groups
                     palette = color_palette, # Custom color palette
                     sorting = "descending",                       # Sort value in descending order
                     add = "segments",                             # Add segments from y = 0 to dots
                     rotate = FALSE,                                # Rotate vertically
                     group = "order",                              # Order by groups
                     dot.size = 10,                                # Large dot size
                     label = df$count,                             # Add count values as dot labels
                     font.label = list(color = "white", size = 9, 
                                       vjust = 0.5, face = 'plain'),               # Adjust label parameters
                     ggtheme = theme_pubr()                        # ggplot2 theme
    )
    
    return(pl)
    
    
}

foo <- te.quants[['brain']]

foo <- foo %>% rename(te_id = te.te.id)

member_counts <- list("order" = te.annotation %>% dplyr::count(order, name = "total_count"),
                   "super_family" = te.annotation %>% dplyr::count(super_family, name = "total_count"),
                   "family" = te.annotation %>% dplyr::count(family, name = "total_count")
                   )

df_prop <- prep_te_proportion_plot(foo, member_counts, te_level = 'super_family', min_count = 10) %>% 
    filter(order %in% orders.of.interest, super_family != 'NA')


plot_h = round(nrow(unique(df_prop["super_family"])) * 0.25)


dotChart(df_prop, te_level = "super_family", y_value = 'count')


dotChart(df_prop, te_level = "super_family", y_value = 'proportion')

