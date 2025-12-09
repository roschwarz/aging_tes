# Categorizing L1 3'ends as a huge set pops up in brain by overlapping TEs with CAGE-peaks

# ================================== Setup ================================================

# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

load_cage_peak_annotation()
load_te_ranges()
load_te_annotation()
aging_tes::load_plotting_env()


# chr1|-|3264193|3264830|LINE|L1|Lx_3end|323|11.00
# ================================= Data =================================================

# Identify transposable elements with their own TSS
# Overlap TEs with CAGE peaks for each tissue
# Subsequently run analyze_te_statistics()
# remove duplicated TEs to avoid biasing the Kimura distribution
indie_TEs <- sapply(names(cageRanges), simplify = FALSE, function(tissue){
    
    cage_peaks <- cageRanges[[tissue]]
    
    overlap_result <- process_overlapping(
        query_ranges = teRanges,
        subject_ranges = cage_peaks,
        col_query_to_subj = "te.id"
    )
    
    df_te <- data.frame(te_id = overlap_result$query_with_hit$te.id,
                        ensembl_gene_id = mcols(overlap_result$query_with_hit)$ensembl_gene_id) %>%
        blackRcloud::splitTEID("te_id") %>%
        filter(order %in% orders.of.interest, !duplicated(te_id))
    
    te_stat_results <- analyze_te_statistics(df_te, entity_col = "te_id", colors = order.color)
    
    res <- list(statistics = te_stat_results,
                overlap_res = overlap_result)
    return(res)
})

#### Combine figures

kimura_tes_w_cage_p <- ggarrange(indie_TEs$brain$statistics$plots$kimura_distribution_histogram +
                                       ggtitle("brain") +
                                       xlim(0,50) +
                                       theme(plot.title = element_text(hjust = 0.5)),
                                   indie_TEs$skin$statistics$plots$kimura_distribution_histogram +
                                       ggtitle("skin") +
                                       xlim(0,50) +
                                       theme(plot.title = element_text(hjust = 0.5)),
                                   indie_TEs$skin$statistics$plots$kimura_distribution_histogram +
                                       ggtitle("blood") +
                                       xlim(0,50) +
                                       theme(plot.title = element_text(hjust = 0.5)),
                                   ncol = 3, nrow = 1,
                                   common.legend = TRUE,
                                   legend = "bottom",
                                   labels = c("A", "B", "C"))

kimura_tes_w_cage_p1 <- annotate_figure(kimura_tes_w_cage_p,
                top = text_grob("Age distribution of TEs overlapping with CAGE-peaks", size = 14))


meta <- list(name = 'age_distribution_tes_with_cage_peaks',
             description = 'Age distribution of TEs overlapping with CAGE peaks. Duplicated TEs are removed',
             tags = c('CAGE-Seq', 'TEs', 'Kimura'),
             parameters = list(tissues = c('brain', 'skin', 'blood')),
             script = 'L1_3_end_analysis.R'
)

fig_index(plot = kimura_tes_w_cage_p1,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 18,
          height = 6,
          dpi = 300,
          format = 'pdf')

# =======================================================================================
# Let's look for 3end L1 elements

# Identify transposable elements with their own TSS
# Overlap TEs with CAGE peaks for each tissue
# Subsequently run analyze_te_statistics()
indie_L1_3end_TEs <- sapply(names(cageRanges), simplify = FALSE, function(tissue){
    
    cage_peaks <- cageRanges[[tissue]]
    
    overlap_result <- process_overlapping(
        query_ranges = teRanges,
        subject_ranges = cage_peaks,
        col_query_to_subj = "te.id"
    )
    
    df_te <- data.frame(te_id = overlap_result$query_with_hit$te.id,
                        ensembl_gene_id = mcols(overlap_result$query_with_hit)$ensembl_gene_id) %>%
        blackRcloud::splitTEID("te_id") %>%
        filter(order %in% orders.of.interest, !duplicated(te_id), order == 'LINE', grepl('3end', family))
    
    te_stat_results <- analyze_te_statistics(df_te, entity_col = "te_id", colors = order.color)
    
    res <- list(statistics = te_stat_results,
                overlap_result)
    return(res)
})

#### Combine figures

kimura_l1_3end_w_cage_p <- ggarrange(indie_L1_3end_TEs$brain$statistics$plots$kimura_distribution_histogram +
                                     ggtitle("brain") +
                                     scale_fill_manual(values = 'black') +
                                     xlim(0,50) +
                                     theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = 'none'),
                                 indie_L1_3end_TEs$skin$statistics$plots$kimura_distribution_histogram +
                                     ggtitle("skin") +
                                     scale_fill_manual(values = 'black') +
                                     xlim(0,50) +
                                     theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = 'none'),
                                 indie_L1_3end_TEs$skin$statistics$plots$kimura_distribution_histogram +
                                     ggtitle("blood") +
                                     scale_fill_manual(values = 'black') +
                                     xlim(0,50) +
                                     theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = 'none'),
                                 ncol = 3, nrow = 1,
                                 common.legend = FALSE,
                                 legend = NULL,
                                 labels = c("A", "B", "C"))

kimura_l1_3end_w_cage_p1 <- annotate_figure(kimura_l1_3end_w_cage_p,
                                        top = text_grob("Age distribution of L1-3ends overlapping with CAGE-peaks", size = 14))


meta <- list(name = 'age_distribution_l1_3end_with_cage_peaks',
             description = 'Age distribution of L1 3-ends overlapping with CAGE-peaks',
             tags = c('CAGE-Seq', 'L1', 'Kimura'),
             parameters = list(tissues = c('brain', 'skin', 'blood')),
             script = 'L1_3_end_analysis.R'
)

fig_index(plot = kimura_l1_3end_w_cage_p1,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 18,
          height = 6,
          dpi = 300,
          format = 'pdf')


# ============================= Analyze =================================================

# Intergenetic TEs with CAGE peak

df_te <- data.frame(te_id = indie_TEs$brain$overlap_res$query_with_hit$te.id,
                    ensembl_gene_id = mcols(indie_TEs$brain$overlap_res$query_with_hit)$ensembl_gene_id) %>% 
    filter(is.na(ensembl_gene_id)) %>%
    blackRcloud::splitTEID("te_id") %>%
    filter(order %in% orders.of.interest)


te_stat_results <- analyze_te_statistics(df_te, entity_col = "te_id", colors = order.color)

te_stat_results$plots$kimura_distribution_histogram

# Intragenetic TEs with CAGE peak

df_te <- data.frame(te_id = indie_TEs$brain$query_with_hit$te.id,
                    ensembl_gene_id = mcols(indie_TEs$brain$query_with_hit)$ensembl_gene_id) %>% 
    filter(!is.na(ensembl_gene_id)) %>%
    blackRcloud::splitTEID("te_id") %>%
    filter(order %in% orders.of.interest)


te_stat_results <- analyze_te_statistics(df_te, entity_col = "te_id", colors = order.color)

te_stat_results$plots$kimura_distribution_histogram

# ========================= Combine data for all tissues ================================


indie_TEs <- do.call('rbind', lapply(names(cageRanges), function(tissue) {
    cage_peaks <- cageRanges[[tissue]]
    df <- blackRcloud::intersectGranger(cage_peaks, teRanges, tab = 'subject') %>% data.frame() %>%
        filter(duplicated(te.id)) %>% 
        blackRcloud::splitTEID("te.id") %>%
        #filter(grepl('3end', family)) %>%
        dplyr::select(te.id, chromosome, strand, start, end, order, super_family, family, id, Kimura, position, ensembl_gene_id, external_gene_name) %>% 
        dplyr::rename(te_id = te.id) %>% 
        mutate(tissue = tissue)
}))



background <- te_annotation %>% 
#        filter(grepl('3end', family)) %>%
        mutate(tissue = 'background')

indie_TEs <- rbind(indie_TEs, background)

indie_TEs$width <- as.numeric(indie_TEs$end) - as.numeric(indie_TEs$start)
indie_TEs$Kimura <- as.numeric(indie_TEs$Kimura)

# =============================== Plots =========================================

# =================== Position of 3end LINE elements ============================

position_count <- indie_TEs %>% 
    group_by(tissue, position) %>% 
    summarise(n = n()) %>% 
    mutate(prop = n/sum(n) * 100)

position_count$tissue <- factor(position_count$tissue, levels = c('background', 'brain', 'skin', 'blood'))

count_text <- position_count %>% 
    group_by(tissue) %>% 
    summarise(total = sum(n)) %>% 
    mutate(text = paste0("n = ", format(total, big.mark = ',')),
           prop = 100)

ggplot(position_count, aes(tissue, prop)) +
    geom_col(aes(fill = position)) +
    geom_text_repel(data = count_text, aes(label = text),
                    size = 8/.pt,
                    vjust = -2.5) +
    labs(y = "Percent of (expressed) L1 3'ends") +
    scale_fill_manual(values = position_color, name = "Genomic location:") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 105)) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          legend.position = 'bottom')

# =================== Age of L1 3'ends (Kimura distance) ========================

breaks <- seq(0, 80, by = 5)

indie_TEs$kimura_group <- cut(indie_TEs$Kimura, breaks = breaks, include.lowest = TRUE, right = FALSE)

indie_TEs_kimura_group <- indie_TEs %>% 
    group_by(tissue, kimura_group) %>% 
    summarise(n = n()) %>% 
    mutate(prop = n/sum(n) * 100)

indie_TEs_kimura_group$tissue <- factor(indie_TEs_kimura_group$tissue,
                                        levels = c('background', 'brain', 'skin', 'blood'))

# Filtered for Age group with a proportion bigger than 0.005    
ggplot(indie_TEs_kimura_group %>% filter(prop >= 0.005), aes(kimura_group, prop, fill = tissue)) +
    geom_col() +
    scale_fill_manual(values = tissue.color) +
    facet_grid(tissue~.) +
    labs(x = "Age group (Kimura distance)",
         y = "Percent of (expressed) L1 3'ends") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
        legend.position = 'bottom')

# =================== Size of L1 3'ends (base pairs) ============================

indie_TEs$tissue <- factor(indie_TEs$tissue, levels = rev(c('background', 'brain', 'skin', 'blood')))

ggplot(indie_TEs, aes(tissue, width, fill = tissue)) +
    geom_violin() +
    scale_fill_manual(values = tissue.color) +
    scale_y_log10() +
    # scale_y_continuous(expand = c(0,0),
                       # limits = c(0, (max(indie_TEs$width)+150))) +
    labs(y = "Size in bp") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = 'None',
          axis.title.y = element_blank())

# ============================= Development =====================================


foo <- te_annotation %>% 
    filter(grepl("L1MdV_", family)) %>% 
    separate(family, into = c("fam", "sub_fam", "domain"), 
             sep = "_", remove = FALSE) %>% 
    mutate(family = paste0(fam, "_", sub_fam))

foo_fam_size <- foo %>% 
    group_by(family) %>% 
    summarise(n = n())

foo <- merge(foo, foo_fam_size, by = 'family')

ggplot(foo, aes(n, as.numeric(Kimura))) +
    geom_point()

anno <- te_annotation %>% 
    separate(family, into = c("fam", "sub_fam", "domain"), 
             sep = "_", remove = FALSE) %>% 
    mutate(family_updated = paste0(fam, "_", sub_fam)) %>% 
    group_by(family) %>% 
    summarise(n = n(),
              kimura_mean = mean(as.numeric(Kimura)))





ggplot(indie_TEs %>% filter(tissue != 'background'), aes(position, fill = tissue)) +
    geom_bar(position = 'dodge') +
    coord_flip()

ggplot(indie_TEs, aes(family, fill = tissue)) +
    geom_bar(position = 'dodge') +
    coord_flip()


ggplot(indie_TEs, aes(Kimura, color = family)) +
    geom_density() +
    facet_grid(tissue~.)

ggplot(indie_TEs, aes(width, color = tissue)) +
    geom_density()


ggplot(indie_TEs, aes(width, Kimura, color = tissue, alpha = 0.05)) +
    geom_point() +
    geom_smooth(method = 'lm')

background_count <- background %>% group_by(family) %>% summarise(background_n = n())
indie_TE_count <- indie_TEs %>% group_by(tissue, family) %>% summarise(n = n())

indie_TEs <- merge(indie_TEs, background_count, by = 'family', all.x = TRUE)

indie_TE_count <- merge(indie_TE_count, background_count, by = 'family', all.x = TRUE)

indie_TE_count$bg_proportion <- indie_TE_count$n/indie_TE_count$background_n*100

ggplot(indie_TE_count %>% filter(tissue != 'background'), aes(family, bg_proportion, fill = tissue)) +
    geom_col(position = 'dodge') +
    coord_flip()


indie_TE_count %>% 
    filter(tissue != 'background') %>% 
    group_by(family) %>% 
    summarise(total = sum(n)) %>% view()
