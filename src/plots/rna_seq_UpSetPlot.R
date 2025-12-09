if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_env()
aging_tes::load_plotting_env()

# ------------------------------------------------------------------------------
# VennDiagrams for all TEs and tissues in male and female
# ------------------------------------------------------------------------------

deseq_te_merged_male <- fread(paste0(table_dir, deseq_results_te_csv_male)) %>% 
    mutate(sex = 'male',
           sex_tissue = paste0(sex, "_",tissue))

deseq_te_merged_female <- fread(paste0(table_dir, deseq_results_te_csv_female)) %>% 
    mutate(sex = 'female',
           sex_tissue = paste0(sex, "_", tissue))

deseq_tes <- rbind(deseq_te_merged_male, deseq_te_merged_female)

expressed_TEs <- sapply(unique(deseq_tes$sex_tissue), function(x){
    
    deseq_tes %>% 
        filter(sex_tissue == x, !is.na(padj)) %>% 
        pull(te_id) %>% 
        unique()
})

up_set_pl <- upset(fromList(expressed_TEs), 
      nsets = length(unique(deseq_tes$sex_tissue)), 
      nintersects = 20,
      order.by = "freq",
      sets.bar.color = tissue.color[c(2,3,4,2,3)],
      main.bar.color = "gray23",
      matrix.color = "gray23",
      point.size = 1.5,
      line.size = 1,
      text.scale = c(1.2, 1, 1.2, 1,
                     1.2, 1),
      mb.ratio = c(0.6, 0.4),
      sets.x.label = "Number of expressed TEs",
      mainbar.y.label = "Number of TEs in intersection",
      keep.order = TRUE
)



meta <- list(name = 'te_instances_upset',
             description = 'UpSet plot of all expressed transposable elements in all tissues and sexes',
             tags = c('expression', 'rna-seq', 'TE'),
             parameters = list(expressed = "!NA", tissues = c('brain', 'skin', 'blood'), sex = c('male', 'female')),
             script = 'rna_seq_UpSetPlot.R')

fig_index(plot = up_set_pl,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 20,
          height = 10,
          dpi = 300,
          format = 'pdf')

# ============================================================================================================
# Just to get the overlap between all conditions.
# ============================================================================================================


library(VennDiagram)

venn <- venn.diagram(
    x = expressed_TEs,
    category.names = names(expressed_TEs),
    
    
    # Circles
    lwd = 2,  
    #fill = tissue.color[2:4], #c('#264653', '#2A9D8F',  '#E9C46A'),
    alpha = c(0.7, 0.7, 0.7, 0.7, 0.7),
    
    
    # Number
    cex = 1, # font size
    fontface = "bold",
    fontfamily = "arial",
    # 
    # # Set names
    cat.cex = 1.5,
    cat.default.pos = "outer",
    cat.fontface = "bold",
    cat.fontfamily = "arial",
    # cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055, 0.055, 0.055, 0.055),
    # main = header,
    scaled = F,
    print.mode = c("raw", "percent"),
    
    # Output
    filename = NULL, #paste0(figures, 'Panel_2C_VennDiagram.svg'),
    imagetype = "svg",
    output = FALSE,
    width = 200,
    height = 500,
    resolution = 300,
    disable.logging = TRUE
    
)

grid::grid.draw(venn)



# -----------------------------------------------------------------------------------------------------------
# UpSet-plot for DETEs
# You can adapt the adjusted p-value to play around
# -----------------------------------------------------------------------------------------------------------

deTEs <- sapply(unique(deseq_tes$sex_tissue), function(x){
    
    deseq_tes %>% 
        filter(sex_tissue == x, padj <= 0.05) %>% 
        pull(te_id) %>% 
        unique()
})

upset(
    fromList(deTEs),
    nsets = length(unique(deseq_tes$sex_tissue)),
    nintersects = 20,
    order.by = "freq",
    sets.bar.color = tissue.color[c(2, 3, 4, 2, 3)],
    main.bar.color = "gray23",
    matrix.color = "gray23",
    point.size = 1.5,
    line.size = 1,
    text.scale = c(2, 1.5, 2, 1.5, 2, 1.5),
    mb.ratio = c(0.6, 0.4),
    sets.x.label = "Number of expressed TEs",
    mainbar.y.label = "Number of TEs in intersection",
    keep.order = TRUE
)

# -----------------------------------------------------------------------------------------------------------
# Correlation of overlapping TEs in brain
# -----------------------------------------------------------------------------------------------------------

brain_overlap <- intersect(deTEs$male_brain, deTEs$female_brain)

# Filter for the specific set of TEs and rearrange the table to have female on x and male on y axis

detes_brain_female <- deseq_tes %>% 
    filter(te_id %in% brain_overlap, sex == 'female', tissue == 'brain') %>% 
    dplyr::select(te_id, log2FoldChange) %>% 
    rename(log2FC_female = log2FoldChange)

detes_brain_male <- deseq_tes %>% 
    filter(te_id %in% brain_overlap, sex == 'male', tissue == 'brain') %>% 
    dplyr::select(te_id, log2FoldChange) %>% 
    rename(log2FC_male = log2FoldChange)

detes_brain_merged <- inner_join(detes_brain_female, detes_brain_male, by = 'te_id')

cor_brain <- cor.test(detes_brain_merged$log2FC_female, detes_brain_merged$log2FC_male, method = 'pearson')

p_brain <- ggplot(detes_brain_merged, aes(x = log2FC_female, y = log2FC_male)) +
    geom_point(color = 'gray23', alpha = 0.7, size = 2) +
    geom_smooth(method = 'lm', color = tissue.color['brain'], fill = tissue.color['brain'], alpha = 0.3) +
    theme_classic() +
    labs(title = paste0('Brain: R = ', round(cor_brain$estimate, 2), ', p-value = ', signif(cor_brain$p.value, 3)),
         x = 'log2FC Female',
         y = 'log2FC Male') +
    theme(text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16))


# Where is the specific set of TEs located in a genomic context?

deseq_tes %>% 
    filter(te_id %in% brain_overlap, sex == 'female', tissue == 'brain') %>% 
    ggplot(aes(x = position, fill = position)) +
    geom_bar()

deseq_tes %>% 
    filter(te_id %in% brain_overlap, sex == 'female', tissue == 'brain') %>% 
    ggplot(aes(x = order, fill = position)) +
    geom_bar()

knitr::kable(
deseq_tes %>% 
    filter(te_id %in% brain_overlap, sex == 'female', tissue == 'brain') %>% 
    dplyr::select(te_id, ensembl_gene_id, external_gene_name) %>% 
    splitTEID('te_id'),
format = "markdown")

# -----------------------------------------------------------------------------------------------------------
# UpSet-plot for DETEs of the public data set
# -----------------------------------------------------------------------------------------------------------

aging_tes::load_rna_seq_public_data_env()

deseq_te_public <- fread(paste0(rna_seq_deseq_dir, "deseq_results_te_instances_public.csv")) %>% 
    mutate(sex_tissue = paste0(sex, "_", tissue))



deTEs <- sapply(unique(deseq_te_public$sex_tissue), function(x){
    
    deseq_te_public %>% 
        filter(sex_tissue == x, padj <= 0.05) %>% 
        pull(te_id) %>% 
        unique()
})

upset(
    fromList(deTEs),
    nsets = length(unique(deseq_te_public$sex_tissue)),
    nintersects = 20,
    order.by = "freq",
    #sets.bar.color = tissue.color[c(2, 3, 4, 2, 3)],
    main.bar.color = "gray23",
    matrix.color = "gray23",
    point.size = 1.5,
    line.size = 1,
    text.scale = c(1.5, 1.2, 1.5, 1.2, 1.5, 1.2),
    #mb.ratio = c(0.6, 0.4),
    sets.x.label = "Number of expressed TEs",
    mainbar.y.label = "Number of TEs in intersection",
    keep.order = TRUE
)


