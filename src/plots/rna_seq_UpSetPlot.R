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
