# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_rna_seq_env()
aging_tes::load_annotations()

load_te_ranges()


# ------------------------------------------------------------------------------
# Step 1: Calculate the proportion of TE classes for all TEs (background)
# ------------------------------------------------------------------------------

deseq_te_merged_male <- fread(paste0(table_dir, deseq_results_te_csv_male)) %>% 
    mutate(sex = 'male',
           sex_tissue = paste0(sex, "_",tissue)) %>% 
    filter(!is.na(padj))

deseq_te_merged_female <- fread(paste0(table_dir, deseq_results_te_csv_female)) %>% 
    mutate(sex = 'female',
           sex_tissue = paste0(sex, "_", tissue)) %>% 
    filter(!is.na(padj))

deseq_tes <- rbind(deseq_te_merged_male, deseq_te_merged_female)

expressed_tes <- deseq_tes %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest) %>% 
    group_by(sex_tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'expressed')

background <- transGrange(teRanges) %>% 
    splitTEID('te.id') %>% 
    filter(order %in% orders.of.interest) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           sex_tissue = 'background',
           group = 'expressed')

# ------------------------------------------------------------------------------
# Step 2: Calculate the proportion of TE classes for expressed and differentially expressed TEs
# ------------------------------------------------------------------------------

dete_up <- deseq_tes %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= FDR, log2FoldChange > 0) %>% 
    group_by(sex_tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'up')

dete_down <- deseq_tes %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= FDR, log2FoldChange < 0) %>% 
    group_by(sex_tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'down')


# ------------------------------------------------------------------------------
# Step 3: Merge the data frames 
# ------------------------------------------------------------------------------

df_classes <- rbind(background, expressed_tes, dete_up, dete_down)

df_classes$group = factor(df_classes$group, levels = c('expressed','down','up'))
df_classes$sex_tissue = factor(df_classes$sex_tissue, levels = c('background', 'male_brain', 'male_skin', 'male_blood', "female_brain", 'female_skin'))


# ------------------------------------------------------------------------------
# Step 4: Create the plot
# ------------------------------------------------------------------------------

df_class_updated <- df_classes %>% 
    mutate(x = paste0(sex_tissue, '_', group))

order_of_x = c('background_expressed',
               'male_brain_expressed',
               'male_skin_expressed',
               'male_blood_expressed',
               'male_brain_down',
               'male_skin_down',
               'male_blood_down',
               'male_brain_up',
               'male_skin_up',
               'male_blood_up',
               'female_brain_expressed',
               'female_skin_expressed',
               'female_brain_down',
               'female_skin_down',
               'female_brain_up',
               'female_skin_up')

df_class_updated$x <- factor(df_class_updated$x, levels = order_of_x)


pl_te_percent_updated <- ggplot(df_class_updated, aes(x, percent, fill = order)) +
    geom_flow(aes(alluvium = order), 
              alpha = .3, 
              color = 'white',
              curve_type = 'linear', 
              width = .3) +
    geom_col(width = .7, color = "white") +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0,0))) +
    scale_x_discrete(expand = c(0,0),labels = c('background_expressed' = 'background',
                                               'male_brain_expressed' = 'brain',
                                               'male_skin_expressed' = 'skin',
                                               'male_blood_expressed' = 'blood',
                                               'male_brain_down' = 'brain',
                                               'male_skin_down' = 'skin',
                                               'male_blood_down' = 'blood',
                                               'male_brain_up' = 'brain',
                                               'male_skin_up' = 'skin',
                                               'male_blood_up' = 'blood',
                                               'female_brain_expressed' = 'brain',
                                               'female_skin_expressed' = 'skin',
                                               'female_brain_down' = 'brain',
                                               'female_skin_down' = 'skin',
                                               'female_brain_up' = 'brain',
                                               'female_skin_up' = 'skin'
                                               
                                               )) +
    scale_fill_manual(values = order.color,
                      name = 'TE class:') +
    labs(y = 'Percent of\n(differentially) expressed TEs') +
    geom_vline(xintercept = c(1.5, 4.5, 7.5, 10.5,12.5,14.5), color = "black", linewidth = .5) +
    theme_rob(base_family = 'arial',
              base_size = 10) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          legend.position = 'bottom',
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1
          ),
          legend.background = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-15,-0, 0,-10),
          legend.key.size = unit(0.5,"line"),
    )

# Save figure with metadata to index
meta <- list(name = 'percentage_expr_te_class',
             description = 'Proportion of TE classes among expressed or differentially expressed TE instances compared to the background. Add bars for female at 7.10.25.',
             tags = c('expression', 'rna-seq'),
             parameters = list(FDR = 0.05, tissues = c('brain', 'skin', 'blood')),
             script = 'rna_seq_te_class_composition.R'
)

fig_index(plot = pl_te_percent_updated,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 15,
          height = 7.5,
          dpi = 300,
          format = 'pdf')


