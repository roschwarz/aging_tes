# --------------------------------- Notes--- -----------------------------------
# 
# This script calculates the proportion of TE classes of (differentially) 
# expressed TEs 
# 
#

# Load Environment
if( !exists("ENVIRONMENT_LOADED") ){
    source('./01_load_environment.R')
} else if( !ENVIRONMENT_LOADED ){
    
    source('./01_load_environment.R')
    
}

#deseqTEs <- read.csv('./results/tables/02_deseq_results_te.csv')

##### R is making issues. You can check the panel1.Rmd file in ../manuscript_results
## as template if anything is working again. Header: Characterization of expressed TEs

if (!exists("deseq.te.merged")) {
    
    deseq.te.merged <- read.csv(paste0(table_dir, '02_deseq_results_te_instances.csv'))
    
}

expressed.tes <- deseq.te.merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'expressed')

background <- transGrange(teRanges) %>% 
    splitTEID('te.id') %>% 
    filter(order %in% orders.of.interest) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           tissue = 'background',
           group = 'expressed')

dete.up <- deseq.te.merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= 0.05, log2FoldChange > 0) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'up')

dete.down <- deseq.te.merged %>% 
    filter(!is.na(padj)) %>% 
    splitTEID('te_id') %>% 
    filter(order %in% orders.of.interest, padj <= 0.05, log2FoldChange < 0) %>% 
    group_by(tissue) %>% 
    dplyr::count(order) %>% 
    mutate(percent = n/sum(n),
           group = 'down')

df.classes <- rbind(background, expressed.tes, dete.up, dete.down)

df.classes$group = factor(df.classes$group, levels = c('expressed','down','up'))
df.classes$tissue = factor(df.classes$tissue, levels = c('background', 'brain', 'skin', 'blood'))

pl.te.percent <- ggplot(df.classes, aes(tissue, percent, fill = order)) +
    geom_col(width = 0.95) +
    scale_fill_manual(values = order.color,
                      name = 'TE class:') +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0,0))) +
    labs(y = 'Percent of\n(differentially) expressed TEs') +
    facet_grid(.~group, scales = 'free_x',space='free') +
    theme_rob(base_family = 'arial',
              base_size = 10) +
    theme(axis.title.x = element_blank(),
          legend.position = 'bottom',
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1
                                     ),
          legend.background = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-15,-10, 0,-10),
          legend.key.size = unit(0.5,"line")
          )

pl.te.percent <- color_strips(pl.te.percent,
                              bg_cols = c('white', direction.color[['down']], direction.color[['up']]), 
                              text_cols = c( "#000000", "#000000","#000000"))


grid.draw(pl.te.percent)

ggsave(pl.te.percent,
       filename = paste0(figure_dir, '5_expr_te_class_4x3_300.pdf'),
       device = cairo_pdf,
       width = 9.2,
       height = 7.5,
       units = "cm",
       dpi = 300
)


library(ggalluvial)

df.class.updated <- df.classes %>% 
    mutate(x = paste0(tissue, '_', group))

order.of.x = c('background_expressed',
               'brain_expressed',
               'skin_expressed',
               'blood_expressed',
               'brain_down',
               'skin_down',
               'blood_down',
               'brain_up',
               'skin_up',
               'blood_up')

df.class.updated$x <- factor(df.class.updated$x, levels = order.of.x)


pl.te.percent.updated <- ggplot(df.class.updated, aes(x, percent, fill = order)) +
    geom_flow(aes(alluvium = order), 
               alpha = .3, 
               color = 'white',
               curve_type = 'linear', 
               width = .3) +
    geom_col(width = .7, color = "white") +
    #facet_grid(.~group, scales = 'free_x',space='free') +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0,0))) +
    scale_x_discrete(expand = c(0,0),labels= c('background_expressed' = 'background',
                              'brain_expressed' = 'brain',
                              'skin_expressed' = 'skin',
                              'blood_expressed' = 'blood',
                              'brain_down' = 'brain',
                              'skin_down' = 'skin',
                              'blood_down' = 'blood',
                              'brain_up' = 'brain',
                              'skin_up' = 'skin',
                              'blood_up' = 'blood')) +
    scale_fill_manual(values = order.color,
                      name = 'TE class:') +
    labs(y = 'Percent of\n(differentially) expressed TEs') +
   # geom_text(aes(y = 1, x = 2, label = c('expressed'))) +
    #theme(plot.margin = margin(20, 2, 2, 2, unit = "pt"))
    geom_vline(xintercept = c(1.5, 4.5, 7.5), color = "black", linewidth = .5) +
    #geom_hline(yintercept = 1, color = "black", linewidth = .5) +
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
          #plot.margin = margin(20, 2, 2, 2, unit = "pt")
    )
    
ggsave(pl.te.percent.updated,
       filename = paste0(figure_dir, '5_expr_te_class_updated_4x3_300.pdf'),
       device = cairo_pdf,
       width = 9.2,
       height = 7.5,
       units = "cm",
       dpi = 300
)


#########
# RESIS #
#########

ggsave(pl.te.percent.updated,
       filename = './manuscripts/nature_aging/resis/Figure1/1C_TE_proporiton.pdf',
       device = cairo_pdf,
       width = 9.2,
       height = 7.5,
       units = "cm",
       dpi = 300
)


write.csv(df.class.updated, file = './manuscripts/nature_aging/resis/Figure1/1C_TE_proporiton.csv')
