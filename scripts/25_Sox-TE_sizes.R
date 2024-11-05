# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}


relative.positions <- read.csv(paste0(table_dir, 
                                      '24_independent_TE_regions_relative_positions.csv'))

sox.auto.TE.region.ids <- relative.positions %>% 
    filter(between(relative.position, 0, 100)) %>% 
    filter(feature == 'sox', region.type == 'body') %>% 
    dplyr::pull(region.id)

te.region.instances <- read.csv(te.region.instances.file)

sox.TE.regions <- te.region.instances %>% 
    filter(te_region_id %in% sox.auto.TE.region.ids)


# Motifs of interest (Sox5 in the current case) are extracted with homer and intersected with TE instances.
# bedtools intersect -a ../../../../data/shared/mm10_transposome_sorted.bed -b Sox5.coordinates.bed > Sox.TEs.bed
sox.intersecting.TEs <- read.csv('./results/cage_RS/brain_downsampled_gclipped/homer/Sox.TEs.bed', sep = '\t', 
                                 col.names = c('chromosome',
                                               'start',
                                               'end',
                                               'TE.ID',
                                               'value',
                                               'strand'))

sox.intersecting.TEs <- sox.intersecting.TEs %>% 
    filter(TE.ID %in% sox.TE.regions$te_id) %>% 
    filter(!duplicated(TE.ID)) %>% 
    dplyr::select(TE.ID)

sox.intersecting.TEs <- splitTEID(sox.intersecting.TEs, 'TE.ID') %>% 
    mutate(te.width = as.numeric(end) - as.numeric(start))

sox.intersecting.TEs$color <- ifelse(between(sox.intersecting.TEs$te.width,  950, 1050), TRUE, FALSE)


pl_length_count <- ggplot(sox.intersecting.TEs, aes(te.width, fill = color)) +
    geom_histogram(binwidth = 10) +
    scale_fill_manual(values = c('TRUE' = '#264653', 'FALSE' = '#a6a6a6')) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,250)) +
    labs(x = "length of sox intersecting TEs [bp]",
         y = "count of TE instances") +
    #theme_rob(12, base_family = 'arial') +
    theme(panel.grid = element_blank(),
          legend.position = 'None',
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8))


sox.intersecting.TEs %>% filter(between(te.width, 950, 1050)) %>% nrow()

pl_family_count <- ggplot(order.TEs(sox.intersecting.TEs %>% filter(between(te.width, 950, 1050)), 'family', decreasing = F), 
       aes(family)) +
    geom_bar(fill = '#264653') +
    scale_y_continuous(expand = c(0.0, 0.0), limits = c(0, 51)) +
    coord_flip() +
    labs(x = 'TE subfamily',
         y = 'count of TE instances') +
    theme(panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          )

pl_with_inset <- ggdraw() +
    draw_plot(pl_length_count) +
    draw_plot(pl_family_count, x = 0.25, y = .13, width = .5, height = 0.86)

ggsave(filename = paste0(figure_dir, '25_te_length_w_inset.pdf'),
       device = cairo_pdf,
       plot = pl_with_inset,
       width = 20,
       height = 9,
       units = 'cm',
       dpi = 300)


#########
# RESIS #
#########

write.csv(sox.intersecting.TEs,
          file = './manuscripts/nature_aging/resis/figures/figure3/3F_soxTEs_length.csv',
          row.names = FALSE)


write.csv(sox.intersecting.TEs %>% filter(between(te.width, 950, 1050)),
          file = './manuscripts/nature_aging/resis/figures/figure3/3F_soxTEs_length_inset.csv',
          row.names = FALSE)

ggsave(filename = './manuscripts/nature_aging/resis/figures/figure3/3F_soxTEs_length.pdf',
       device = cairo_pdf,
       plot = pl_with_inset,
       width = 20,
       height = 9,
       units = 'cm',
       dpi = 300)


######### Under Construnction ###########
#
# Basic idea extract the coordinates of the specific set of TE and store them
# as bed file. Use this .bed-file to run Homer

sox.intersecting.TEs.bed <- sox.intersecting.TEs %>% 
    filter(between(te.width, 950, 1050)) %>% 
    dplyr::select(chromosome, start, end, TE.ID, te.width, strand)

write.table(sox.intersecting.TEs.bed,
            file = paste0(data_dir, '25_sox.intersecting.TEs.bed'),
            sep = '\t',
            col.names = FALSE,
            quote = FALSE,
            row.names = FALSE)

