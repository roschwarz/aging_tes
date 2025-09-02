# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()

# =============================== Functions ====================================

read.region.feature.intersection <- function(file, region_type, feature_type){
    
    #browser()
    df <- read.csv(file)
    
    df <- df %>% mutate(feature_pos = ifelse(strand == '-', region.end - feature.start, feature.start - region.start),
                        relative_position = feature_pos*100/region.width,
                        relative_position_group = round(relative_position),
                        feature = feature_type,
                        region_type = region_type)
    
    return(df)
} 

makeRatioDataFrame <- function(df){
    
    sox <- df %>% 
        filter(feature == 'sox') %>% 
        group_by(relative_position_group) %>% 
        dplyr::count(name = 'sox_relative_pos_n')
    
    tss <- df %>% 
        filter(feature == 'tss') %>% 
        group_by(relative_position_group) %>% 
        dplyr::count(name = 'tss_relative_pos_n')
    
    ratio <- merge(tss, sox, by = 'relative_position_group', all.x = T) %>% 
        filter(!is.na(relative_position_group)) %>% 
        mutate(sox_relative_pos_n = ifelse(is.na(sox_relative_pos_n), 0, sox_relative_pos_n),
               ratio = sox_relative_pos_n/tss_relative_pos_n,
               feature = 'ratio') %>% 
        dplyr::rename( relative_position = relative_position_group )
    
    return(ratio)
    
    
}


# ==============================================================================
# ================================ Brain =======================================
# ==============================================================================

annotations <- "./results/cage_seq_gclipped/brain_downsampled_gclipped/annotations/"

indie_te_island_feature_intersection_files <- list(tss.body = paste0(annotations, 'brain_downsampled_independent_TE_regions.tss.intersection.bed'),
                                          tss.down = paste0(annotations, 'brain_downsampled_TE_region_down_stream.tss.intersection.bed'),
                                          tss.up =   paste0(annotations, 'brain_downsampled_TE_region_up_stream.tss.intersection.bed'),
                                          sox.body = paste0(annotations, 'brain_downsampled_independent_TE_regions.sox.intersection.bed' ),
                                          sox.down = paste0(annotations, 'brain_downsampled_TE_region_down_stream.sox.intersection.bed'),
                                          sox.up = paste0(annotations, 'brain_downsampled_TE_region_up_stream.sox.intersection.bed' ))

indie_te_island_feature_rel_pos <- do.call('rbind', sapply(names(indie_te_island_feature_intersection_files), simplify = F, function(x){
    
    feature_type = str_split(x, "[.]")[[1]][1]
    region_type = str_split(x, "[.]")[[1]][2]
    
    df <-  read.region.feature.intersection(indie_te_island_feature_intersection_files[[x]], 
                                            region_type = region_type, 
                                            feature_type = feature_type)
    
    if (region_type == 'down') {
        df <- df %>%
            mutate(relative_position = relative_position + 100,
                   relative_position_group = relative_position_group + 100)
    }else if (region_type == 'up') {
        
        df <- df %>%
            mutate(relative_position = relative_position - 100,
                   relative_position_group = relative_position_group - 100)
    }
    
    return(df)
    
    
}))

indie_te_island_feature_rel_pos$region_type <- 
    factor(indie_te_island_feature_rel_pos$region_type, levels = c('up', 'body', 'down'))


write.table(indie_te_island_feature_rel_pos,
            file = paste0(table_dir, 'brain_independent_TE_island_relative_sox_positions.csv'),
            row.names = F,
            sep = ',',
            quote = F)

# ==============================================================================
# Plot - Ratio of Sox and transcription start sites (TSS)
# ==============================================================================

y_max = 1100
y_text = 1000
text_size = 8
geom_text_size = 8/.pt #0.36 * text_size

my_theme <-  theme(legend.position = "none",
                   legend.title = element_blank(),
                   legend.background = element_blank(),
                   legend.direction = 'horizontal',
                   plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"),
                   panel.grid = element_blank(),
                   axis.text = element_text(family = 'arial',
                                            size = 8),
                   axis.title = element_text(family = 'arial',
                                             size = 8,
                                             face = 'bold'),
                   panel.background = element_blank())

#############
# body plot #
#############

p3 <- indie_te_island_feature_rel_pos %>% 
    filter(dplyr::between(relative_position, 0, 100)) %>% 
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, 
                   position = 'dodge', 
                   boundary = 0, 
                   closed = "left")

p3 <- p3 +
    geom_segment(data = indie_te_island_feature_rel_pos %>% 
                     filter(feature == 'sox', dplyr::between(relative_position, 0, 100)) %>% 
                     filter(!duplicated(as.factor(round(relative_position, 1)))),
                 aes(x = relative_position, y = -10,
                     xend = relative_position, yend = -40),
                 size = 0.5) 

ratio <- makeRatioDataFrame(indie_te_island_feature_rel_pos)

p3 <- p3 +
    geom_path(data = ratio %>% filter(dplyr::between(relative_position, 0, 100)), aes(y = ratio*3520), size = 1) +
    scale_y_continuous(sec.axis = sec_axis(~./3520, name = "# Sox/# TSS"),
                       expand = c(0.01, 0.01), limits = c(-50,y_max)) + 
    scale_x_continuous(expand = c(0.002,0.002), limits = c(0,100))



#Saving ratio axis to grob
g <- ggplotGrob(p3)
g_sec_axis <- gtable_filter(g, 'axis-r|ylab-r', trim = T)

p3 <- p3 + 
    theme(legend.position = "none", 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.title.x = element_blank(),
    ) 

text.p3 <- data.frame(x = 50, y = y_text, label = c('TE island'))

p3 <- p3 + 
    scale_color_manual(values = c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values = c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p3, 
              aes(x = x, y = y, label = label), 
              inherit.aes = F, 
              size = geom_text_size, 
              family = 'arial') +
    labs(y = 'count', x = 'relative position in TE island [%]') +
    my_theme + 
    theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"),
          legend.position = c(0.3, 0.9),
          legend.direction = 'vertical')

##################
# up-stream plot #
##################

text.p1 <- data.frame(x = -50, y = y_text, label = c('up-stream'))

p1 <- indie_te_island_feature_rel_pos %>% 
    filter(dplyr::between(-1*relative_position,0.5 , 100)) %>% 
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left")  + 
    theme(legend.position = "none", panel.grid.minor = element_blank(), axis.title.x = element_blank()) + 
    scale_x_continuous(expand = c(0,0), breaks = c(-25,-50,-75,-100)) + 
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(-50,y_max))

p1 <- p1 +  
    geom_segment(data = indie_te_island_feature_rel_pos %>% 
                     filter(dplyr::between(-1*relative_position, 0.5 , 100)) %>% 
                     filter(feature == 'sox') %>% 
                     filter(!duplicated(as.factor(round(relative_position, 1)))), 
                 aes(x = relative_position, y = -10, 
                     xend = relative_position, yend = -40), 
                 size = 0.09) 

p1 <- p1 + 
    scale_color_manual(values = c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values = c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p1, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count',
         x = 'relative position in TE island [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 10), "pt"))

#p1

####################
# down-stream plot #
####################

text.p2 <- data.frame(x = 150,
                      y = y_text,
                      label = c('down-stream'))


p2 <- indie_te_island_feature_rel_pos %>% 
    filter(dplyr::between(relative_position,100.5 , 200)) %>% 
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left") + 
    theme(legend.position = "none", 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.x = element_blank(), 
          panel.grid.minor = element_blank()) + 
    scale_x_continuous(expand = c(0,0), breaks = c(125,150,175,200)) + 
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(-50,y_max))


p2 <- p2 + 
    geom_segment(data = (indie_te_island_feature_rel_pos %>% 
                             filter(dplyr::between(relative_position,100.5 , 200)) %>% 
                             filter(feature == 'sox') %>% 
                             filter(!duplicated(as.factor(round(relative_position, 1))))) , 
                 aes(x = relative_position, y = -10, 
                     xend = relative_position, yend = -40), 
                 size = 0.09) 

p2  <- p2 +    
    scale_color_manual(values = c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values = c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p2, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count',
         x = 'relative position in TE-region [%]') +
    my_theme + 
    theme(plot.margin = unit(c(t = 0, r = 35, b = 0, l = 0), "pt"),
          axis.text.y.left = element_blank(),
          axis.title.y.left = element_blank(),
          axis.line.y.left = element_blank(),
          axis.ticks.y.left = element_blank())


#Add axis
p2 <- p2 + 
    annotation_custom(g_sec_axis, xmin = 211)

#p2

cowplot::plot_grid(p1, p3, p2, rel_widths = c(1.5, 3, 1.5), labels = c("","",""), align = "h", nrow = 1)

ggsave(
    filename = paste0(figure_dir, "panel_3_brain_indie_te_island_tss_sox_ratio.pdf"),
    device = cairo_pdf,
    plot = last_plot(),
    width = 20,
    height = 15,
    units = "cm",
    dpi = 300)



# ==============================================================================
# Plot - Counts of TE instances overlapping with Sox
# ==============================================================================

aging_tes::load_annotations()
load_te_ranges() 

te_island_5_prime_extended <- readGeneric('data/processed/annotation/mm10_TE_island_5prime_extended.bed',
                                          strand = 6,
                                          meta.cols = list(names = 4))


te_island_instances <- blackRcloud::intersectGranger(te_island_5_prime_extended, teRanges, tab = 'all')                 
                                                                                                                        
# TEs can intersect with multiple genes, so that they occur multiple times in                                           
# the teRange object. Therefore a unique is applied at the end of the pipeline                                          
te_island_instances <- te_island_instances %>%                                                                          
       dplyr::select(query.names, subject.te.id, subject.position) %>%                                                     
       filter(!duplicated(subject.te.id)) %>%                                                                              
       unique()                                                                                                            
                                                                                                                         
names(te_island_instances) <- c("te_island_id", "te_id", "te_position")      

indie_te_island_sox_overlap <- indie_te_island_feature_rel_pos %>% 
    filter(dplyr::between(relative_position, 0, 100)) %>% 
    filter(feature == 'sox', region_type == 'body') %>% 
    dplyr::pull(region.id)

sox_te_islands <- te_island_instances %>% 
    filter(te_island_id %in% indie_te_island_sox_overlap)

# Motifs of interest (Sox5 in the current case) are extracted with homer and intersected with TE instances.
# bedtools intersect -a ../../../../data/shared/mm10_transposome_sorted.bed -b Sox5.coordinates.bed > Sox.TEs.bed
sox_intersecting_tes <- read.csv('./results/cage_seq_gclipped//brain_downsampled_gclipped/homer/Sox.TEs.bed', sep = '\t', 
                                 col.names = c('chromosome',
                                               'start',
                                               'end',
                                               'te_id',
                                               'value',
                                               'strand'))

sox_intersecting_tes <- sox_intersecting_tes %>% 
    filter(te_id %in% sox_te_islands$te_id) %>% 
    filter(!duplicated(te_id)) %>% 
    dplyr::select(te_id)

sox_intersecting_tes <- splitTEID(sox_intersecting_tes, 'te_id') %>% 
    mutate(te_width = as.numeric(end) - as.numeric(start))

sox_intersecting_tes$color <- ifelse(dplyr::between(sox_intersecting_tes$te_width,  950, 1050), TRUE, FALSE)

pl_length_count <- ggplot(sox_intersecting_tes, aes(te_width, fill = color)) +
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

sox_intersecting_tes %>% filter(dplyr::between(te_width, 950, 1050)) %>% nrow()

pl_family_count <- ggplot(order.TEs(sox_intersecting_tes %>% filter(dplyr::between(te_width, 950, 1050)), 'family', decreasing = F), 
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

ggsave(filename = paste0(figure_dir, 'panel_3_length_of_sox_intersecting_tes_w_inset.pdf'),
       device = cairo_pdf,
       plot = pl_with_inset,
       width = 20,
       height = 9,
       units = 'cm',
       dpi = 300)


# ==============================================================================
# ================================ Blood =======================================
# ==============================================================================

annotations <- "./results/cage_seq_gclipped//blood_gclipped/annotations/"

blood_region_feature_intersection_files <- list(tss.body = paste0(annotations, 'blood_independent_TE_regions.tss.intersection.bed'),
                                          tss.down = paste0(annotations, 'blood_TE_region_down_stream.tss.intersection.bed'),
                                          tss.up = paste0(annotations, 'blood_TE_region_up_stream.tss.intersection.bed'),
                                          sox.body = paste0(annotations, 'blood_independent_TE_regions.sox.intersection.bed' ),
                                          sox.down = paste0(annotations, 'blood_TE_region_down_stream.sox.intersection.bed'),
                                          sox.up = paste0(annotations, 'blood_TE_region_up_stream.sox.intersection.bed' ))


blood_indie_te_island_feature_rel_pos <- do.call('rbind', sapply(names(blood_region_feature_intersection_files), simplify = F, function(x){
    
    feature_type = str_split(x, "[.]")[[1]][1]
    region_type = str_split(x, "[.]")[[1]][2]
    
    df <-  read.region.feature.intersection(blood_region_feature_intersection_files[[x]],
                                            region_type = region_type,
                                            feature_type = feature_type)
    
    if (region_type == 'down') {
        df <- df %>%
            mutate(relative_position = relative_position + 100,
                   relative_position_group = relative_position_group + 100)
    }else if (region_type == 'up'){
        
        df <- df %>%
            mutate(relative_position = relative_position - 100,
                   relative_position_group = relative_position_group - 100)
    }
    
    
    
    return(df)
    
    
}))

blood_indie_te_island_feature_rel_pos$region_type <-
    factor(blood_indie_te_island_feature_rel_pos$region_type, levels = c('up', 'body', 'down'))


write.table(blood_indie_te_island_feature_rel_pos,
            file = paste0(table_dir, 'blood_independent_TE_island_relative_sox_positions.csv'),
            row.names = F,
            sep = ',',
            quote = F)


# ==============================================================================
# Plot
# ==============================================================================

y_max = 150
y_text = 100
text_size = 12
geom_text_size = 0.36 * text_size

my_theme <-  theme(legend.position = "none",
                   legend.title = element_blank(),
                   legend.background = element_blank(),
                   legend.direction = 'horizontal',
                   plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"),
                   panel.grid = element_blank(),
                   axis.text = element_text(family = 'arial',
                                            size = 10),
                   axis.title = element_text(family = 'arial',
                                             size = 12,
                                             face = 'bold'),
                   panel.background = element_blank())

#############
# body plot #
#############

p3 <- blood_indie_te_island_feature_rel_pos %>%
    filter(dplyr::between(relative_position, 0, 100)) %>%
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left")

p3 <- p3 +
    geom_segment(data = blood_indie_te_island_feature_rel_pos %>%
                     filter(feature == 'sox', dplyr::between(relative_position, 0, 100)) %>%
                     filter(!duplicated(as.factor(round(relative_position, 1)))),
                 aes(x = relative_position, y = -1,
                     xend = relative_position, yend = -6),
                 size = 0.5)

ratio <- makeRatioDataFrame(blood_indie_te_island_feature_rel_pos)

p3 <- p3 +
    geom_path(data = ratio %>% filter(dplyr::between(relative_position, 0, 100)), aes(y = ratio*500), size = 1)+
    scale_y_continuous(sec.axis = sec_axis(~./500, name = "# Sox/# TSS"), expand = c(0.01, 0.01), limits=c(-7,y_max)) +
    scale_x_continuous(expand = c(0.002,0.002), limits=c(0,100))


p3

#Saving ratio axis to grob

g <- ggplotGrob(p3)
g_sec_axis <- gtable_filter(g, 'axis-r|ylab-r', trim=T)

p3 <- p3 +
    theme(legend.position = "none",
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())

text.p3 <- data.frame(x = 50, y = y_text, label = c('TE island'))

p3 <- p3 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p3, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count', x = 'relative position in TE island [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"))

##################
# up-stream plot #
##################

text.p1 <- data.frame(x = -50, y = y_text, label = c('up-stream'))

dummy_sox <- blood_indie_te_island_feature_rel_pos %>% filter(dplyr::between(relative_position,100.5 , 200), feature == 'sox') %>% 
    mutate(relative_position = -0.6,
           relative_position_group = -1,
           region_type = 'up')

p1 <- rbind(dummy_sox, blood_indie_te_island_feature_rel_pos) %>%
    filter(dplyr::between(-1*relative_position,0.5 , 100)) %>%
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth=2.5, position = 'dodge', boundary = 0, closed = "left")  +
    theme(legend.position = "none", panel.grid.minor = element_blank(), axis.title.x = element_blank()) +
    scale_x_continuous(expand = c(0,0), breaks=c(-25,-50,-75,-100)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits=c(-7,y_max))

p1 <- p1 +
    geom_segment(data = blood_indie_te_island_feature_rel_pos %>%
                     filter(dplyr::between(-1*relative_position, 0.5 , 100))%>%
                     filter(feature == 'sox') %>%
                     filter(!duplicated(as.factor(round(relative_position, 1)))),
                 aes(x = relative_position, y = -1,
                     xend = relative_position, yend = -6),
                 size = 0.09)

p1 <- p1 + scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p1, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count',
         x = 'relative position in TE-region [%]') +
    my_theme + theme(legend.position = c(0.3, 0.9))

p1

####################
# down-stream plot #
####################

text.p2 <- data.frame(x = 150,
                      y=y_text,
                      label = c('down-stream'))


p2 <- blood_indie_te_island_feature_rel_pos %>%
    filter(dplyr::between(relative_position,100.5 , 200)) %>%
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    # histogram plus vertical lines to separate different regions
    geom_histogram(binwidth=2.5, position = 'dodge', boundary = 0, closed = "left") +
    theme(legend.position = "none",
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0,0), breaks=c(125,150,175,200)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits=c(-7,y_max))


p2 <- p2 +
    geom_segment(data = (blood_indie_te_island_feature_rel_pos %>%
                             filter(dplyr::between(relative_position,100.5 , 200)) %>%
                             filter(feature == 'sox') %>%
                             filter(!duplicated(as.factor(round(relative_position, 1))))) ,
                 aes(x = relative_position, y = -1,
                     xend = relative_position, yend = -6),
                 size = 0.09)


p2  <- p2 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p2, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count',
         x = 'relative position in TE island [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 35, b = 0, l = 0), "pt"),
          axis.text.y.left = element_blank(),
          axis.title.y.left = element_blank(),
          axis.line.y.left = element_blank(),
          axis.ticks.y.left = element_blank())


#Add axis
p2 <- p2 +
    annotation_custom(g_sec_axis, xmin=211)
p2

cowplot::plot_grid(p1, p3, p2, rel_widths = c(1.8, 3, 1.8), labels = c("","",""), align = "h", nrow = 1)

ggsave(
    filename = paste0(figure_dir, "supplemental_panel_4_blood_indie_te_island_tss_sox_ratio.svg") ,
    plot = last_plot(),
    width = 30,
    height = 20,
    units = "cm",
    dpi = 300)

# ==============================================================================
# ================================ Skin ========================================
# ==============================================================================

skin_annotations <- "./results/cage_seq_gclipped/skinII_gclipped/annotations/"

skin_region_feature_intersection_files <- list(tss.body = paste0(skin_annotations, 'skinII_independent_TE_regions.tss.intersection.bed'),
                                          tss.down = paste0(skin_annotations, 'skin_TE_region_down_stream.tss.intersection.bed'),
                                          tss.up =   paste0(skin_annotations, 'skin_TE_region_up_stream.tss.intersection.bed'),
                                          sox.body = paste0(skin_annotations, 'skinII_independent_TE_regions.sox.intersection.bed' ),
                                          sox.down = paste0(skin_annotations, 'skin_TE_region_down_stream.sox.intersection.bed'),
                                          sox.up = paste0(skin_annotations, 'skin_TE_region_up_stream.sox.intersection.bed' ))
# 
skin_indie_te_island_feature_rel_pos <- do.call('rbind', sapply(names(skin_region_feature_intersection_files), simplify = F, function(x){
    
    feature_type = str_split(x, "[.]")[[1]][1]
    region_type = str_split(x, "[.]")[[1]][2]
    
    df <-  read.region.feature.intersection(skin_region_feature_intersection_files[[x]],
                                            region_type = region_type,
                                            feature_type = feature_type)
    
    if(region_type == 'down'){
        df <- df %>%
            mutate(relative_position = relative_position + 100,
                   relative_position_group = relative_position_group + 100)
    }else if(region_type == 'up'){
        
        df <- df %>%
            mutate(relative_position = relative_position - 100,
                   relative_position_group = relative_position_group - 100)
    }
    
    return(df)
}))

skin_indie_te_island_feature_rel_pos$region_type <-
    factor(skin_indie_te_island_feature_rel_pos$region_type, levels = c('up', 'body', 'down'))


write.table(skin_indie_te_island_feature_rel_pos,
            file = paste0(table_dir, 'skin_independent_TE_island_relative_sox_positions.csv'),
            row.names = F,
            sep = ',',
            quote = F)


# ==============================================================================
# Plot
# ==============================================================================

y_max = 250
y_text = 180
text_size = 12
geom_text_size = 0.36 * text_size

my_theme <-  theme(legend.position = "none",
                   legend.title =element_blank(),
                   legend.background = element_blank(),
                   legend.direction = 'horizontal',
                   plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"),
                   panel.grid = element_blank(),
                   axis.text = element_text(family = 'arial',
                                            size = 10),
                   axis.title = element_text(family = 'arial',
                                             size = 12,
                                             face = 'bold'),
                   panel.background = element_blank())

#############
# body plot #
#############

p3 <- skin_indie_te_island_feature_rel_pos %>%
    filter(dplyr::between(relative_position, 0, 100)) %>%
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left")

p3 <- p3 +
    geom_segment(data = skin_indie_te_island_feature_rel_pos %>%
                     filter(feature == 'sox', dplyr::between(relative_position, 0, 100)) %>%
                     filter(!duplicated(as.factor(round(relative_position, 1)))),
                 aes(x = relative_position, y = -1,
                     xend = relative_position, yend = -6),
                 size = 0.5)

ratio <- makeRatioDataFrame(skin_indie_te_island_feature_rel_pos)

p3 <- p3 +
    geom_path(data = ratio %>% filter(dplyr::between(relative_position, 0, 100)), aes(y=ratio*750), size = 1)+
    scale_y_continuous(sec.axis = sec_axis(~./750, name = "# Sox/# TSS"), expand = c(0.01, 0.01), limits=c(-7,y_max)) +
    scale_x_continuous(expand = c(0.002,0.002), limits=c(0,100))

p3


#Saving ratio axis to grob

g <- ggplotGrob(p3)
g_sec_axis <- gtable_filter(g, 'axis-r|ylab-r', trim=T)

p3 <- p3 +
    theme(legend.position = "none",
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())

text.p3 <- data.frame(x = 50, y = y_text, label = c('TE island'))

p3 <- p3 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p3, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count', x = 'relative position in TE island [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"))

##################
# up-stream plot #
##################

text.p1 <- data.frame(x = -50, y = y_text, label = c('up-stream'))

p1 <- skin_indie_te_island_feature_rel_pos %>%
    filter(dplyr::between(-1*relative_position,0.5 , 100)) %>%
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth=2.5, position = 'dodge', boundary = 0, closed = "left")  +
    theme(legend.position = "none", panel.grid.minor = element_blank(), axis.title.x = element_blank()) +
    scale_x_continuous(expand = c(0,0), breaks=c(-25,-50,-75,-100)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits=c(-7,y_max))

p1 <- p1 +
    geom_segment(data = skin_indie_te_island_feature_rel_pos %>%
                     filter(dplyr::between(-1*relative_position, 0.5 , 100))%>%
                     filter(feature == 'sox') %>%
                     filter(!duplicated(as.factor(round(relative_position, 1)))),
                 aes(x = relative_position, y = -1,
                     xend = relative_position, yend = -6),
                 size = 0.09)

p1 <- p1 + scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p1, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count',
         x = 'relative_position in TE-region [%]') +
    my_theme + theme(legend.position = c(0.3, 0.9))

p1

####################
# down-stream plot #
####################

text.p2 <- data.frame(x = 150,
                      y = y_text,
                      label = c('down-stream'))


dummy_sox <- skin_indie_te_island_feature_rel_pos %>% filter(dplyr::between(-1*relative_position,60 , 100), feature == 'sox') %>% 
    mutate(relative_position = 101,
           relative_position_group = 101,
           region_type = 'down')

p2 <- rbind(skin_indie_te_island_feature_rel_pos, dummy_sox) %>%
    filter(dplyr::between(relative_position,100.5 , 200)) %>%
    ggplot(aes(relative_position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left") +
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0,0), breaks = c(125,150,175,200)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(-7,y_max))


p2 <- p2 +
    geom_segment(data = (skin_indie_te_island_feature_rel_pos %>%
                             filter(dplyr::between(relative_position,100.5 , 200)) %>%
                             filter(feature == 'sox') %>%
                             filter(!duplicated(as.factor(round(relative_position, 1))))) ,
                 aes(x = relative_position, y = -1,
                     xend = relative_position, yend = -6),
                 size = 0.09)

p2  <- p2 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p2, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count',
         x = 'relative_position in TE-region [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 35, b = 0, l = 0), "pt"),
          axis.text.y.left = element_blank(),
          axis.title.y.left = element_blank(),
          axis.line.y.left = element_blank(),
          axis.ticks.y.left = element_blank())


#Add axis
p2 <- p2 +
    annotation_custom(g_sec_axis, xmin = 211)
p2

cowplot::plot_grid(p1, p3, p2, rel_widths = c(1.8, 3, 1.8), labels = c("","",""), align = "h", nrow = 1)

ggsave(
    filename = paste0(figure_dir, "supplemental_panel_4_skin_indie_te_island_tss_sox_ratio.svg") ,
    plot = last_plot(),
    width = 30,
    height = 20,
    units = "cm",
    dpi = 300)


# ==============================================================================
# Count of TE Families overlapping with sox motifs
# ==============================================================================

##### Figure C & D ########
#
# bedtools intersect -a <(bedtools intersect -a data/shared/mm10_transposome.bed \
#       -b data/shared/blood_independent_TE_regions.bed ) \
#   -b results/cage_RS/blood_gclipped/homer/RBP1.coordinates.bed > results/cage_RS/blood_gclipped/homer/blood.sox.tes.bed

skin_sox10_TE_intersection <- read.csv('./results/cage_seq_gclipped//skinII_gclipped/homer/skin.sox.tes.bed',
                                       head = F,
                                       sep = '\t')

skin_fam_count <- skin_sox10_TE_intersection %>% 
    filter(!duplicated(V4)) %>% 
    splitTEID('V4') %>%  
    filter(super_family != 'NA') %>% 
    count(family)


skin_fam_count$family <- factor(skin_fam_count$family, levels = skin_fam_count[order(skin_fam_count$n), 'family'])


skin_fam_count_pl <- ggplot(skin_fam_count, aes(family, n)) +
    geom_col(fill = '#264653') +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) + 
    labs(title = 'skin',
         y = "count of sox intersecting TEs",
         x = "TE family") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 12))

ggsave(
    filename = paste0(figure_dir, 'fam_count_sox_intersection_skin.svg'),
    plot = skin_fam_count_pl,
    width = 7,
    height = 5,
    units = "cm",
    dpi = 300)

#write.csv(skin_fam_count, 
#          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_4/S4C_fam_count_skin.csv')

blood_sox4_TE_intersection <- read.csv('./results/cage_seq_gclipped/blood_gclipped/homer/blood.sox.tes.bed',
                                       head = F,
                                       sep = '\t')

blood_fam_count <- blood_sox4_TE_intersection %>% 
    filter(!duplicated(V4)) %>% 
    splitTEID('V4') %>%  
    filter(super_family != 'NA') %>% 
    count(family)


blood_fam_count$family <- factor(blood_fam_count$family, levels = blood_fam_count[order(blood_fam_count$n), 'family'])


blood_fam_count_pl <- ggplot(blood_fam_count, aes(family, n)) +
    geom_col(fill = '#264653') +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0,0.2))) + 
    labs(title = 'blood',
         y = "count of sox intersecting TEs",
         x = "TE family") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 12))

ggsave(
    filename = paste0(figure_dir, 'fam_count_sox_intersection_blood.svg'),
    plot = blood_fam_count_pl,
    width = 7,
    height = 5,
    units = "cm",
    dpi = 300)

write.csv(blood_fam_count, 
          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_4/S4D_fam_count_blood.csv')

##### Figure  E #####

sox_TE_intersection <- read.csv('./results/cage_seq_gclipped/brain_downsampled_gclipped/homer/Sox.TEs.bed',
                                head = F,
                                sep = '\t')


super_fam_count <- sox_TE_intersection %>% 
    filter(!duplicated(V4)) %>% 
    splitTEID('V4') %>% 
    filter(super_family != 'NA') %>% 
    count(super_family) %>% 
    slice_max(order_by = n, n = 20) %>% 
    mutate(super_family = case_when(super_family == "Alu" ~ "B1", .default = super_family))

super_fam_count$super_family <- factor(super_fam_count$super_family, levels = super_fam_count[order(super_fam_count$n), 'super_family'])

ggplot(super_fam_count, aes(super_family, n)) +
    geom_col() +
    scale_y_log10(expand = c(0,0)) +
    coord_flip() +
    labs(y = "count of sox intersecting TEs",
         x = "TE super family") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))


ggsave(
    filename = paste0(figure_dir, 'sox_motif_TE_intersection_brain.svg'),
    plot = last_plot(),
    width = 10,
    height = 8,
    units = "cm",
    dpi = 300)

write.csv(super_fam_count,
          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_4/S4E_Sox.motif.TE.intersection.csv')

