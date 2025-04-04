# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}


read.region.feature.intersection <- function(file, region.type, feature.type){
    
    df <- read.csv(file)
    
    df <- df %>% mutate(feature.pos = ifelse(strand == '-', region.end - feature.start, feature.start - region.start),
                        relative.position = feature.pos*100/region.width,
                        relative.position.group = round(relative.position),
                        feature = feature.type,
                        region.type = region.type)
} 

#
# Fence plots
# 
# You find the code for the bed files in motifAnalysis.sh
#
# ================================= Blood ======================================
# 

annotations <- "./results/cage_RS/blood_gclipped/annotations/"

region.feature.intersection.files <- list(tss.body = paste0(annotations, 'blood_independent_TE_regions.tss.intersection.bed'),
                                          tss.down = paste0(annotations, 'blood_TE_region_down_stream.tss.intersection.bed'),
                                          tss.up = paste0(annotations, 'blood_TE_region_up_stream.tss.intersection.bed'),
                                          sox.body = paste0(annotations, 'blood_independent_TE_regions.sox.intersection.bed' ),
                                          sox.down = paste0(annotations, 'blood_TE_region_down_stream.sox.intersection.bed'),
                                          sox.up = paste0(annotations, 'blood_TE_region_up_stream.sox.intersection.bed' ))


auto.feature.rel.pos <- do.call('rbind', sapply(names(region.feature.intersection.files), simplify = F, function(x){
    
    feature.type = str_split(x, "[.]")[[1]][1]
    region.type = str_split(x, "[.]")[[1]][2]
    
    df <-  read.region.feature.intersection(region.feature.intersection.files[[x]],
                                            region.type = region.type,
                                            feature.type = feature.type)
    
    if(region.type == 'down'){
        df <- df %>%
            mutate(relative.position = relative.position + 100,
                   relative.position.group = relative.position.group + 100)
    }else if(region.type == 'up'){
        
        df <- df %>%
            mutate(relative.position = relative.position - 100,
                   relative.position.group = relative.position.group - 100)
    }
    
    
    
    return(df)
    
    
}))

auto.feature.rel.pos$region.type <-
    factor(auto.feature.rel.pos$region.type, levels = c('up', 'body', 'down'))


write.table(auto.feature.rel.pos,
            file = paste0(table_dir, 's04_blood_independent_TE_regions_relative_positions.csv'),
            row.names = F,
            sep = ',',
            quote = F)

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

p3 <- auto.feature.rel.pos %>%
    filter(between(relative.position, 0, 100)) %>%
    ggplot(aes(relative.position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left")

p3 <- p3 +
    geom_segment(data = auto.feature.rel.pos %>%
                     filter(feature == 'sox', between(relative.position, 0, 100)) %>%
                     filter(!duplicated(as.factor(round(relative.position, 1)))),
                 aes(x = relative.position, y = -1,
                     xend = relative.position, yend = -6),
                 size = 0.5)

ratio <- makeRatioDataFrame(auto.feature.rel.pos)

p3 <- p3 +
    geom_path(data = ratio %>% filter(between(relative.position, 0, 100)), aes(y=ratio*500), size = 1)+
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
          axis.text.y = element_blank(),
          axis.title.x=element_blank())

text.p3 <- data.frame(x = 50, y = y_text, label = c('TE region'))

p3 <- p3 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p3, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count', x = 'relative position in TE-region [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"))

##################
# up-stream plot #
##################

text.p1 <- data.frame(x = -50, y = y_text, label = c('up-stream'))

dummy_sox <- auto.feature.rel.pos %>% filter(between(relative.position,100.5 , 200), feature == 'sox') %>% 
    mutate(relative.position = -0.6,
           relative.position.group = -1,
           region.type = 'up')

p1 <- rbind(dummy_sox, auto.feature.rel.pos) %>%
    filter(between(-1*relative.position,0.5 , 100)) %>%
    ggplot(aes(relative.position, fill = feature, color = feature)) +
    geom_histogram(binwidth=2.5, position = 'dodge', boundary = 0, closed = "left")  +
    theme(legend.position = "none", panel.grid.minor = element_blank(), axis.title.x = element_blank()) +
    scale_x_continuous(expand = c(0,0), breaks=c(-25,-50,-75,-100)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits=c(-7,y_max))

p1 <- p1 +
    geom_segment(data = auto.feature.rel.pos %>%
                     filter(between(-1*relative.position, 0.5 , 100))%>%
                     filter(feature == 'sox') %>%
                     filter(!duplicated(as.factor(round(relative.position, 1)))),
                 aes(x = relative.position, y = -1,
                     xend = relative.position, yend = -6),
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


p2 <- auto.feature.rel.pos %>%
    filter(between(relative.position,100.5 , 200)) %>%
    ggplot(aes(relative.position, fill = feature, color = feature)) +
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
    geom_segment(data = (auto.feature.rel.pos %>%
                             filter(between(relative.position,100.5 , 200)) %>%
                             filter(feature == 'sox') %>%
                             filter(!duplicated(as.factor(round(relative.position, 1))))) ,
                 aes(x = relative.position, y = -1,
                     xend = relative.position, yend = -6),
                 size = 0.09)


p2  <- p2 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
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
    annotation_custom(g_sec_axis, xmin=211)
p2

cowplot::plot_grid(p1, p3, p2, rel_widths = c(1.8, 3, 1.8), labels = c("","",""), align = "h", nrow = 1)

ggsave(
    filename = paste0(figure_dir, "s04_blood_TSS_Sox_in_independent_TE_regions.svg") ,
    plot = last_plot(),
    width = 30,
    height = 20,
    units = "cm",
    dpi = 300)
# 
# 
# # ================================= skin ======================================
# 
annotations <- "./results/cage_RS/skinII_gclipped/annotations/"

region.feature.intersection.files <- list(tss.body = paste0(annotations, 'skinII_independent_TE_regions.tss.intersection.bed'),
                                          tss.down = paste0(annotations, 'skin_TE_region_down_stream.tss.intersection.bed'),
                                          tss.up =   paste0(annotations, 'skin_TE_region_up_stream.tss.intersection.bed'),
                                          sox.body = paste0(annotations, 'skinII_independent_TE_regions.sox.intersection.bed' ),
                                          sox.down = paste0(annotations, 'skin_TE_region_down_stream.sox.intersection.bed'),
                                          sox.up = paste0(annotations, 'skin_TE_region_up_stream.sox.intersection.bed' ))
# 
auto.feature.rel.pos <- do.call('rbind', sapply(names(region.feature.intersection.files), simplify = F, function(x){
    
    feature.type = str_split(x, "[.]")[[1]][1]
    region.type = str_split(x, "[.]")[[1]][2]
    
    df <-  read.region.feature.intersection(region.feature.intersection.files[[x]],
                                            region.type = region.type,
                                            feature.type = feature.type)
    
    if(region.type == 'down'){
        df <- df %>%
            mutate(relative.position = relative.position + 100,
                   relative.position.group = relative.position.group + 100)
    }else if(region.type == 'up'){
        
        df <- df %>%
            mutate(relative.position = relative.position - 100,
                   relative.position.group = relative.position.group - 100)
    }
    
    
    
    return(df)
    
    
}))

auto.feature.rel.pos$region.type <-
    factor(auto.feature.rel.pos$region.type, levels = c('up', 'body', 'down'))


write.table(auto.feature.rel.pos,
            file = paste0(table_dir, 's04_skin_independent_TE_regions_relative_positions.csv'),
            row.names = F,
            sep = ',',
            quote = F)

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

p3 <- auto.feature.rel.pos %>%
    filter(between(relative.position, 0, 100)) %>%
    ggplot(aes(relative.position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left")

p3 <- p3 +
    geom_segment(data = auto.feature.rel.pos %>%
                     filter(feature == 'sox', between(relative.position, 0, 100)) %>%
                     filter(!duplicated(as.factor(round(relative.position, 1)))),
                 aes(x = relative.position, y = -1,
                     xend = relative.position, yend = -6),
                 size = 0.5)

ratio <- makeRatioDataFrame(auto.feature.rel.pos)

p3 <- p3 +
    geom_path(data = ratio %>% filter(between(relative.position, 0, 100)), aes(y=ratio*750), size = 1)+
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
          axis.text.y = element_blank(),
          axis.title.x=element_blank())

text.p3 <- data.frame(x = 50, y = y_text, label = c('TE region'))

p3 <- p3 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p3, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count', x = 'relative position in TE-region [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "pt"))

##################
# up-stream plot #
##################

text.p1 <- data.frame(x = -50, y = y_text, label = c('up-stream'))

p1 <- auto.feature.rel.pos %>%
    filter(between(-1*relative.position,0.5 , 100)) %>%
    ggplot(aes(relative.position, fill = feature, color = feature)) +
    geom_histogram(binwidth=2.5, position = 'dodge', boundary = 0, closed = "left")  +
    theme(legend.position = "none", panel.grid.minor = element_blank(), axis.title.x = element_blank()) +
    scale_x_continuous(expand = c(0,0), breaks=c(-25,-50,-75,-100)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits=c(-7,y_max))

p1 <- p1 +
    geom_segment(data = auto.feature.rel.pos %>%
                     filter(between(-1*relative.position, 0.5 , 100))%>%
                     filter(feature == 'sox') %>%
                     filter(!duplicated(as.factor(round(relative.position, 1)))),
                 aes(x = relative.position, y = -1,
                     xend = relative.position, yend = -6),
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


dummy_sox <- auto.feature.rel.pos %>% filter(between(-1*relative.position,60 , 100), feature == 'sox') %>% 
    mutate(relative.position = 101,
           relative.position.group = 101,
           region.type = 'down')

p2 <- rbind(auto.feature.rel.pos, dummy_sox) %>%
    filter(between(relative.position,100.5 , 200)) %>%
    ggplot(aes(relative.position, fill = feature, color = feature)) +
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
    geom_segment(data = (auto.feature.rel.pos %>%
                             filter(between(relative.position,100.5 , 200)) %>%
                             filter(feature == 'sox') %>%
                             filter(!duplicated(as.factor(round(relative.position, 1))))) ,
                 aes(x = relative.position, y = -1,
                     xend = relative.position, yend = -6),
                 size = 0.09)

p2  <- p2 +
    scale_color_manual(values =c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values =c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
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
p2

cowplot::plot_grid(p1, p3, p2, rel_widths = c(1.8, 3, 1.8), labels = c("","",""), align = "h", nrow = 1)

ggsave(
    filename = paste0(figure_dir, "s04_skin_TSS_Sox_in_independent_TE_regions.svg") ,
    plot = last_plot(),
    width = 30,
    height = 20,
    units = "cm",
    dpi = 300)


##### Figure C & D ########
#
# bedtools intersect -a <(bedtools intersect -a data/shared/mm10_transposome.bed \
#       -b data/shared/blood_independent_TE_regions.bed ) \
#   -b results/cage_RS/blood_gclipped/homer/RBP1.coordinates.bed > results/cage_RS/blood_gclipped/homer/blood.sox.tes.bed

skin.sox10.TE.intersection <- read.csv('./results/cage_RS/skinII_gclipped/homer/skin.sox.tes.bed',
                                head = F,
                                sep = '\t')

skin_fam_count <- skin.sox10.TE.intersection %>% 
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
    filename = paste0(figure_dir, 's04_fam_count_skin.svg'),
    plot = skin_fam_count_pl,
    width = 7,
    height = 5,
    units = "cm",
    dpi = 300)

write.csv(skin_fam_count, 
          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_4/S4C_fam_count_skin.csv')
    
blood.sox4.TE.intersection <- read.csv('./results/cage_RS/blood_gclipped/homer/blood.sox.tes.bed',
                                head = F,
                                sep = '\t')

blood_fam_count <- blood.sox4.TE.intersection %>% 
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
    filename = paste0(figure_dir, 's04_fam_count_blood.svg'),
    plot = blood_fam_count_pl,
    width = 7,
    height = 5,
    units = "cm",
    dpi = 300)

write.csv(blood_fam_count, 
          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_4/S4D_fam_count_blood.csv')

##### Figure  E #####

sox.TE.intersection <- read.csv('./results/cage_RS/brain_downsampled_gclipped/homer/Sox.TEs.bed',
                                head = F,
                                sep = '\t')


super_fam_count <- sox.TE.intersection %>% 
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
    filename = paste0(figure_dir, 's04_Sox.motif.TE.intersection.svg'),
    plot = last_plot(),
    width = 10,
    height = 8,
    units = "cm",
    dpi = 300)

write.csv(super_fam_count,
          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_4/S4E_Sox.motif.TE.intersection.csv')
