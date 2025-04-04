# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}

#system('bash ./prepFencePlot_V2.sh', wait = TRUE)

# =============================== Functions ====================================

read.region.feature.intersection <- function(file, region.type, feature.type){
    
    df <- read.csv(file)
    
    df <- df %>% mutate(feature.pos = ifelse(strand == '-', region.end - feature.start, feature.start - region.start),
                        relative.position = feature.pos*100/region.width,
                        relative.position.group = round(relative.position),
                        feature = feature.type,
                        region.type = region.type)
} 

makeRatioDataFrame <- function(df){
    
    sox <- df %>% 
        filter(feature == 'sox') %>% 
        group_by(relative.position.group) %>% 
        dplyr::count(name = 'sox.relative.pos.n')
    
    tss <- df %>% 
        filter(feature == 'tss') %>% 
        group_by(relative.position.group) %>% 
        dplyr::count(name = 'tss.relative.pos.n')
    
    ratio <- merge(tss, sox, by = 'relative.position.group', all.x = T) %>% 
        filter(!is.na(relative.position.group)) %>% 
        mutate(sox.relative.pos.n = ifelse(is.na(sox.relative.pos.n), 0, sox.relative.pos.n),
               ratio = sox.relative.pos.n/tss.relative.pos.n,
               feature = 'ratio') %>% 
        dplyr::rename( relative.position = relative.position.group )
    
    return(ratio)
    
    
}

# ================================ Brain =======================================

annotations <- "./results/cage_RS/brain_downsampled_gclipped/annotations/"

region.feature.intersection.files <- list(tss.body = paste0(annotations, 'brain_downsampled_independent_TE_regions.tss.intersection.bed'),
                                          tss.down = paste0(annotations, 'brain_downsampled_TE_region_down_stream.tss.intersection.bed'),
                                          tss.up =   paste0(annotations, 'brain_downsampled_TE_region_up_stream.tss.intersection.bed'),
                                          sox.body = paste0(annotations, 'brain_downsampled_independent_TE_regions.sox.intersection.bed' ),
                                          sox.down = paste0(annotations, 'brain_downsampled_TE_region_down_stream.sox.intersection.bed'),
                                          sox.up = paste0(annotations, 'brain_downsampled_TE_region_up_stream.sox.intersection.bed' ))

auto.feature.rel.pos <- do.call('rbind', sapply(names(region.feature.intersection.files), simplify = F, function(x){
    
    feature.type = str_split(x, "[.]")[[1]][1]
    region.type = str_split(x, "[.]")[[1]][2]
    
    df <-  read.region.feature.intersection(region.feature.intersection.files[[x]], 
                                            region.type = region.type, 
                                            feature.type = feature.type)
    
    if (region.type == 'down') {
        df <- df %>%
            mutate(relative.position = relative.position + 100,
                   relative.position.group = relative.position.group + 100)
    }else if (region.type == 'up') {
        
        df <- df %>%
            mutate(relative.position = relative.position - 100,
                   relative.position.group = relative.position.group - 100)
    }
    
    return(df)
    
    
}))

auto.feature.rel.pos$region.type <- 
    factor(auto.feature.rel.pos$region.type, levels = c('up', 'body', 'down'))


write.table(auto.feature.rel.pos,
            file = paste0(table_dir, '24_independent_TE_regions_relative_positions.csv'),
            row.names = F,
            sep = ',',
            quote = F)


write.table(auto.feature.rel.pos,
            file = './manuscripts/nature_aging/resis/figures/figure3/3E_brain_TSS_Sox_in_independent_TE_regions.csv',
            row.names = F,
            sep = ',',
            quote = F)

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

p3 <- auto.feature.rel.pos %>% 
    filter(between(relative.position, 0, 100)) %>% 
    ggplot(aes(relative.position, fill = feature, color = feature)) +
        geom_histogram(binwidth = 2.5, 
                       position = 'dodge', 
                       boundary = 0, 
                       closed = "left")

p3 <- p3 +
    geom_segment(data = auto.feature.rel.pos %>% 
                     filter(feature == 'sox', between(relative.position, 0, 100)) %>% 
                     filter(!duplicated(as.factor(round(relative.position, 1)))),
                 aes(x = relative.position, y = -10,
                     xend = relative.position, yend = -40),
                 size = 0.5) 

ratio <- makeRatioDataFrame(auto.feature.rel.pos)

p3 <- p3 +
    geom_path(data = ratio %>% filter(between(relative.position, 0, 100)), aes(y = ratio*3520), size = 1) +
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

p1 <- auto.feature.rel.pos %>% 
    filter(between(-1*relative.position,0.5 , 100)) %>% 
    ggplot(aes(relative.position, fill = feature, color = feature)) +
    geom_histogram(binwidth = 2.5, position = 'dodge', boundary = 0, closed = "left")  + 
    theme(legend.position = "none", panel.grid.minor = element_blank(), axis.title.x = element_blank()) + 
    scale_x_continuous(expand = c(0,0), breaks = c(-25,-50,-75,-100)) + 
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(-50,y_max))

p1 <- p1 +  
    geom_segment(data = auto.feature.rel.pos %>% 
                     filter(between(-1*relative.position, 0.5 , 100)) %>% 
                     filter(feature == 'sox') %>% 
                     filter(!duplicated(as.factor(round(relative.position, 1)))), 
                 aes(x = relative.position, y = -10, 
                     xend = relative.position, yend = -40), 
                 size = 0.09) 

p1 <- p1 + 
    scale_color_manual(values = c(sox = '#264653', tss = 'gray80', ratio = '#9E2A2B')) +
    scale_fill_manual(values = c(sox = '#264653',  tss = 'gray80', ratio = '#9E2A2B')) +
    geom_text(data = text.p1, aes(x = x, y = y, label = label), inherit.aes = F, size = geom_text_size, family = 'arial') +
    labs(y = 'count',
         x = 'relative position in TE-region [%]') +
    my_theme +
    theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 10), "pt"))

#p1

####################
# down-stream plot #
####################

text.p2 <- data.frame(x = 150,
                      y = y_text,
                      label = c('down-stream'))


p2 <- auto.feature.rel.pos %>% 
    filter(between(relative.position,100.5 , 200)) %>% 
    ggplot(aes(relative.position, fill = feature, color = feature)) +
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
    geom_segment(data = (auto.feature.rel.pos %>% 
                             filter(between(relative.position,100.5 , 200)) %>% 
                             filter(feature == 'sox') %>% 
                             filter(!duplicated(as.factor(round(relative.position, 1))))) , 
                 aes(x = relative.position, y = -10, 
                     xend = relative.position, yend = -40), 
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

ggsave(filename = paste0(figure_dir, "24_brain_TSS_Sox_in_independent_TE_regions.pdf"),
       last_plot(),
       width = )

ggsave(
    filename = paste0(figure_dir, "24_brain_TSS_Sox_in_independent_TE_regions_10x7.5.pdf"),
    device = cairo_pdf,
    plot = last_plot(),
    width = 20,
    height = 15,
    units = "cm",
    dpi = 300)



