# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}


# Collect TE instances of TE-regions to get information about the composition of the TE
# regions later in the analysis.
#TE.region.instances <- intersectGranger(TE.region, teRanges, tab = 'all')
TE.region.instances <- intersectGranger(TEregionRanges_5prime, teRanges, tab = 'all')

# TEs can intersect with multiple genes, so that they occur multiple times in
# the teRange objects. Therefore a unique is applied at the end of the pipeline
TE.region.instances <- TE.region.instances %>% 
    dplyr::select(query.names, subject.te.id, subject.position) %>% 
    filter(!duplicated(subject.te.id)) %>% 
    unique()

names(TE.region.instances) <- c("te_region_id", "te_id", "te_position")



auto.TE.cluster.files <- list(brain = "data/shared/brain_downsampled_independent_TE_regions.bed",
                              blood = "data/shared/blood_independent_TE_regions.bed",
                              skin = "data/shared/skinII_independent_TE_regions.bed"
)

auto_TE_regions <- sapply(tissues, simplify = F, function(x){
        df <- read.csv(auto.TE.cluster.files[[x]], 
                       sep = '\t', 
                       header = F) %>% 
            dplyr::pull(V4)
    })







#############################
# Composition of TE regions #
#############################
# Take the TE region instance table and filter for all TE regions that are 
# autonomously expressed. Determine the position of the instance within the TE
# region and categorize the TE regions into single (consists of one instance),
# double (consists of two instances) and multiple (consists of more than 2 instance)

auto_TE_region_composition <- sapply(names(auto_TE_regions), simplify = F, function(x){
    
    auto.TE.region.instances <- TE.region.instances %>% 
        filter(te_region_id %in% auto_TE_regions[[x]])
    
    
    names(auto.TE.region.instances) <- c('te_region_id', 'te_id', 'te_genomic_pos')
    
    auto.TE.region.instances <- splitTEID(auto.TE.region.instances, 'te_id')
    
    auto.TE.region.instances <- auto.TE.region.instances %>% 
        group_by(te_region_id) %>% 
        mutate(member = n(),
               position = ifelse(strand == '-', rev(1:n()), 1:n()),
               region.type = ifelse(member == 1, 'single', 
                                    ifelse(member == 2, 'double', 'multiple')))
    
    
    return(auto.TE.region.instances)
})



auto_region_types <- sapply(c('brain', 'skin', 'blood'), simplify = F, function(x){
    
    df <- auto_TE_region_composition[[x]]
    
    df <- df %>% filter(!duplicated(te_region_id))
    
    df$region.type <- factor(df$region.type, levels = c('single', 'double', 'multiple'))
    
    print(paste(x, ":", nrow(df)))
    
    pl <- ggplot(df, aes(region.type)) +
        geom_bar() +
        labs(title = x,
             x = 'type of TE island') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
              axis.text.y = element_text(size = 6),
              axis.title = element_text(size = 8),
              plot.title = element_text(size = 10)) 
    
    return(pl)
    
})

auto_region_cat <- gridExtra::grid.arrange(grobs = auto_region_types, nrow = 1)


auto_single_region <- sapply(names(auto_TE_region_composition), simplify = F, function(x){
    
    df <- auto_TE_region_composition[[x]]
    
    df <-  df %>% 
        filter(member == 1) %>%
        mutate(super_family = case_when(super_family == "Alu" ~ "B1", .default = super_family)) %>% 
        order.TEs('super_family', decreasing = F) %>% 
        group_by(order) %>% 
        count(super_family) %>%
        ungroup() %>% 
        slice_max(order_by = n, n = 15) 
    
    pl <- ggplot(df, aes(super_family, n)) +
        geom_col() +
        coord_flip() +
        labs(title = x, 
             y = 'count') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
              axis.title = element_text(size = 8),
              plot.title = element_text(size = 10)) 
    
    
    return(pl)
})

single_te_island_compo_pl <- gridExtra::grid.arrange(grobs = auto_single_region, nrow = 1)


auto_double_region <- sapply(names(auto_TE_region_composition), simplify = F, function(x){
    
    df <- auto_TE_region_composition[[x]]
    
    df <-  df %>% 
        filter(member == 2) %>%
        mutate(super_family = case_when(super_family == "Alu" ~ "B1", .default = super_family)) %>% 
        order.TEs('super_family', decreasing = F) %>% 
        group_by(order) %>% 
        count(super_family) %>%
        ungroup() %>% 
        slice_max(order_by = n, n = 15) 
    
    pl <- ggplot(df, aes(super_family, n)) +
        geom_col() +
        coord_flip() +
        labs(title = x, 
             y = 'count') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
              axis.text.y = element_text(size = 6),
              axis.title = element_text(size = 8),
              plot.title = element_text(size = 10)) 
    
    return(pl)
})

double_te_island_compo_pl <- gridExtra::grid.arrange(grobs = auto_double_region, nrow = 1)


auto_multi_region <- sapply(names(auto_TE_region_composition), simplify = F, function(x){
    
    df <- auto_TE_region_composition[[x]]
    
    df <-  df %>% 
        filter(member > 2) %>%
        mutate(super_family = case_when(super_family == "Alu" ~ "B1", .default = super_family)) %>% 
        order.TEs('super_family', decreasing = F) %>% 
        group_by(order) %>% 
        count(super_family) %>%
        ungroup() %>% 
        slice_max(order_by = n, n = 15) 
    
    pl <- ggplot(df, aes(super_family, n)) +
        geom_col() +
        coord_flip() +
        labs(title = x, 
             y = 'count') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
              axis.text.y = element_text(size = 6),
              axis.title = element_text(size = 8),
              plot.title = element_text(size = 10)) 
    
    
    return(pl)
})

multi_te_island_compo_pl <- gridExtra::grid.arrange(grobs = auto_multi_region, nrow = 1)


cairo_pdf(filename = paste0(figure_dir, 's03_characterizing_te_island.pdf'), width = 7.9, height = 4.3)
ggarrange(auto_region_cat, 
          single_te_island_compo_pl,
          double_te_island_compo_pl,
          multi_te_island_compo_pl,
          font.label = list(size = 12),
          ncol = 2,
          nrow = 2,
          labels = c("A", "B", "C", "D"),
          align = 'hv'
             )
dev.off()


#########
# RESIS #
#########

ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_3/S3A_TE_island_types.pdf',
    plot = auto_region_cat,
    device = cairo_pdf,
    width = 10,
    height = 5,
    units = "cm",
    dpi = 300)


auto_region_types_csv <- do.call('rbind', 
                             sapply(c('brain', 'skin', 'blood'), simplify = F, function(x){
    
    df <- auto_TE_region_composition[[x]]
    
    df <- df %>% filter(!duplicated(te_region_id))
    
    df$region.type <- factor(df$region.type, levels = c('single', 'double', 'multiple'))
    
    print(paste(x, ":", nrow(df)))
    
    
    df$tissue <- x
    
    
    return(df)
    
    })

)

write.csv(auto_region_types_csv, 
          file = './manuscripts/nature_aging/resis/figures/supplemental_figure_3/S3_TE_island_types.csv')

ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_3/S3B_TE_island_single.pdf',
    plot = single_te_island_compo_pl ,
    device = cairo_pdf,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 300)

ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_3/S3C_TE_island_double.pdf',
    plot = double_te_island_compo_pl ,
    device = cairo_pdf,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 300)

ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_3/S3D_TE_island_multi.pdf',
    plot = multi_te_island_compo_pl ,
    device = cairo_pdf,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 300)


 
###### Write TE island annotation of the category multi #####

TEregionAnnotation <- data.frame(TEregionRanges_5prime)

sapply(names(auto_TE_region_composition), simplify = F, function(x){
    te_islands_list <- auto_TE_region_composition[[x]] %>% 
        filter(member > 2) %>% 
        pull(te_region_id) %>% 
        unique()
    
    te_islands <- TEregionAnnotation %>% 
        filter(names %in% te_islands_list) %>% 
        dplyr::select(seqnames, start, end, names, width, strand)
    
    writeBedFile(te_islands,
                 file.name = paste0("s03_", x, "_multi_TE_island.bed"), directory = data_dir)
})
