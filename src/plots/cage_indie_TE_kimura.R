
# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_cage_seq_env()
aging_tes::load_annotations()
aging_tes::load_plotting_env()

cageRanges <- load_cage_peak_annotation()
teRanges <- load_te_ranges()
te_annotation <- load_te_annotation()

# ------------------------------------------------------------------------------
# Identify transposable elements with their own TSS
# ------------------------------------------------------------------------------

te_cages <- sapply(names(cageRanges), simplify = F, function(x) {
    
    df <- intersectGranger(cageRanges[[x]], teRanges, 'all')
    
    names(df) <- str_replace_all(names(df), c('query' = 'cage',
                                              'subject' = 'te')) #makes te.te.id, is a bit annoying
    
    return(df)
    
})

cage_TEs <- do.call('rbind', sapply(names(te_cages), simplify = F, function(x){
    
    te_cages[[x]] %>% 
        dplyr::select('te.te.id', 'cage.names') %>% 
        splitTEID("te.te.id") %>% 
        #filter(order %in% orders.of.interest) %>% 
        mutate(tissue = x) %>% 
        dplyr::rename(te_id = te.te.id,
                      peak_id = cage.names)
    
}))


cage_TEs %>% 
    filter(super_family == "L1", grepl("3end", family)) %>% 
    ggplot() +
    geom_segment(aes(x = Kimura, xend = Kimura+0.1, y = family, yend = family), size = 10) +
    facet_grid(.~tissue)


l1md <- cage_TEs %>% filter(tissue == "brain", super_family == "L1", grepl("L1Md", family)) %>% 
    mutate(type = case_when(grepl("3end", family) ~ '3end',
                            grepl("5end", family) ~ '5end', .default = 'body')
           )

l1md %>% 
    ggplot() +
    #geom_segment(aes(x = Kimura, xend = Kimura+0.1, y = 0, yend = 1), size = 0.5) +
    stat_interval(aes(Kimura, family), .width = c(.1, .25, .5, .75, 1), height = 5, show.legend = F) +
    stat_halfeye(aes(Kimura, family), .width = 0, fill = "#F28705", size = 0.5, point_color = "white", )

te_annotation %>% filter(grepl("L1Md", family)) %>% 
    ggplot(aes(as.numeric(Kimura))) +
    geom_density()
