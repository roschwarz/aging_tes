if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}


GO.plot <- function(df){
    
    x.max <- max(-log10(df$p.adjust)) + 1
    
    p1 <- ggplot(df, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
        geom_col() +
        geom_segment(aes(x = (-log10(p.adjust)), 
                         xend = x.max, 
                         y = Description, 
                         yend = Description), 
                     color = 'gray50',
                     linetype  = 2) +
        theme_rob(10, base_family = 'arial') +
        facet_grid(ONTOLOGY~., scales = 'free_y', switch = "y") +
        labs(x = expression(log[10](FDR))) +
        geom_text(aes(label = Count, hjust = 1.5), color = 'white', family = 'arial', size = 7*0.36) + # size from pt into mm 8 (pt) * 0.36
        scale_fill_manual(values = c('#98afba','#647a85','#264653')) +
        scale_x_continuous(expand = expansion(mult = c(0, 0)),
                           #limits = c(0, 50)) +
                           limits = c(0, round(x.max))) +
        theme(axis.title.y = element_blank(),
              panel.grid.major = element_blank(), 
              panel.spacing.y = unit(0.25, "lines"),
              panel.border = element_blank(),
              axis.line = element_line(size = 0.5),
              strip.text.y = element_text(angle = 360),
              legend.position = 'None',
              strip.placement = "outside",
              strip.background = element_rect(color = 'white'),
              plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0),'pt'))
    
    return(p1)
}


brain <- readGeneric('data/shared/brain_can_TE_island.bed', strand = 6, meta.cols = list(names = 4))
skin <- readGeneric('data/shared/skin_can_TE_island.bed', strand = 6, meta.cols = list(names = 4))
blood <- readGeneric('data/shared/blood_can_TE_island.bed', strand = 6, meta.cols = list(names = 4))


# -------------------------- Background genes ----------------------------------

deseq.gene <- read.csv(paste0(table_dir, '02_deseq_results_gene.csv'))

expressed.genes <- sapply(tissues, simplify = F, function(x){
    
    deseq.gene %>% filter(tissue == x) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
    
})

target.genes <- blackRcloud::fullIntersectionResult(brain, teRanges) %>% data.frame() %>% 
    dplyr::pull(subject.ensembl_gene_id) %>% unique()

go.res.brain <- goOver(target.genes, expressed.genes$brain)

goBarFacetUpdate(go.res.brain@result)


data_clean_brain <- go.res.brain@result %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_brain$Description <- factor(data_clean_brain$Description,
                                       levels = data_clean_brain[order(data_clean_brain$p.adjust,
                                                                       decreasing = T),
                                                                 'Description'])

pl_brain <- GO.plot(data_clean_brain)

pl_brain


#Skin

target.genes <- blackRcloud::fullIntersectionResult(skin, teRanges) %>% data.frame() %>% 
    dplyr::pull(subject.ensembl_gene_id) %>% unique()

go.res.skin <- goOver(target.genes, expressed.genes$skin)

goBarFacetUpdate(go.res.skin@result)


data_clean_skin <- go.res.skin@result %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_skin$Description <- factor(data_clean_skin$Description,
                                       levels = data_clean_skin[order(data_clean_skin$p.adjust,
                                                                       decreasing = T),
                                                                 'Description'])

pl_skin <- GO.plot(data_clean_skin)

pl_skin

#blood

target.genes <- blackRcloud::fullIntersectionResult(blood, teRanges) %>% data.frame() %>% 
    dplyr::pull(subject.ensembl_gene_id) %>% unique()

go.res.blood <- goOver(target.genes, expressed.genes$blood)

goBarFacetUpdate(go.res.blood@result)


data_clean_blood <- go.res.blood@result %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_blood$Description <- factor(data_clean_blood$Description,
                                       levels = data_clean_blood[order(data_clean_blood$p.adjust,
                                                                       decreasing = T),
                                                                 'Description'])

pl_blood <- GO.plot(data_clean_blood)

pl_blood
##### Explore


target.genes.brain <- blackRcloud::fullIntersectionResult(brain, teRanges) %>% data.frame() %>% 
    dplyr::pull(subject.external_gene_name) %>% unique()

target.genes.skin <- blackRcloud::fullIntersectionResult(skin, teRanges) %>% data.frame() %>% 
    dplyr::pull(subject.external_gene_name) %>% unique()

target.genes.blood <- blackRcloud::fullIntersectionResult(blood, teRanges) %>% data.frame() %>% 
    dplyr::pull(subject.external_gene_name) %>% unique()


uni.target.symb <- intersect(target.genes.blood, intersect(target.genes.brain, target.genes.skin))

target.genes.brain <- blackRcloud::fullIntersectionResult(brain, teRanges) %>% data.frame() %>%
    filter(!is.na(subject.ensembl_gene_id)) %>% 
    dplyr::pull(subject.ensembl_gene_id) %>% unique()

target.genes.skin <- blackRcloud::fullIntersectionResult(skin, teRanges) %>% data.frame() %>% 
    filter(!is.na(subject.ensembl_gene_id)) %>% 
    dplyr::pull(subject.ensembl_gene_id) %>% unique()

target.genes.blood <- blackRcloud::fullIntersectionResult(blood, teRanges) %>% data.frame() %>% 
    filter(!is.na(subject.ensembl_gene_id)) %>% 
    dplyr::pull(subject.ensembl_gene_id) %>% unique()


bg.genes <- intersect(expressed.genes$brain, intersect(expressed.genes$skin, expressed.genes$skin))
uni.target <- intersect(target.genes.blood, intersect(target.genes.brain, target.genes.skin))

go.res.brain <- goOver(uni.target, bg.genes)

data_clean_brain <- go.res.brain@result %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_brain$Description <- factor(data_clean_brain$Description,
                                       levels = data_clean_brain[order(data_clean_brain$p.adjust,
                                                                       decreasing = T),
                                                                 'Description'])

pl_brain <- GO.plot(data_clean_brain)

pl_brain


# aging genes

aging_genes <- do.call('rbind', sapply(list.files('data/shared/', pattern = 'aging'), simplify = FALSE, function(x){
    
    df <- read.csv(file = paste0('data/shared/', x))
    
    row.names(df) <- NULL
    
    return(df)
    
}))
    
intersect(uni.target.symb, aging_genes$Symbol)
