# Load Environment
if( !exists("ENVIRONMENT_LOADED") ){
    
    source('./01_load_environment.R')
    
} else if( !ENVIRONMENT_LOADED ){
    
    source('./01_load_environment.R')
    
}

go2gene <- load.GO.Term.Annotation() # load goTerms + Ensembl ids

intronic.tes <- transGrange(teRanges) %>% 
    filter(position == 'intronic') %>% 
    group_by(ensembl_gene_id) %>% 
    dplyr::count(name = 'n.intronicTEs')

gene.intronic.tes <- 
    merge(transGrange(geneRanges), 
          intronic.tes, by = 'ensembl_gene_id', all.x = T)

gene.intronic.tes$n.intronicTEs <- ifelse(is.na(gene.intronic.tes$n.intronicTEs), 
                                          0 , gene.intronic.tes$n.intronicTEs)

target.wo.TEs <- gene.intronic.tes %>% 
    filter(n.intronicTEs == 0) %>% 
    pull(ensembl_gene_id)

back <- gene.intronic.tes %>% 
    pull(ensembl_gene_id)


go.enrichments.res.wo.TEs <- do.call('rbind', sapply(names(go2gene), 
                                                     simplify = F, 
                                                     function(ontology){
                                                         
                                                         onto.go.set <- go2gene[[ontology]]
                                                         
                                                         onto.results <- applyFisherGo(target.wo.TEs, 
                                                                                       back, 
                                                                                       onto.go.set, 
                                                                                       ontology, 
                                                                                       minGenes = 10, 
                                                                                       maxGenes = 500)
                                                         
                                                         onto.results <- filter(onto.results, p.adjust <= FDR)
                                                         
                                                         goterms <- as.data.frame(GOTERM)[2:3] %>% 
                                                             filter(!duplicated(go_id)) %>% 
                                                             rename(Description = Term)
                                                         
                                                         
                                                         onto.results <- merge(onto.results, as.data.frame(goterms), by='go_id', all.x = T)
                                                         
                                                         return(onto.results)
                                                         
                                                     }))

data_clean_wo_TE_enriched <- go.enrichments.res.wo.TEs %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_wo_TE_enriched$Description <- factor(data_clean_wo_TE_enriched$Description,
                                                level = data_clean_wo_TE_enriched[order(data_clean_wo_TE_enriched$p.adjust, decreasing = T), 'Description'])

#GO.plot(data_clean_wo_TE_enriched)

data_clean_wo_TE_depleted <- go.enrichments.res.wo.TEs %>% 
    filter(enriched == 'depleted') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_wo_TE_depleted$Description <- factor(data_clean_wo_TE_depleted$Description,
                                                level = data_clean_wo_TE_depleted[order(data_clean_wo_TE_depleted$p.adjust, decreasing = T), 'Description'])

#GO.plot(data_clean_wo_TE_depleted)


data_clean_wo_TE <- rbind(data_clean_wo_TE_depleted, rev(data_clean_wo_TE_enriched))

data_clean_wo_TE$log.p.adjust <- ifelse(data_clean_wo_TE$enriched == 'depleted', 
                                        log10(data_clean_wo_TE$p.adjust), 
                                        -log10(data_clean_wo_TE$p.adjust))

ggplot(data_clean_wo_TE, aes(log.p.adjust, Description, fill = ONTOLOGY)) +
    geom_col() +
    theme_rob(10, base_family = 'arial') +
    facet_grid(ONTOLOGY~., scales = 'free_y', switch = "y") +
    labs(x = expression(log[10](FDR))) +
    geom_text(aes(label = paste0(target.count, '/', go.count), hjust = -0.2), color = 'white', family = 'arial', size = 7*0.36) + # size from pt into mm 8 (pt) * 0.36
    geom_text(aes(label = paste0(target.count, '/', go.count), hjust = +1.2), color = 'white', family = 'arial', size = 7*0.36) + # size from pt into mm 8 (pt) * 0.36
    scale_fill_manual(values = c('#98afba','#647a85','#264653')) +
    scale_x_break(c(-5, -60)) +
    theme(axis.title.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.spacing.y = unit(0.25, "lines"),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.5),
          strip.text.y = element_text(angle = 360),
          legend.position = 'None',
          strip.placement = "outside",
          strip.background = element_rect(color = 'white'),
          plot.margin = unit(c(t=0, r = 0, b = 0, l = 0),'pt'))


# The following would be for resis, but I copied it from the old RESIS directory.
# 
# ggsave(
#     filename = paste0('../figures/28_GO_gene_wo_TEs.svg'),
#     plot = last_plot(),
#     width = 20,
#     height = 20,
#     units = "cm",
#     dpi = 300)
# 
# 
# write.table(data_clean_wo_TE,
#             file = '../../submission/Resis/Supplemental_figure5/Supplemental_figure5.csv', 
#             col.names = F, 
#             row.names = F,
#             quote = F,
#             sep = ',')