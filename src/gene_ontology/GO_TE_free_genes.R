# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

load_analysis_env()
load_annotations()
load_plotting_env()

load_te_ranges()
load_gene_ranges()
go2gene <- load_GO_Term_Annotation()

intronic_tes <- transGrange(teRanges) %>% 
    filter(position == 'intronic') %>% 
    group_by(ensembl_gene_id) %>% 
    dplyr::count(name = 'n.intronicTEs')

gene_intronic_tes <- 
    merge(transGrange(geneRanges), 
          intronic_tes, by = 'ensembl_gene_id', all.x = T)

gene_intronic_tes$n.intronicTEs <- ifelse(is.na(gene_intronic_tes$n.intronicTEs), 
                                          0 , gene_intronic_tes$n.intronicTEs)

target_wo_TEs <- gene_intronic_tes %>% 
    filter(n.intronicTEs == 0) %>% 
    pull(ensembl_gene_id)

back <- gene_intronic_tes %>% 
    pull(ensembl_gene_id)


go_enrichments_res_wo_TEs <- do.call('rbind', sapply(names(go2gene), 
                                                     simplify = F, 
                                                     function(ontology){
                                                         
                                                         onto_go_set <- go2gene[[ontology]]
                                                         
                                                         onto_results <- applyFisherGo(target_wo_TEs, 
                                                                                       back, 
                                                                                       onto_go_set, 
                                                                                       ontology, 
                                                                                       minGenes = 10, 
                                                                                       maxGenes = 500)
                                                         
                                                         onto_results <- filter(onto_results, p.adjust <= FDR)
                                                         
                                                         goterms <- as.data.frame(GOTERM)[2:3] %>% 
                                                             filter(!duplicated(go_id)) %>% 
                                                             rename(Description = Term)
                                                         
                                                         
                                                         onto_results <- merge(onto_results, as.data.frame(goterms), by='go_id', all.x = T)
                                                         
                                                         return(onto_results)
                                                         
                                                     }))

data_clean_wo_TE_enriched <- go_enrichments_res_wo_TEs %>% 
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

data_clean_wo_TE_depleted <- go_enrichments_res_wo_TEs %>% 
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

capture.output(sessionInfo(),
               file = "./src/gene_ontology/GO_TE_free_genes_sessionInfo.txt"
)
