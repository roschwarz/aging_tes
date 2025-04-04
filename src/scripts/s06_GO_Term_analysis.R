
# Load Environment
if( !exists("ENVIRONMENT_LOADED") ){
    
    source('./01_load_environment.R')
    
} else if( !ENVIRONMENT_LOADED ){
    
    source('./01_load_environment.R')
    
}


# ======================== Functions ===========================================
GO.plot.deg <- function(df, tissue){
    
    x.max <- round(max(-log10(df$p.adjust))) + 1
    #x.max <- round(x.max, -1) + 1
    
    p1 <- ggplot(df, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
        # geom_segment(aes(x=0, 
        #                  xend = x.max, 
        #                  y=Description, 
        #                  yend=Description), 
        #              color = 'gray50',
        #              linetype  = 2) +
        geom_col() +
        #theme_rob(8) +
        facet_grid(ONTOLOGY~., scales = 'free_y', switch = "y") +
        labs(x = expression(log[10](FDR)),
             title = tissue) +
        geom_text(aes(label = paste0(target.count, '/', go.count), 
                      hjust = 1.2), 
                  color = 'white', 
                  size = 6/.pt) +
        scale_fill_manual(values = c('#98afba','#647a85','#264653')) +
        scale_x_continuous(expand = expansion(mult = c(0, 0)),
                           limits = c(0, round(x.max))) +
        theme(axis.title.y = element_blank(),
              panel.grid.major = element_blank(), 
              panel.spacing.y = unit(0.25, "lines"),
              panel.border = element_blank(),
              axis.line = element_line(size = 0.5),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 8),
              strip.text.y = element_text(angle = 360),
              legend.position = 'None',
              strip.placement = "outside",
              strip.background = element_rect(color = 'white'),
              strip.text = element_text(size = 8),
              #plot.margin = unit(c(t=0, r = 0, b = 0, l =0),'pt')
              )
    
    return(p1)
}

# ==================== Genes with independent TE ===============================

go2gene <- load.GO.Term.Annotation() # load goTerms + Ensembl ids

deseq.gene <- read.csv(paste0(table_dir, '02_deseq_results_gene.csv'))

expressed_genes <- sapply(tissues, simplify = F, function(x){
    
    deseq.gene %>% filter(tissue == x) %>% 
        dplyr::pull(ensembl_gene_id) %>% 
        unique()
    
})

# ====================   GO enrichment for DEGs ================================

# ========================= Brain ==============================================

deg_targets <- sapply(tissues, simplify = F, function(x){
    
    deseq.gene %>% filter(tissue == x, padj <= FDR) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

deg_target_brain <- deg_targets$brain
background_brain <- expressed_genes$brain

brain_deg_go.enrichments.res <- do.call('rbind', sapply(names(go2gene), simplify = F, function(ontology){
    
    onto.go.set <- go2gene[[ontology]]
    
    onto.results <- applyFisherGo(deg_target_brain, 
                                  background_brain, 
                                  onto.go.set, 
                                  ontology, 
                                  minGenes = 10, 
                                  maxGenes = 500)
    
    onto.results <- filter(onto.results, p.adjust <= FDR)
    
    goterms <- as.data.frame(GOTERM)[2:3] %>% 
        filter(!duplicated(go_id)) %>% 
        rename(Description = Term)
    
    
    onto.results <- merge(onto.results, as.data.frame(goterms), by = 'go_id', all.x = T)
    
    return(onto.results)
    
})
)

go_deg_clean_brain <-  brain_deg_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, 
                     n = 10, 
                     with_ties = FALSE) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

go_deg_clean_brain$Description <- factor(go_deg_clean_brain$Description,
                                      level = go_deg_clean_brain[order(go_deg_clean_brain$p.adjust, decreasing = T), 'Description'])

go.deg.pl.brain <- GO.plot.deg(go_deg_clean_brain, 'brain')


#########
# RESIS #
#########

write.csv(go_deg_clean_brain, './manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6A_brain_GO_Term_diff_genes.csv')


ggsave(
    filename ='./manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6A_brain_GO_Term_diff_genes.pdf',
    plot = go.deg.pl.brain,
    device = cairo_pdf,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 300)
# ========================= skin ==============================================

deg_targets <- sapply(tissues, simplify = F, function(x){
    
    deseq.gene %>% filter(tissue == x, padj <= FDR) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

deg_target_skin <- deg_targets$skin
background_skin <- expressed_genes$skin

skin_deg_go.enrichments.res <- do.call('rbind', sapply(names(go2gene), simplify = F, function(ontology){
    
    onto.go.set <- go2gene[[ontology]]
    
    onto.results <- applyFisherGo(deg_target_skin, 
                                  background_skin, 
                                  onto.go.set, 
                                  ontology, 
                                  minGenes = 10, 
                                  maxGenes = 500)
    
    onto.results <- filter(onto.results, p.adjust <= FDR)
    
    goterms <- as.data.frame(GOTERM)[2:3] %>% 
        filter(!duplicated(go_id)) %>% 
        rename(Description = Term)
    
    
    onto.results <- merge(onto.results, as.data.frame(goterms), by = 'go_id', all.x = T)
    
    return(onto.results)
    
}))

go_deg_clean_skin <-  skin_deg_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

go_deg_clean_skin$Description <- factor(go_deg_clean_skin$Description,
                                         level = go_deg_clean_skin[order(go_deg_clean_skin$p.adjust, decreasing = T), 'Description'])

go.deg.pl.skin <- GO.plot.deg(go_deg_clean_skin, 'skin')

# ========================= blood ==============================================

deg_targets <- sapply(tissues, simplify = F, function(x){
    
    deseq.gene %>% filter(tissue == x, padj <= FDR) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

deg_target_blood <- deg_targets$blood
background_blood <- expressed_genes$blood

blood_deg_go.enrichments.res <- do.call('rbind', sapply(names(go2gene), simplify = F, function(ontology){
    
    onto.go.set <- go2gene[[ontology]]
    
    onto.results <- applyFisherGo(deg_target_blood, 
                                  background_blood, 
                                  onto.go.set, 
                                  ontology, 
                                  minGenes = 10, 
                                  maxGenes = 500)
    
    onto.results <- filter(onto.results, p.adjust <= FDR)
    
    goterms <- as.data.frame(GOTERM)[2:3] %>% 
        filter(!duplicated(go_id)) %>% 
        rename(Description = Term)
    
    
    onto.results <- merge(onto.results, as.data.frame(goterms), by = 'go_id', all.x = T)
    
    return(onto.results)
    
}))

go_deg_clean_blood <-  blood_deg_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

go_deg_clean_blood$Description <- factor(go_deg_clean_blood$Description,
                                         level = go_deg_clean_blood[order(go_deg_clean_blood$p.adjust, decreasing = T), 'Description'])

go.deg.pl.blood <- GO.plot.deg(go_deg_clean_blood, 'blood')


go.deg.pls <- ggarrange(go.deg.pl.brain, 
                        go.deg.pl.skin, 
                        go.deg.pl.blood,
                        labels = c("A", "B", "C"),
                        nrow = 3)

ggsave(filename = paste0(figure_dir, "s06_GO_term_analysis_deg.pdf"),
       plot = go.deg.pls,
       device = cairo_pdf,
       width = 10,
       height = 25,
       units = c("cm"),
       dpi = 300)


#########
# RESIS #
#########

write.csv(go_deg_clean_brain, './manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6A_brain_GO_Term_diff_genes.csv')


ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6A_brain_GO_Term_diff_genes.pdf',
    plot = go.deg.pl.brain,
    device = cairo_pdf,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 300)

write.csv(go_deg_clean_skin, './manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6B_skin_GO_Term_diff_genes.csv')


ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6B_skin_GO_Term_diff_genes.pdf',
    plot = go.deg.pl.skin,
    device = cairo_pdf,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 300)

write.csv(go_deg_clean_blood, './manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6C_blood_GO_Term_diff_genes.csv')


ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_6/S6C_blood_GO_Term_diff_genes.pdf',
    plot = go.deg.pl.blood,
    device = cairo_pdf,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 300)