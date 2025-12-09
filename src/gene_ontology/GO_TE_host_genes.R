# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

load_analysis_env()
load_annotations()
load_plotting_env()

load_te_ranges()
load_gene_ranges()


# ======================== Functions ===========================================

GO.plot <- function(df){
    
    x.max <- round(max(-log10(df$p.adjust))) + 1
    x.max <- round(x.max, -1) + 5
    
    p1 <- ggplot(df, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
        geom_segment(aes(x=0, 
                         xend = x.max, 
                         y=Description, 
                         yend=Description), 
                     color = 'gray50',
                     linetype  = 2) +
        geom_col() +
        #theme_rob(8) +
        facet_grid(ONTOLOGY~., scales = 'free_y', switch = "y") +
        labs(x = expression(log[10](FDR))) +
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
              axis.text = element_text(size = 8),
              axis.title = element_text(size = 10),
              strip.text.y = element_text(angle = 360),
              legend.position = 'None',
              strip.placement = "outside",
              strip.background = element_rect(color = 'white'),
              strip.text = element_text(size = 10),
              plot.margin = unit(c(t=0, r = 0, b = 0, l =0),'pt'))
    
    return(p1)
}


intron.plot <- function(df){
    
    ggplot(df, aes(category, log(ratio), fill = group)) +
        geom_boxplot( outlier.colour="gray90",
                      outlier.alpha = 0.7,
                      outlier.fill="white",
                      outlier.size=0.2,
                      notch = T) +
        geom_hline(yintercept = 0) +
        facet_grid(ontology~., scales = 'free_y') +
        coord_flip() +
        scale_fill_manual(values = c("#e63946", "gray90")) + 
        #theme_rob(8) +
        theme(panel.grid = element_blank(),
              axis.title.y = element_blank(),
              #axis.title.x = element_text(face = 'plain', size = 8),
              panel.spacing.y = unit(0.25, "lines"),
              panel.border = element_blank(),
              axis.text = element_text(size = 8),
              axis.title = element_text(size = 10),
              axis.line = element_line(size = 0.5),
              legend.position = 'None',
              strip.background = element_blank(),
              strip.text = element_blank(),
              axis.text.y = element_blank(),
              plot.margin = unit(c(t=0, r = 0, b = 3, l =0),'pt'),
        )
    
} 

# ==================== Genes with independent TE ===============================

go2gene <- load_GO_Term_Annotation() # load goTerms + Ensembl ids

deseq.gene <- read.csv(paste0(table_dir, '02_deseq_results_genes.csv'))


expressed_genes <- sapply(tissues, simplify = F, function(x){
    
    deseq.gene %>% filter(tissue == x) %>% 
        dplyr::pull(ensembl_gene_id) %>% 
        unique()
    
})

data <- read.csv(paste0(table_dir, 'indie_te_island_host_gene_correlation.csv'))

target.genes <- sapply(tissues, simplify = F, function(x){
    
    data %>% filter(tissue == x, ensembl_gene_id %in% expressed_genes[[x]]) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

# ========================== Brain =============================================

target_brain <- target.genes$brain
background_brain <- expressed_genes$brain


brain_go.enrichments.res <- do.call('rbind', sapply(names(go2gene), simplify = F, function(ontology){
    
    onto.go.set <- go2gene[[ontology]]
    
    onto.results <- applyFisherGo(target_brain, 
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
    
}))

data_clean_brain <-  brain_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_brain$Description <- factor(data_clean_brain$Description,
                                       level = data_clean_brain[order(data_clean_brain$p.adjust, decreasing = T), 'Description'])

go.pl.brain <- GO.plot(data_clean_brain)


# ==================== TE enrichment in introns of genes ========================

intronic.tes <- transGrange(teRanges) %>% 
    filter(position == 'intronic') %>% 
    group_by(ensembl_gene_id) %>% 
    dplyr::count(name = 'n.intronicTEs')

gene.intronic.tes <- 
    merge(transGrange(geneRanges), 
          intronic.tes, by = 'ensembl_gene_id', all.x = T)

gene.intronic.tes$n.intronicTEs <- ifelse(is.na(gene.intronic.tes$n.intronicTEs), 
                                          0 , gene.intronic.tes$n.intronicTEs)

go.terms.brain <- brain_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    filter(p.adjust <= 0.05) %>% 
    dplyr::slice_min(order_by = p.adjust, n =10, with_ties = FALSE) %>% 
    ungroup()


brain_go.term.intron.counts <- do.call('rbind', sapply(unique(go.terms.brain$go_id), simplify = F, function(x){
    
    cat =  go.terms.brain %>% 
        filter(go_id == x) %>% 
        pull(Description)
    
    
    ont =  go.terms.brain %>% 
        filter(go_id == x) %>% 
        pull(ONTOLOGY)
    
    genes <- go2gene[[ont]][[x]]
    
    df <- gene.intronic.tes %>% 
        filter(ensembl_gene_id %in% genes) %>% 
        mutate(go_id = x,
               category = cat,
               ontology = ont)
    
    return(df)
    
}))

brain_background_set <- gene.intronic.tes %>% 
    filter(ensembl_gene_id %in% expressed_genes$brain)


brain_bootstrap.GO.result <- do.call('rbind',
                                     
                                     sapply(unique(brain_go.term.intron.counts$go_id), simplify = F, function(x){
                                         
                                         target <- brain_go.term.intron.counts %>% 
                                             filter(go_id == x)
                                         
                                         category = unique(target$category)
                                         ontology = unique(target$ontology)
                                         
                                         results.target <- bootStrapping(target, 
                                                                         brain_background_set, 
                                                                         1000)
                                         
                                         results.target$go_id <- x 
                                         results.target$category <- category
                                         results.target$ontology <- ontology
                                         results.target$group <- "GO term genes"
                                         
                                         results.random <- bootStrapping(sample_n(brain_background_set, nrow(target)),
                                                                         brain_background_set,
                                                                         1000)
                                         
                                         results.random$go_id <- x 
                                         results.random$category <- category
                                         results.random$ontology <- unique(target$ontology)
                                         results.random$group <- "Randomly sampled expressed genes"
                                         
                                         
                                         
                                         return(rbind(results.target, results.random))
                                         
                                     }))


brain_bootstrap.GO.result$category <- factor(brain_bootstrap.GO.result$category, 
                                             levels = data_clean_brain[order(data_clean_brain$p.adjust, decreasing = T), 'Description'])



intron.pl.brain <- intron.plot(brain_bootstrap.GO.result)

pl_brain <- ggarrange(go.pl.brain, 
          intron.pl.brain, 
          widths = c(0.7, 0.3),
          align = 'h')

meta <- list(name = 'GO_TE_intron_enrichment_brain',
             description = 'GO term enrichment analysis for genes containing independently expressed TE islands in their introns in brain tissue. Enrichment of TEs in introns.',
             tags = c('GO enrichment', 'TE islands', 'intronic TEs', 'brain'),
             parameters = list(FDR = 0.05, tissues = c('brain'), sex = c('male')),
             script = 'GO_TE_host_genes.R'
)

fig_index(plot = pl_brain,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 18,
          height = 12,
          dpi = 300,
          format = 'pdf')


# ========================== Skin =============================================

target_skin <- target.genes$skin
background_skin <- expressed_genes$skin


skin_go.enrichments.res <- do.call('rbind', sapply(names(go2gene), simplify = F, function(ontology){
    
    onto.go.set <- go2gene[[ontology]]
    
    onto.results <- applyFisherGo(target_skin, 
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

data_clean_skin <- skin_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_skin$Description <- factor(data_clean_skin$Description,
                                      level = data_clean_skin[order(data_clean_skin$p.adjust, decreasing = T), 'Description'])

go.pl.skin <- GO.plot(data_clean_skin)


go.terms.skin <- skin_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    filter(p.adjust <= 0.05) %>% 
    dplyr::slice_min(order_by = p.adjust, n =10, with_ties = FALSE) %>% 
    ungroup()


skin_go.term.intron.counts <- do.call('rbind', sapply(unique(go.terms.skin$go_id), simplify = F, function(x){
    
    cat =  go.terms.skin %>% 
        filter(go_id == x) %>% 
        pull(Description)
    
    
    ont =  go.terms.skin %>% 
        filter(go_id == x) %>% 
        pull(ONTOLOGY)
    
    genes <- go2gene[[ont]][[x]]
    
    df <- gene.intronic.tes %>% 
        filter(ensembl_gene_id %in% genes) %>% 
        mutate(go_id = x,
               category = cat,
               ontology = ont)
    
    return(df)
    
}))

skin_background_set <- gene.intronic.tes %>% 
    filter(ensembl_gene_id %in% expressed_genes$skin)


skin_bootstrap.GO.result <- do.call('rbind',
                                    
                                    sapply(unique(skin_go.term.intron.counts$go_id), simplify = F, function(x){
                                        
                                        target <- skin_go.term.intron.counts %>% 
                                            filter(go_id == x)
                                        
                                        category = unique(target$category)
                                        ontology = unique(target$ontology)
                                        
                                        results.target <- bootStrapping(target, 
                                                                        skin_background_set, 
                                                                        1000)
                                        
                                        results.target$go_id <- x 
                                        results.target$category <- category
                                        results.target$ontology <- ontology
                                        results.target$group <- "GO term genes"
                                        
                                        results.random <- bootStrapping(sample_n(skin_background_set, nrow(target)),
                                                                        skin_background_set,
                                                                        1000)
                                        
                                        results.random$go_id <- x 
                                        results.random$category <- category
                                        results.random$ontology <- unique(target$ontology)
                                        results.random$group <- "Randomly sampled expressed genes"
                                        
                                        
                                        
                                        return(rbind(results.target, results.random))
                                        
                                    }))


skin_bootstrap.GO.result$category <- factor(skin_bootstrap.GO.result$category, 
                                            levels = data_clean_skin[order(data_clean_skin$p.adjust, decreasing = T), 'Description'])



intron.pl.skin <- intron.plot(skin_bootstrap.GO.result)


supplemental_skin <- annotate_figure(ggarrange(go.pl.skin, 
                                               intron.pl.skin, 
                                               widths = c(0.7, 0.3),
                                               align = 'h'),
                                     top = text_grob("skin", size = 12))

supplemental_skin

# ========================== blood =============================================

target_blood <- target.genes$blood
background_blood <- expressed_genes$blood


blood_go.enrichments.res <- do.call('rbind', sapply(names(go2gene), simplify = F, function(ontology){
    
    onto.go.set <- go2gene[[ontology]]
    
    onto.results <- applyFisherGo(target_blood, 
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

data_clean_blood <-  blood_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(ONTOLOGY = case_when(ONTOLOGY == 'BP' ~ 'Biological\nProcess',
                                ONTOLOGY == 'CC' ~ 'Cellular\nComponent',
                                TRUE ~ 'Molecular\nFunction'))

data_clean_blood$Description <- factor(data_clean_blood$Description,
                                       level = data_clean_blood[order(data_clean_blood$p.adjust, decreasing = T), 'Description'])

go.pl.blood <- GO.plot(data_clean_blood)


# ==================== TE enrichment in intron of genes ========================

# intronic.tes <- transGrange(teRanges) %>% 
#     filter(position == 'intronic') %>% 
#     group_by(ensembl_gene_id) %>% 
#     dplyr::count(name = 'n.intronicTEs')
# 
# gene.intronic.tes <- 
#     merge(transGrange(geneRanges), 
#           intronic.tes, by = 'ensembl_gene_id', all.x = T)
# 
# gene.intronic.tes$n.intronicTEs <- ifelse(is.na(gene.intronic.tes$n.intronicTEs), 
#                                           0 , gene.intronic.tes$n.intronicTEs)

go.terms.blood <- blood_go.enrichments.res %>% 
    filter(enriched == 'enriched') %>% 
    group_by(ONTOLOGY) %>% 
    filter(p.adjust <= 0.05) %>% 
    dplyr::slice_min(order_by = p.adjust, n =10, with_ties = FALSE) %>% 
    ungroup()


blood_go.term.intron.counts <- do.call('rbind', sapply(unique(go.terms.blood$go_id), simplify = F, function(x){
    
    cat =  go.terms.blood %>% 
        filter(go_id == x) %>% 
        pull(Description)
    
    
    ont =  go.terms.blood %>% 
        filter(go_id == x) %>% 
        pull(ONTOLOGY)
    
    genes <- go2gene[[ont]][[x]]
    
    df <- gene.intronic.tes %>% 
        filter(ensembl_gene_id %in% genes) %>% 
        mutate(go_id = x,
               category = cat,
               ontology = ont)
    
    return(df)
    
}))

blood_background_set <- gene.intronic.tes %>% 
    filter(ensembl_gene_id %in% expressed_genes$blood)


blood_bootstrap.GO.result <- do.call('rbind',
                                     
                                     sapply(unique(blood_go.term.intron.counts$go_id), simplify = F, function(x){
                                         
                                         target <- blood_go.term.intron.counts %>% 
                                             filter(go_id == x)
                                         
                                         category = unique(target$category)
                                         ontology = unique(target$ontology)
                                         
                                         results.target <- bootStrapping(target, 
                                                                         blood_background_set, 
                                                                         1000)
                                         
                                         results.target$go_id <- x 
                                         results.target$category <- category
                                         results.target$ontology <- ontology
                                         results.target$group <- "GO term genes"
                                         
                                         results.random <- bootStrapping(sample_n(blood_background_set, nrow(target)),
                                                                         blood_background_set,
                                                                         1000)
                                         
                                         results.random$go_id <- x 
                                         results.random$category <- category
                                         results.random$ontology <- unique(target$ontology)
                                         results.random$group <- "Randomly sampled expressed genes"
                                         
                                         
                                         
                                         return(rbind(results.target, results.random))
                                         
                                     }))


blood_bootstrap.GO.result$category <- factor(blood_bootstrap.GO.result$category, 
                                             levels = data_clean_blood[order(data_clean_blood$p.adjust, decreasing = T), 'Description'])



intron.pl.blood <- intron.plot(blood_bootstrap.GO.result)


supplemental_blood <- annotate_figure(ggarrange(go.pl.blood, 
                                                intron.pl.blood, 
                                                widths = c(0.7, 0.3),
                                                align = 'h'),
                                      top = text_grob("blood", size = 12))

supplemental_go_te_containing_ges <- ggarrange(supplemental_skin,
                            supplemental_blood,
                            labels = c("A","B"),
                            nrow = 2,
                            ncol = 1)

supplemental_meta <- list(name = 'GO_TE_intron_enrichment_skin_blood',
             description = 'GO term enrichment analysis for genes containing independently expressed TE islands in their introns in skin and blood tissue. Enrichment of TEs in introns.',
             tags = c('GO enrichment', 'TE islands', 'intronic TEs', 'skin', 'blood'),
             parameters = list(FDR = 0.05, tissues = c('skin', 'blood'), sex = c('male')),
             script = 'GO_TE_host_genes.R'
)

fig_index(plot = supplemental_go_te_containing_ges,
          outdir = figure_dir,
          meta = supplemental_meta,
          index_file = 'figure_index.tsv',
          width = 18,
          height = 28,
          dpi = 300,
          format = 'pdf')

capture.output(sessionInfo(), file = paste0(figure_dir, 'GO_TE_host_genes_sessionInfo.txt'))
