# --------------------------------- Notes --------------------------------------
# 
# The purpose of that script is to do a GO-Term analysis for host genes of
# autonomous TEs.
#
#

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    source('./01_load_environment.R')
}

# ----------------------------- Functions --------------------------------------

GO.plot <- function(df){
    
    x.max <- max(-log10(df$p.adjust)) + 5
    
    p1 <- ggplot(df, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
        geom_col() +
        geom_segment(aes(x = (-log10(p.adjust)), 
                         xend = round(x.max), 
                         y = Description, 
                         yend = Description), 
                     color = 'gray50',
                     linetype  = 2) +
        theme_rob(10, base_family = 'arial') +
        facet_grid(ONTOLOGY~., scales = 'free_y', switch = "y") +
        labs(x = expression(log[10](FDR))) +
        geom_text(aes(label = Count, hjust = 1.5), color = 'white', family = 'arial', size = 6/.pt) +#*0.36) + # size from pt into mm 8 (pt) * 0.36
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
 
# -------------------------- Background genes ----------------------------------

deseq.gene <- read.csv(paste0(table_dir, '02_deseq_results_gene.csv'))

expressed.genes <- sapply(tissues, simplify = F, function(x){
 
    deseq.gene %>%
        filter(tissue == x, baseMean >= 1) %>% 
        dplyr::pull(ensembl_gene_id) %>%
        unique()
    
})

# ---------------------------- Target genes ------------------------------------
#
# We need the TE regions that do have a CAGE peak and RNA-Seq signal and 
# their associated genes. Therefore, the auto_TE_host_correlation table
# is loaded as they is based on the original definition of independent TE 
# region. The host genes are filtered by their expression, which means if they 
# are in the background or not.

data <- read.csv(paste0(table_dir, '10_auto_TE_Host_correlation.csv'))

target.genes <- sapply(tissues, simplify = F, function(x){
    
    data %>%
        filter(tissue == x, ensembl_gene_id %in% expressed.genes[[x]]) %>% 
        dplyr::pull(ensembl_gene_id) %>%
        unique()
    
})

# ------------------------------  Go Run ---------------------------------------

# ==============================   Brain =======================================

go.res.brain <- goOver(target.genes$brain, expressed.genes$brain)

goBarFacetUpdate(go.res.brain@result)

write.table(go.res.brain@result,
            file = paste0(table_dir, '11_Brain_Host_GO_results.tsv'),
            col.names = T, 
            row.names = F,
            quote = F,
            sep = '\t')

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


# --------------------- Count TEs in introns of genes --------------------------

# Determine the number of intronic TEs for each gene independent of their 
# expression signals. teRanges contains the information where each TE is 
# located (intronic, intergenic, or exonic). Hence, filter for intronic genes 
# and a subsequent count of ensembl id occurrences gives you the number of 
# intronic TEs per gene.

intronic.tes <- transGrange(teRanges) %>% 
    filter(position == 'intronic') %>% 
    group_by(ensembl_gene_id) %>% 
    dplyr::count(name = 'n.intronicTEs')

gene.intronic.tes <- 
    merge(transGrange(geneRanges), intronic.tes, by = 'ensembl_gene_id', all.x = T)

gene.intronic.tes$n.intronicTEs <- ifelse(is.na(gene.intronic.tes$n.intronicTEs), 
                                          0 , gene.intronic.tes$n.intronicTEs)


go.Genes.brain <- data_clean_brain %>% 
    group_by(ONTOLOGY) %>% 
    filter(p.adjust <= 0.05) %>% 
    dplyr::slice_min(order_by = p.adjust, n = 10) %>% 
    ungroup() %>% 
    as.data.frame()

df.intron.counts <- do.call('rbind', sapply(unique(go.Genes.brain$ID), simplify = F, function(x){
    
    cat =  go.Genes.brain %>% filter(ID == x) %>% pull(Description)
    ont =  go.Genes.brain %>% filter(ID == x) %>% pull(ONTOLOGY)
    genes <- go.Genes.brain %>% filter(ID == x) %>% pull(geneID) 
    genes <- unique(unlist(strsplit(genes, '/')))
    
    df <- gene.intronic.tes %>% 
        filter(external_gene_name %in% genes) %>% 
        mutate(category = cat,
               ontology = ont)
    
    return(df)
    
}))

backgroundGenes <- sapply(tissues, simplify = F, function(x){
    
    exp.genes <- deseq.gene %>% 
        filter(tissue == x, baseMean >= 1) %>% 
        dplyr::pull(ensembl_gene_id) 
    
    gene.intronic.tes %>% filter(ensembl_gene_id %in% exp.genes)
    
})


background <- backgroundGenes$brain

bootstrap.result <- do.call('rbind',
                            
                            sapply(unique(df.intron.counts$category), simplify = F, function(x){
                                
                                target <- df.intron.counts %>% filter(category == x)
                                
                                
                                results.target <- bootStrapping(target, 
                                                                background, 
                                                                1000)
                                
                                results.target$category <- x 
                                results.target$ontology <- unique(target$ontology)
                                results.target$group <- "Top 10 Go-Term genes"
                                
                                
                                results.random <- bootStrapping(sample_n(background, nrow(target)),
                                                                background,
                                                                1000)
                                
                                results.random$category <- x 
                                results.random$ontology <- unique(target$ontology)
                                results.random$group <- "Randomly drawn expressed genes"
                                
                                return(rbind(results.target, results.random))
                                
                            }))


intron_enrichment_brain <- ggplot(bootstrap.result, aes(category, log(ratio), fill = group)) +
    geom_boxplot( outlier.colour="gray90",
                  outlier.alpha = 0.7,
                  outlier.fill="white",
                  outlier.size=0.2,
                  notch = T) +
    geom_hline(yintercept = 0) +
    facet_grid(ontology~., scales = 'free_y') +
    coord_flip() +
    scale_fill_manual(values = c("gray90", "#e63946")) + 
    theme_rob(10, base_family = 'arial') +
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face = 'plain', size = 10),
          panel.spacing.y = unit(0.25, "lines"),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.5),
          legend.position = 'None',
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(t=0, r = 0, b = 3, l =0),'pt'),
    )

intron_enrichment_brain

ggarrange(pl_brain, 
          intron_enrichment_brain,
          widths = c(0.7, 0.3),
          align = "h")
#save(pl_brain, file = '../figures/rdata/GO.brain.Rdata')

# ==============================   skin =======================================

go.res.skin <- goOver(target.genes$skin, expressed.genes$skin)

goBarFacetUpdate(go.res.skin@result)

write.table(go.res.skin@result,
            file = paste0(table_dir, '11_Skin_Host_GO_results.tsv'),
            col.names = T, 
            row.names = F,
            quote = F,
            sep = '\t')

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
#save(pl_skin, file = '../figures/rdata/GO.skin.Rdata')

# ==============================   blood =======================================

go.res.blood <- goOver(target.genes$blood, expressed.genes$blood)

goBarFacetUpdate(go.res.blood@result)

write.table(go.res.blood@result,
            file = paste0(table_dir, '11_Blood_Host_GO_results.tsv'),
            col.names = T, 
            row.names = F,
            quote = F,
            sep = '\t')

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

#save(pl_blood, file = '../figures/rdata/GO.blood.Rdata')


# -------------------  Go Run with adapted targets -----------------------------

# ######################## Both up-regulated ###################################

target.genes <- sapply(tissues, simplify = F, function(x){
    
    data %>% filter(tissue == x, rna_quadrant == 1, 
                    ensembl_gene_id %in% expressed.genes[[x]]) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

# ==============================   Brain =======================================

go.res.brain <- goOver(target.genes$brain, expressed.genes$brain)

goBarFacetUpdate(go.res.brain@result)

# ==============================   skin =======================================

go.res.skin <- goOver(target.genes$skin, expressed.genes$skin)

goBarFacetUpdate(go.res.skin@result)

# ==============================   blood =======================================

go.res.blood <- goOver(target.genes$blood, expressed.genes$blood)

goBarFacetUpdate(go.res.blood@result)

# -------------------  Go Run with adapted targets -----------------------------
# ####################### Both down-regulated ##################################

target.genes <- sapply(tissues, simplify = F, function(x){
    
    data %>% filter(tissue == x, rna_quadrant == 3, 
                    ensembl_gene_id %in% expressed.genes[[x]]) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

# ==============================   Brain =======================================

go.res.brain <- goOver(target.genes$brain, expressed.genes$brain)

goBarFacetUpdate(go.res.brain@result)

# ==============================   skin =======================================

go.res.skin <- goOver(target.genes$skin, expressed.genes$skin)

goBarFacetUpdate(go.res.skin@result)

# ==============================   blood =======================================

go.res.blood <- goOver(target.genes$blood, expressed.genes$blood)

goBarFacetUpdate(go.res.blood@result)

# -------------------  Go Run with adapted targets -----------------------------
# ######################## Gene up & TE down ###################################

target.genes <- sapply(tissues, simplify = F, function(x){
    
    data %>% filter(tissue == x, rna_quadrant == 2, 
                    ensembl_gene_id %in% expressed.genes[[x]]) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

# ==============================   Brain =======================================

go.res.brain <- goOver(target.genes$brain, expressed.genes$brain)

goBarFacetUpdate(go.res.brain@result)

# ==============================   skin =======================================

go.res.skin <- goOver(target.genes$skin, expressed.genes$skin)

goBarFacetUpdate(go.res.skin@result)

# ==============================   blood =======================================

go.res.blood <- goOver(target.genes$blood, expressed.genes$blood)

goBarFacetUpdate(go.res.blood@result)

# -------------------  Go Run with adapted targets -----------------------------
# ######################## Gene down & TE up ###################################

target.genes <- sapply(tissues, simplify = F, function(x){
    
    data %>% filter(tissue == x, rna_quadrant == 4, 
                    ensembl_gene_id %in% expressed.genes[[x]]) %>% 
        dplyr::pull(ensembl_gene_id) %>% unique()
    
})

# ==============================   Brain =======================================

go.res.brain <- goOver(target.genes$brain, expressed.genes$brain)

goBarFacetUpdate(go.res.brain@result)

# ==============================   skin =======================================

go.res.skin <- goOver(target.genes$skin, expressed.genes$skin)

goBarFacetUpdate(go.res.skin@result)

# ==============================   blood =======================================

go.res.blood <- goOver(target.genes$blood, expressed.genes$blood)

goBarFacetUpdate(go.res.blood@result)
