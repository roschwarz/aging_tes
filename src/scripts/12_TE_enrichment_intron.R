# --------------------------------- Notes --------------------------------------
# 
# The purpose of that script is...

# Load Environment
if( !exists("ENVIRONMENT_LOADED") ){
    source('./01_load_environment.R')
} else if( !ENVIRONMENT_LOADED ){
    source('./01_load_environment.R')
}

# ------------------------------- Functions ------------------------------------

# bootStrapping <- function(n.intronic.tes.target, n.intronic.tes.background, moves = 1000){
#     
#     random.target <- sample_n(n.intronic.tes.target, moves, replace = T)[['n.intronicTEs']] 
#     
#     random.background <- sample_n(n.intronic.tes.background, moves, replace = T)[['n.intronicTEs']]
#     
#     
#     r <- (random.target+0.1)/(random.background+0.1)
#     
#     
#     result <- data.frame(n.target = random.target, 
#                          n.background = random.background, 
#                          ratio = r)
#     
#     return(result)
# }

# ------------------------  Load expressed genes -------------------------------

expressed.genes <- read.csv(paste0(table_dir, '02_deseq_results_gene.csv'))

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

# -------------------- Collect genes of interest -------------------------------

# ============================== Brain =========================================

go.Genes.brain <- read.csv(paste0(table_dir, '11_Brain_Host_GO_results.tsv'), sep ='\t')

go.Genes.brain <- go.Genes.brain %>% 
    group_by(ONTOLOGY) %>% 
    filter(p.adjust <= 0.05) %>% 
    dplyr::slice_min(order_by = p.adjust, n =10) %>% 
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
   
    exp.genes <- expressed.genes %>% 
        filter(tissue == x) %>% 
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


pl <- ggplot(bootstrap.result, aes(category, log(ratio), fill = group)) +
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
    
    

save(pl, file = '../figures/rdata/TE.intron.enrichment.brain.Rdata')

# #################### for whole GO-Term set ###################################

go.term.intron.counts <- do.call('rbind', sapply(unique(go.Genes.brain$ID), simplify = F, function(x){
    
    cat =  go.Genes.brain %>% filter(ID == x) %>% pull(Description)
    ont =  go.Genes.brain %>% filter(ID == x) %>% pull(ONTOLOGY)
    genes <- getGoGenes(x)
    #genes <- unique(unlist(strsplit(genes, '/')))
    
    df <- gene.intronic.tes %>% 
        filter(ensembl_gene_id %in% genes) %>% 
        mutate(category = cat,
               ontology = ont)
    
    return(df)
    
}))

bootstrap.GO.result <- do.call('rbind',
                            
                            sapply(unique(go.term.intron.counts$category), simplify = F, function(x){
                                
                                target <- go.term.intron.counts %>% filter(category == x)
                                
                                
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


pl.GO <- ggplot(bootstrap.GO.result, aes(category, log(ratio), fill = group)) +
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

# =============================== Skin =========================================

go.Genes.skin <- read.csv(paste0(table_dir, '11_Skin_Host_GO_results.tsv'), sep ='\t')

go.Genes.skin <- go.Genes.skin %>% 
    group_by(ONTOLOGY) %>% 
    filter(p.adjust <= 0.05) %>% 
    dplyr::slice_min(order_by = p.adjust, n =10) %>% 
    ungroup() %>% 
    as.data.frame()

df.intron.counts <- do.call('rbind', sapply(unique(go.Genes.skin$ID), simplify = F, function(x){
    
    cat =  go.Genes.skin %>% filter(ID == x) %>% pull(Description)
    ont =  go.Genes.skin %>% filter(ID == x) %>% pull(ONTOLOGY)
    genes <- go.Genes.skin %>% filter(ID == x) %>% pull(geneID) 
    genes <- unique(unlist(strsplit(genes, '/')))
    
    df <- gene.intronic.tes %>% 
        filter(external_gene_name %in% genes) %>% 
        mutate(category = cat,
               ontology = ont)
    
    return(df)
    
}))

backgroundGenes <- sapply(tissues, simplify = F, function(x){
   
    exp.genes <- expressed.genes %>% 
        filter(tissue == x) %>% 
        dplyr::pull(ensembl_gene_id) 
    
    gene.intronic.tes %>% filter(ensembl_gene_id %in% exp.genes)
    
})

    
background <- backgroundGenes$skin
    
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


pl_skin <- ggplot(bootstrap.result, aes(category, log(ratio), fill = group)) +
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
    
    

save(pl_skin, file = '../figures/rdata/TE.intron.enrichment.skin.Rdata')

# =============================== Blood ========================================

go.Genes.blood <- read.csv(paste0(table_dir, '11_Blood_Host_GO_results.tsv'), sep ='\t')

go.Genes.blood <- go.Genes.blood %>% 
    group_by(ONTOLOGY) %>% 
    filter(p.adjust <= 0.05) %>% 
    dplyr::slice_min(order_by = p.adjust, n =10) %>% 
    ungroup() %>% 
    as.data.frame()

df.intron.counts <- do.call('rbind', sapply(unique(go.Genes.blood$ID), simplify = F, function(x){
    
    cat =  go.Genes.blood %>% filter(ID == x) %>% pull(Description)
    ont =  go.Genes.blood %>% filter(ID == x) %>% pull(ONTOLOGY)
    genes <- go.Genes.blood %>% filter(ID == x) %>% pull(geneID) 
    genes <- unique(unlist(strsplit(genes, '/')))
    
    df <- gene.intronic.tes %>% 
        filter(external_gene_name %in% genes) %>% 
        mutate(category = cat,
               ontology = ont)
    
    return(df)
    
}))

backgroundGenes <- sapply(tissues, simplify = F, function(x){
   
    exp.genes <- expressed.genes %>% 
        filter(tissue == x) %>% 
        dplyr::pull(ensembl_gene_id) 
    
    gene.intronic.tes %>% filter(ensembl_gene_id %in% exp.genes)
    
})

    
background <- backgroundGenes$blood
    
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


pl_blood <- ggplot(bootstrap.result, aes(category, log(ratio), fill = group)) +
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
    
    

save(pl_blood, file = '../figures/rdata/TE.intron.enrichment.blood.Rdata')
