library(org.Mm.eg.db)
library(GO.db)
        

getBackgroundGenes <- function(deseq.res){
    # Takes a deseq result table and filters for expressed genes
    # A gene is considered as expressed when the adjusted p value is not NA
    background.genes <- as.data.frame(deseq.res)
    background.genes <- background.genes[!is.na(background.genes$padj),]
    background.genes <- unique(rownames(background.genes))
    
    return(background.genes)
}


goOver <- function(target, background, onto = 'ALL'){
    
    require(clusterProfiler)
    if(length(target) == 0 || length(background) == 0){
        warning('Return Null as target gene or background gene set are empty!')
        return(NULL)
    }
    
    enrichGO(
        gene = target,
        universe = background,
        OrgDb = 'org.Mm.eg.db',
        ont = onto,
        keyType = "ENSEMBL",
        pAdjustMethod = 'BH',
        #pvalueCutoff = 0.05,
        #qvalueCutoff = 0.05,
        readable = TRUE
    )
    
}


getGoGenes <- function(goTerm){
    # Takes a GO-Term and returns a list of genes that are associated to this
    # set. CURRENTLY only mouse specific.
    library(org.Mm.eg.db)
    
    goAnnotation <- AnnotationDbi::select(org.Mm.eg.db, 
                                       keytype="GOALL", 
                                       keys=goTerm, 
                                       columns=c("ENSEMBL"))
    
    genes <- goAnnotation %>% dplyr::pull(ENSEMBL) %>% unique()
    
    return(genes)
    
}


functionalGO <- function(genes_of_interest, 
                         background_genes, 
                         species = "Mus musculus",
                         ontologie = "BP"){
    
    
    if(class(genes_of_interest) != "character" || class(background_genes) != "character") {
        stop("Your gene lists are not in the correct format. A vector of characters is needed.")
    }
    
    dbs <- list("Mus musculus" = "org.Mm.eg.db", 
                "Homo sapiens" = "org.Hs.eg.db")
    
    if(!(species %in% names(dbs))){
        stop("Data base of your species is not contained in the dbs list. Add
             the specific db to your list in functionalGO().")
    }
    
    require(clusterProfiler, quietly = T)
    
    switch (species,
            "Mus musculus" = require(org.Mm.eg.db),
            "Homo sapiens" = require(org.Hs.eg.db)
    )
    
    
    
    enrichGO(gene = genes_of_interest, 
             universe = background_genes,
             keyType = "ENSEMBL",
             OrgDb = dbs[species][[1]], 
             ont = ontologie, 
             pAdjustMethod = "BH", 
             qvalueCutoff = 0.05, 
             readable = TRUE)
    
}

# dotplot arise error if there is no significant data point contained
# I catch this error before with a statement that makes sense.
goPlots <- function(ego, res = ".", prefix = "", sig = 0.05){
    
    require(enrichplot, quietly = T)
    require(tidyverse)
    
    if(!file.exists(paste0(res, "/figures"))){
        dir.create(paste0(res, "/figures"))
    }
    
    number.of.sig <- ego@result %>% filter(p.adjust <= sig) %>% nrow()
    
    if(number.of.sig == 0){
        stop('No significant data point contained.')
    }
    
    png(filename = paste0(res, "/figures/", prefix, ".go.dotplot.png"), width = 1500, height = 1200)
    show(dotplot(ego, showCategory=50))
    dev.off()
    
    ego <- enrichplot::pairwise_termsim(ego)
    
    png(filename = paste0(res, "/figures/", prefix, ".go.relationship.png"), width = 1800, height = 1600)
    show(emapplot(ego, showCategory=50))
    dev.off()
    
    
}

enrichGO_for_list <- function(list.of.target.genes, background.genes){
    
    enrichtGOresult = list()
    
    for(ontology in c('MF', 'BP', 'CC')){
        
        enrichtGOresult[[ontology]] <- sapply(tissues, 
                                              simplify = F, 
                                              function(x){ 
                                                  
                                                  ego1 <- gseGO(geneList = te.associated.genes.sorted[[x]],
                                                                OrgDb = 'org.Mm.eg.db',
                                                                ont = ontology,
                                                                keyType = "ENTREZID",
                                                                minGSSize    = 10,
                                                                maxGSSize    = 500,
                                                                pvalueCutoff = 0.05,
                                                                verbose      = FALSE)
                                                  
                                                  if(nrow(ego1) == 0){
                                                      return(NULL)}
                                                  
                                                  return(ego1)
                                                  
                                              }
                                              
        )
        
    }
    
    return(enrichtGOresult)
    
    
}

gOTEassoGenes <- function(teRanges, geneRanges, cageRanges, teDict){
    
   #teDict is a table with information about association to genes 
    
   # TEs that share a peak with an exon are removed!!!
    
    
   # intersect cage peaks with gene and TE annotation
   # determine target and background gene set
   # run GO-Term analysis
   # return the GO term result
    
    expressedGenes <- transGrange(intersectGranger(geneRanges,
                                                   cageRanges,
                                                   tab = 'all'))
    
    # get TEs that intersect with at least one CAGE peak
    expressedTEs <- transGrange(intersectGranger(teRanges, 
                                                 cageRanges, 
                                                 tab='all'))
    
    # remove duplicated TE ids
    expressedTEs <- expressedTEs %>% 
        #filter(!(subject.peak_name %in% expressedGenes$subject.peak_name)) %>% # remove peaks that are shared with exons
        dplyr::select(query.te.id) %>% 
        filter(!duplicated(query.te.id)) %>% 
        pull(query.te.id)
    
    target.genes <- teDict %>% 
        filter(TE.ID %in% expressedTEs, distance <= 10000) %>% #filter for expressed TEs; gene needs to be within 10kb
        dplyr::select(TE.ID, gene.ID) %>% 
        filter(gene.ID %in% expressedGenes$query.ensembl_gene_id) %>%  #filter TEs with expressed gene
        filter(!duplicated(gene.ID)) %>% 
        dplyr::pull(gene.ID)
    
    background.genes <- expressedGenes %>% 
        dplyr::pull(query.ensembl_gene_id) %>% 
        unique()
    
    res <- goOver(target.genes, background.genes)
    
    return(res) 
}


# ========================== Own approach ======================================
# Approach is based in Fisher's test.

# loads the go terms associated ensembl ids separately for the three ontologies.
load_GO_Term_Annotation <- function(column='ENSEMBL', ontology = c('BP', 'MF', 'CC')){
    
    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    
    
    goTerm_genes <- sapply(ontology, simplify = F, function(x){
        
        
        terms <- goterms[goterms == x]   
        
        go2gene <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                         keys=names(terms),
                                         column=column,
                                         keytype = 'GOALL',
                                         multiVals = 'list')
        
        return(go2gene)
        
    }) 
    
    return(goTerm_genes)
}


# Function that calculates a contingency matrix that is needed for the fisher's 
# test
# the following function takes a target gene set, background gene set and 
# the set of genes contained in a specific go term.
create_GO_contingency_tab <- function(target, background, go_gene_set){
    
    target_gene_in_go <- length(intersect(go_gene_set, target))
    target_gene_not_in_go <- length(target) - target_gene_in_go
    
    # Is that correct to remove the target genes from the background?
    # I think so, considering the contingency matrix it doesn't make sense
    # to be in both groups.
    background_gene_set <- setdiff(background, target) 
    
    background_in_go = length(intersect(go_gene_set, background_gene_set))
    background_not_in_go = length(background_gene_set) - background_in_go
    
    contingency_tab <- matrix(c(target_gene_in_go, target_gene_not_in_go, 
                               background_in_go, background_not_in_go), 
                             nrow = 2) 
    
    return(contingency_tab)
}


# Needs as input the contingency matrix and the go term.
# The function distinguish between not significant and significant, but the 
# significant parts are separated into enriched and depleted.
# returns a data frame that contains following values:
# - go_term
# - p-value
# - odds_ration
# - enrichment information

GO_fisher <- function(contingency_tab, go_term_id){
    
    res <- fisher.test(contingency_tab)
    
    if(res$p.value > 0.05){
        
        go_res <- data.frame(go_id = go_term_id,
                             p_value = res$p.value, 
                             odds_ratio = res$estimate,
                             enriched = NA_character_)
        
    }else{
        
        # less then expected
        res_depleted <- fisher.test(contingency_tab, alternative = 'less')
        
        if(res_depleted$p.value <= 0.05){
            
            go_res <- data.frame(go_id = go_term_id, 
                                 p_value = res_depleted$p.value, 
                                 odds_ratio = res_depleted$estimate,
                                 enriched = 'depleted')
            
        }else{
            
            # more then expected
            res_enriched <- fisher.test(contingency_tab, alternative = 'greater')
            
            go_res <- data.frame(go_id = go_term_id, 
                                 p_value = res_enriched$p.value, 
                                 odds_ratio = res_enriched$estimate,
                                 enriched = 'enriched')
            
        }
    }
    
    return(go_res)
    
}

applyFisherGo <- function(target_gene_set, background_gene_set, go_set, ontology, minGenes = 10, maxGenes = 500){
    
    onto.results <- do.call('rbind', sapply(names(go_set), simplify = F, function(go_term){
        
        if(dplyr::between(length(go_set[[go_term]]), minGenes,maxGenes)){
            
            go_gene_set <- go_set[[go_term]]
            
            contingency_tab <- create_GO_contingency_tab(target_gene_set, background_gene_set, go_gene_set)
            
            df <- GO_fisher(contingency_tab, go_term)
            
            df$ONTOLOGY <- ontology
            df$target.count <- contingency_tab[1,1]
            df$go.count <- length(go_gene_set)
            
            return(df)
            
        }
        
        return(NULL)
    }))
    
    onto.results$p.adjust <- p.adjust(onto.results$p_value, method = "fdr")
    
    return(onto.results)
}
