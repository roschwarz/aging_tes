# Get a set of DEGs that are associated with expressed TEs in their 
# environment (restricted by distance; default 10 kb). 
# Genes are duplicated in 
# the data frame, since multiple TEs are associated with the same gene.
# To each gene the information is added how much TEs are expressed in their 
# environment and the mean and median log2FC of these. The information to the 
# respective genes are stored in the columns te.median.log2FC, te.mean.log2FC, 
# number.of.TE. Additionally, the mean and median distance is added to the 
# table.
# Finally, only one line per gene is returned.
#
# Updates:
#
# The distance and gene information are already added to the TE deseq result
# table, so that the merge with the closest table isn't anymore necessary.
#
# The distance.threshold is now set by the distance argument, as there was no
# filtering when I used the same argument name as the column name.
# filter(distance <= distance) seems not to work.
#
# Add the TE.IDs of the associated TEs in a column (TEs) separated by a ';'.

geneTEAsso <- function(te.dres, gene.dres, closest=NULL, distance = 10000){
    
    require(tidyverse)
    
    distance.threshold <- distance
    
    if(anyNA(te.dres$padj)){
        te.dres <- te.dres %>% filter(!is.na(padj))
    }
    
    if(class(te.dres) != 'data.frame'){
        te.dres <- as.data.frame(te.dres)
    }
    
    if(class(gene.dres) != 'data.frame'){
        gene.dres <- as.data.frame(gene.dres)
    }
    
    gene.te.asso <- te.dres %>% 
        filter(distance <= distance.threshold) %>%
        group_by(ensembl_gene_id) %>% 
        mutate(te.median.log2FC = mean(log2FoldChange), 
               te.mean.log2FC = median(log2FoldChange),
               te.median.distance = median(distance),
               te.mean.distance = mean(distance),
               number.of.TE = n(),
               TEs = paste(TE.ID, collapse = ';')) %>% 
        dplyr::select(ensembl_gene_id, 
                      te.median.log2FC, 
                      te.mean.log2FC, 
                      te.median.distance,
                      te.mean.distance,
                      number.of.TE,
                      TEs) %>% 
        unique()
    
    degs <- gene.dres %>% filter(padj <= 0.05)
    
    if(!'ensembl_gene_id' %in% names(gene.dres)){
        degs <- degs %>% rownames_to_column('ensembl_gene_id')
    }
    
    
    gene.te.asso <- merge(degs, gene.te.asso, by = 'ensembl_gene_id')
    
    gene.te.asso$gene.exp.direction <-  ifelse(gene.te.asso$log2FoldChange < 0, 'down', 'up')
    
    gene.te.asso <- filter(gene.te.asso, !is.na(ensembl_gene_id))
    
    return(gene.te.asso)
    
}
