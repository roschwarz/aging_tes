# --------------------------------- Notes--- -----------------------------------
# This script is used for the RNA-Seq quantification
#
# Input:
#
# - Count tables of the respective tissue (EXPR.csv) in ./results/rna_RS/salmonTE/<tissue>_deduplicated/
# - te region file for deseq analysis for TE islands 
#
# Output data (../data/rna; ../tables):
#
# TE instances
#
# dds.TE.Salmon.Rdata - DESeq2 dds object for TEs base on SalmonTE
# deseq.TE.SalmonTE.Rdata - DESeq2 results for TEs based on SalmonTE counts
# 02_deseq_results_te.csv - tissue merged csv file of DESeq results of TE
#   instances
#
# TE region
#
# dds.TE.region.Salmon.Rdata - DESeq2 dds object for TE regions base on SalmonTE
# deseq.TE.region.SalmonTE.Rdata - DESeq2 results for TE regions based on
#       SalmonTE counts
# 02_deseq_results_te_region.csv - tissue merged csv file of DESeq results of TE regions. The coordinates of
# each TE regions is added and the information if a TE region is differentially expressed and in which region.
# The .csv table can be easy filter for DETEs, e.g. with grep <tissue>,TRUE,up 02_deseq_results_re_region.csv
#
#
# Genes
#
# A transcriptome was added to the salmonTE index, hence the salmonTe count
# tablecontains also counts for the gene isoforms. For the gene part the counts
# of each isoform are summed up.
#
# dds.gene.Salmon.Rdata - DESeq2 dds object for genes base on SalmonTE
# deseq.genes.SalmonTE.Rdata - DESeq2 results for genes based on SalmonTE counts
# 02_deseq_results_genes.csv - tissue merged csv file of DESeq results of genes
#
# Gene isoforms
#
# dds.gene.isoform.Salmon.Rdata - DESeq2 dds object for gene isoforms base on SalmonTE
# deseq.gene.isoform.SalmonTE.Rdata - DESeq2 results for gene isoforms based on SalmonTE counts
# 02_deseq_results_gene_isoform.csv - tissue merged csv file of DESeq results of gene isoforms
#

# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {

    source("./01_load_environment.R")

}

getConditions <- function(file_names) {
# Project specific. It uses the naming convention of the fastq file to
# determine the tissue and age of each sample.
    
    
    samples <- data.frame(
        SampleID = file_names,
        condition = file_names)
    
    # define the conditions for each sample
    if (grepl("skin", file_names[1], ignore.case = T)) {
        samples <- samples %>%
            separate(condition,
                     c("ID", "age"),
                     sep = "[_]") %>%
            mutate(age = str_replace(age, "[0-9]*Skin", "")) %>%
            mutate(condition = as.factor(case_when(
                age == "MO" ~ "old",
                age == "MY" ~ "young"
            ))) %>%
            dplyr::select(SampleID, condition) %>% 
        column_to_rownames("SampleID")
        
        return(samples)
    }
    
    samples <- samples %>%
        separate(condition,
                 c("ID", "org", "age", "tissue", "number"),
                 sep = "[_]") %>%
        mutate(condition = as.factor(case_when(age == "o" ~ "old",
                                               age == "y" ~ "young"))) %>%
        dplyr::select(SampleID, condition) %>% 
        column_to_rownames("SampleID")
    
    
    return(samples)
    
    
}

# ----- Quantification of TE expression during aging in different tissues ------

if (!file.exists(paste0(data_dir, "rna_RS_deseq_TE_instances_SalmonTE.Rdata"))) {

    # read the count tables and store them in a list
    # path to count tables can be found in 01_load_environment
    count.tables.rna <- sapply(names(counts_rna_RS), simplify = FALSE,
                               function(x) {
                                   
        load_salmonTE_counts(counts_rna_RS[[x]],
                             "instance")[["instance"]] %>% 
                                       tibble::column_to_rownames("TE")


    })

    # assign the sample names to a condition (age)
    conditions <- sapply(names(counts_rna_RS), simplify = FALSE,
                         function(x) {

            getConditions(names(count.tables.rna[[x]]))

    })

    # runDESeq
    dds.te <- sapply(names(count.tables.rna), simplify = FALSE, function(x) {

        doDEseq(count.matrix = round(count.tables.rna[[x]]),
                col_data = conditions[[x]],
                reference = "young",
                target = "te"
        )
    })

    # collect DESeq results
    deseq.te <- sapply(names(dds.te), simplify = FALSE, function(x) {

        getDEseqResults(dds.te[[x]],
                        coefficient = "condition_old_vs_young",
                        FDR.filter = FALSE)

    })

    save(dds.te, file = paste0(data_dir,
                               "rna_RS_dds_TE_instances_SalmonTE.Rdata"))
    save(deseq.te, file = paste0(data_dir,
                                 "rna_RS_deseq_TE_instances_SalmonTE.Rdata"))

    rm(dds.te)

}else{

    deseq.te <- loadRdata(paste0(data_dir, "rna_RS_deseq_TE_instances_SalmonTE.Rdata"))
    
}


# Merge the results of the different tissues to get one huge table with all results
deseq.te.merged <- do.call("rbind",
                           sapply(names(deseq.te), simplify = F, function(x) {
                               deseq.te[[x]] %>%
                                   rownames_to_column(var = "te_id") %>%
                                   mutate(tissue = x)

                           }))

deseq.te.merged <- merge(deseq.te.merged, te.annotation, by = 'te_id')

write.table(deseq.te.merged,
            file = paste0(table_dir, '02_deseq_results_te_instances.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

rm(deseq.te.merged, deseq.te)

# ----- Quantification of summed TE region expression during aging in different tissues ---

if (!file.exists(paste0(data_dir, "rna_RS_deseq_TE_region_SalmonTE.Rdata"))) {
    
    TE.region.instances <- read.csv(te.region.instances.file)
    
    count.tables.rna <-
        sapply(names(counts_rna_RS), simplify = FALSE,
               function(x) {
                   load_salmonTE_counts(counts_rna_RS[[x]],
                                        "region",
                                        TE.region.instances)[["region"]]
                   
               })
    
    
    
    # assign the sample names to a condition (age)
    conditions <- sapply(names(counts_rna_RS), simplify = FALSE,
                         function(x) {

            getConditions(names(count.tables.rna[[x]]))

    })
    
    # runDESeq
    dds.te.region <- sapply(names(count.tables.rna), simplify = F,  function(x){
        
        doDEseq(count.matrix = round(count.tables.rna[[x]]),
                col_data = conditions[[x]],
                reference = "young",
                target = "all"
        )
    })
    
    # collect DESeq results
    deseq.te.region <- sapply(names(dds.te.region), simplify = F, function(x){
        
        getDEseqResults(dds.te.region[[x]],
                        coefficient = "condition_old_vs_young",
                        FDR.filter = F)
        
    })
    
    
    save(dds.te.region, file = paste0(data_dir,
                               "rna_RS_dds_TE_region_SalmonTE.Rdata"))
    
    save(deseq.te.region, file = paste0(data_dir,
                                 "rna_RS_deseq_TE_region_SalmonTE.Rdata"))

    rm(dds.te.region)
    
}else{
    
    deseq.te.region <- loadRdata(paste0(data_dir,
                                 "rna_RS_deseq_TE_region_SalmonTE.Rdata"))
    
} 

teRegionCoordinates <- data.frame(teRegionRanges)

deseq.te.region.merged <- 
    do.call("rbind",
            sapply(names(deseq.te.region),
                   simplify = F,
                   function(x) {
                       
                       tmp_df <- deseq.te.region[[x]] %>%
                           rownames_to_column(var = "te_region_id") %>%
                           mutate(tissue = x,
                                  diff = case_when(padj <= FDR ~ TRUE, .default = FALSE),
                                  exp_dir = case_when(log2FoldChange < 0 ~ "down", 
                                                      log2FoldChange > 0 ~ "up",
                                                      .default = "None"))
                       
                       tmp_df <- merge(tmp_df,
                                       teRegionCoordinates, 
                                       by.x = "te_region_id", by.y = 'names')
                       
                       tmp_df <- tmp_df %>% 
                           mutate(te_region_coordinates = paste0(seqnames, ":", start, "-", end))
                       
                       return(tmp_df)
                       
                       
                   }))

write.table(deseq.te.region.merged, 
            file = paste0(table_dir, "02_deseq_results_te_region.csv"), 
            col.names = TRUE, 
            row.names = FALSE,
            quote = FALSE,
            sep = ',')

rm(deseq.te.region.merged, deseq.te.region, teRegionCoordinates)


# ----- Quantification of gene expression during aging in different tissues ----

if (!file.exists(paste0(data_dir, "rna_RS_deseq_gene_SalmonTE.Rdata"))) {
    
    count.tables.rna <-
        sapply(names(counts_rna_RS), simplify = FALSE,
               function(x) {
                   load_salmonTE_counts(counts_rna_RS[[x]],
                                        "gene")[["gene"]]
                   
               })
    
    
    # assign the sample names to a condition (age)
    conditions <- sapply(names(counts_rna_RS), simplify = FALSE,
                         function(x) {

            getConditions(names(count.tables.rna[[x]]))

    })
    
    # runDESeq
    dds.gene <- sapply(names(count.tables.rna), simplify = F,  function(x){
        
        doDEseq(count.matrix = round(count.tables.rna[[x]]),
                col_data = conditions[[x]],
                reference = 'young',
                target = 'all'
        )
    })
    
    # collect DESeq results
    deseq.gene <- sapply(names(dds.gene), simplify = F, function(x){
        
        getDEseqResults(dds.gene[[x]], coefficient = 'condition_old_vs_young', FDR.filter = F)
        
    })
    
    save(dds.gene,
         file = paste0(data_dir, "rna_RS_dds_gene_SalmonTE.Rdata"))
    save(deseq.gene,
         file = paste0(data_dir, "rna_RS_deseq_gene_SalmonTE.Rdata"))
    
    rm(dds.gene)
    
}else{
    
    deseq.gene <- loadRdata(paste0(data_dir,
                                   "rna_RS_deseq_gene_SalmonTE.Rdata"))
}


deseq.gene.merged <- do.call('rbind', sapply(names(deseq.gene), simplify = F, function(x){
    
    deseq.gene[[x]] %>% 
        rownames_to_column(var = "ensembl_gene_id") %>% 
        mutate(tissue = x)
    
}))

deseq.gene.merged <- merge(deseq.gene.merged,
                           transGrange(geneRanges),
                           by = 'ensembl_gene_id')

write.table(deseq.gene.merged, 
            file = paste0(table_dir, '02_deseq_results_gene.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

rm(deseq.gene.merged, deseq.gene)
