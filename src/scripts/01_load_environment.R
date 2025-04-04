# Changelog
#
# 31012024
#
# Renamed mm10_TE_region_extended.bed to mm10_TE_region_5prime_extended.bed to clarify in
# which direction the TE region was extended. That is important because the Quant data is
# now integrated.


# ------------------------------- Packages -------------------------------------

setwd("../../")

library(blackRcloud) # own simple package
library(genomation) # package that allows to read bed files as a grange object
library(ggrepel)
library(DESeq2)
library(ggbreak)
library("ggh4x")
library(gtable)
library(pheatmap)
library(ggpubr)
library(ggpp)
library(ggrastr)
library(gridExtra)
library(ggvenn)
library(GO.db)
library(org.Mm.eg.db)
library(extrafont)
library(ggpmisc)
library(VennDiagram)
library(ggvenn)
library(ggalluvial) # for expressed te classes
library(cowplot)

# for heatmaps
library(circlize)
library(ComplexHeatmap)

# -------------------------- private functions ---------------------------------

source("src/lib/goAnalysis.R")
source("src/lib/plots.R")
source("src/lib/styles.R")
source("src/lib/deseq.R")
source("src/lib/dataHandling.R")
source("src/lib/helperFunctions.R")

# -------------------------- common settings -----------------------------------

chr_of_interest <- paste0('chr', c(seq(1,19), "X", "Y"))
tissues = c("brain", "skin", "blood")
orders.of.interest <- c("LINE", "SINE", "LTR", "DNA")
FDR = 0.05
extrafont::loadfonts()

volcano_colors <- setNames(c("#e63946", "gray", "#457b9d"), 
                           c("Up", "Not Sig", "Down")) 

data_dir = "data/processed/"
figure_dir = "results/figures/"
table_dir = "results/tables/"
common_data = "data/shared/"

if (!dir.exists(data_dir)) {
    dir.create(data_dir)
}

if (!dir.exists(figure_dir)) {
    dir.create(figure_dir)
}

if (!dir.exists(table_dir)) {
    dir.create(table_dir)
}
if (!dir.exists(common_data)) {
    dir.create(common_data)
}
# ------------------------------ common data -----------------------------------

## ============================= Count Tables ==================================

## #############################    RNA       ##################################     

counts_rna_RS <- list(brain = "results/rna_RS/salmonTE/brain_deduplicated/EXPR.csv",
                      skin = "results/rna_RS/salmonTE/skinII_deduplicated/EXPR.csv",
                      blood = "results/rna_RS/salmonTE/blood_deduplicated/EXPR.csv")

# counts_rna_RS <- list(brain = "results/rna_phd/salmonTE/brain_deduped/EXPR.csv",
#                       skin = "results/rna_phd/salmonTE/skinII_deduped/EXPR.csv",
#                       blood = "results/rna_phd/salmonTE/blood_deduped/EXPR.csv")

# count.files.te.regions <- list(brain = "../../rna/salmonTE/te.regions/brain_deduped/EXPR.csv",
#                          skin = "../../rna/salmonTE/te.regions/skinII_deduped/EXPR.csv",
#                          blood = "../../rna/salmonTE/te.regions/blood_deduped/EXPR.csv")

## #############################    CAGE       #################################     

# cage.counts <- "../../cage/rippchen/results_Gclipped/counted/"
# 
# cage.files <- list(
#     brain = paste0(cage.counts, "brain_down_sampled/brain.count.matrix.csv"),
#     blood = paste0(cage.counts, "blood.count.matrix.csv"),
#     skin = paste0(cage.counts,"skinII.count.matrix.csv")
# )

cage.files <- list(
    brain = 'results/cage_RS/brain_downsampled_gclipped/counts/brain_downsampled_peak_counts.csv',
    blood = 'results/cage_RS/blood_gclipped/counts/blood_peak_counts.csv',
    skin = 'results/cage_RS/skinII_gclipped/counts/skinII_peak_counts.csv')


## #############################    QUANT       ################################     

# cage.counts <- "../../cage/rippchen/results_Gclipped/counted/"
# 
# cage.files <- list(
#     brain = paste0(cage.counts, "brain_down_sampled/brain.count.matrix.csv"),
#     blood = paste0(cage.counts, "blood.count.matrix.csv"),
#     skin = paste0(cage.counts,"skinII.count.matrix.csv")
# )

quant.files <- list(
    brain = 'results/quant_RS/brain/counts/brain_peak_counts.csv',
    blood = 'results/quant_RS/blood/counts/blood_peak_counts.csv',
    skin = 'results/quant_RS/skinII/counts/skinII_peak_counts.csv')

# ------------------------------ annotations -----------------------------------

if (!exists("geneRanges")) {
    geneRanges <- blackRcloud::getGeneRanges()
    gene.exonRanges <- blackRcloud::getGeneRanges(type = "exon")
}


if (!exists("teRanges")) {
    

  #buildTEannotation(bedFile = "../../common_data/mm10/mm10_transposome.bed",
  #                  organism = "mmusculus",
  #                  gene.version = 102,
  #                  te.version = 1)
 
    teRanges <- getRanges("mmusculus", 102, "te", 1)
}


if (!file.exists(paste0(common_data, "te.annotation.Rdata"))) {
    
    te.annotation <- transGrange(teRanges) %>% 
        dplyr::rename(te_id = te.id) %>% 
        filter(!duplicated(te_id)) %>% 
        dplyr::select("te_id", "position", "ensembl_gene_id", "external_gene_name") %>% 
        splitTEID("te_id")
    
    
    save(te.annotation, file = paste0(common_data, "te.annotation.Rdata"))
    
}else{
    
    te.annotation <- loadRdata(paste0(common_data, "te.annotation.Rdata"))
    
}

# =============================== TE region ====================================
#
# Create and load TE region annotations. TEs within 500 bps are merged to TE
# regions, which are extended by 500 bp in either 5'-prime, 3'-prime, or both (full)
# directions. When the annotation files are not available than threy are created 
# with the respective shell script.


chrom.size <- paste0(common_data, "GRCm38.p6.chrom.sizes")
te.region.file <- paste0(common_data, "mm10_TE_region_5prime_extended.bed") # get rid of this file, however, you have to adapt the whole code base

# Updated version of TE region files for the integration of Quant data
te_region_file = 'data/shared/mm10_TE_region.bed'
te_region_5prime_file = 'data/shared/mm10_TE_region_5prime_extended.bed'
te_region_3prime_file = 'data/shared/mm10_TE_region_3prime_extended.bed'

te.annotation.file <- paste0(common_data, "mm10_transposome_sorted.bed")


if (!file.exists(te.region.file)) {
    
    system(paste("bash src/scripts/buildTE-region.annotation.sh",
                 te.annotation.file, 
                 chrom.size,
                 common_data), 
           wait = TRUE)
    
    teregionRanges <- readGeneric(te.region.file, 
                                  strand = 6, 
                                  meta.cols = list(names = 4))
    
    
    teRegionRanges <-
        readGeneric(te_region_file,
                    strand = 6,
                    meta.cols = list(names = 4))
    
    
    TEregionRanges_5prime <-
        readGeneric(te_region_5prime_file,
                    strand = 6,
                    meta.cols = list(names = 4)) 
    
    #te_region_3prime_file = '../../data/shared/mm10_TE_region_3prime_extended.bed'
    TEregionRanges_3prime <-
        readGeneric(te_region_3prime_file,
                    strand = 6,
                    meta.cols = list(names = 4)) 
    
}else{
    
    teregionRanges <- readGeneric(te.region.file,
                                  strand = 6,
                                  meta.cols = list(names = 4)) 
    
    teRegionRanges <-
        readGeneric(te_region_file,
                    strand = 6,
                    meta.cols = list(names = 4))
    
    TEregionRanges_5prime <-
        readGeneric(te_region_5prime_file,
                    strand = 6,
                    meta.cols = list(names = 4)) 
    
    TEregionRanges_3prime <-
        readGeneric(te_region_3prime_file,
                    strand = 6,
                    meta.cols = list(names = 4)) 
}

te.region.instances.file = paste0(common_data, "mm10_TE_region_instances.csv")

if (!file.exists(te.region.instances.file)) {
    
    TE.region.instances <- intersectGranger(teregionRanges, teRanges, tab = "all")

    # TEs can intersect with multiple genes, so that they occur multiple times in
    # the teRange objects. Therefore a unique is applied at the end of the pipeline
    TE.region.instances <- TE.region.instances %>% 
        dplyr::select(query.names, subject.te.id, subject.position) %>% 
        filter(!duplicated(subject.te.id)) %>% 
        unique()
    
    names(TE.region.instances) <- c("te_region_id", "te_id", "te_position")
    
    write.table(TE.region.instances, 
                file = te.region.instances.file, 
                col.names = T, 
                row.names = F,
                quote = F,
                sep = ",")
}

# ------------------------- TE-Gene-Relationships -------------------------------
#

te_region_gene_relation <- intersectGranger(teregionRanges, 
                                            geneRanges, 
                                            "all") %>% 
    dplyr::select("query.names", # TE region id 
                  "subject.ensembl_gene_id", 
                  "subject.external_gene_name")

names(te_region_gene_relation) <- c("te_region_id",
                                    "ensembl_gene_id",
                                    "external_gene_name")



ENVIRONMENT_LOADED <- TRUE


