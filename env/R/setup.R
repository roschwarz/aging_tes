###########
# general #
###########

# Helper: Timestamped logging
logmsg <- function(msg) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", msg, "\n")
}

# ------------------------------------------------------------------------------
# Load
# ------------------------------------------------------------------------------

load_project_fonts <- function() {
    message("Loading fonts...")
    extrafont::loadfonts()
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(DESeq2)
    library(blackRcloud)
    library(extrafont)
   
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

common_data = 'data/shared/'
table_dir = './results/tables/'
figure_dir = 'results/figures/'

create_project_dirs <- function() {
    dirs <- c("data/processed",
              "results/figures",
              "results/tables",
              "data/shared")
    for (d in dirs) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    }
}

# ------------------------------------------------------------------------------
# biological settings
# ------------------------------------------------------------------------------

tissues = c("brain", "skin", "blood")
orders.of.interest <- c("LINE", "SINE", "LTR", "DNA")
chr_of_interest <- paste0("chr", c(1:19, "X", "Y"))
FDR <- 0.05

# ------------------------------------------------------------------------------
# colors
# ------------------------------------------------------------------------------

# color codes for specific scenarios
tissue.color <- c(background = 'black',
                  brain = '#58B2AA',
                  skin = '#EFA081',
                  blood = '#685299')

tissue.text.color <- c( brain = "#ffffff",
                        skin = "#000000",
                        blood = "#ffffff")

order.color = c(DNA = "#E49400",
                LINE = "#831335",
                LTR = "#14273F",
                SINE = "#8493AE")

direction.color = c(up = '#E63846',
                    down = '#467B9D',
                    Up = '#E63846',
                    Down = '#467B9D',
                    equal = 'grey80',
                    "Not Sig" = 'grey80' )


# ------------------------------------------------------------------------------
# Analysis Stuff
# ------------------------------------------------------------------------------
#' @export
load_analysis_env <- function(){
    message("→ Loading analysis environment...")
    base_path <- paste0(getwd(), "/env/R/analysis")
    
    analysis_files <- c("deseq.R")
    
    for (f in analysis_files){
        source(file.path(base_path, f))
    }
    
    
}

# ------------------------------------------------------------------------------
# Annotations
# ------------------------------------------------------------------------------
#' @export
load_annotations <- function(){
    message("→ Loading annotations...")
    base_path <- paste0(getwd(), "/env/R/data")
    
    annotations <- c("annotations.R")
    
    for (f in annotations){
        source(file.path(base_path, f))
    }
    
    
}


# ------------------------------------------------------------------------------
# Plotting Stuff
# ------------------------------------------------------------------------------
#' @export
load_plotting_env <- function(){
    
    suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(ggalluvial)
    })
    
    message("→ Loading plotting settings...") 
    
    base_path <- paste0(getwd(), "/env/R/plots")
    
    plot_files <- c("theme_rob.R",
                    "plot_heatmaps.R",
                    "plot_volcano.R")
    
    for (f in plot_files){
        source(file.path(base_path, f))
    }
    
    # # set the theme of interest
    if (exists("theme_rob")) {
        message("Set standard theme to theme_rob")
        ggplot2::theme_set(theme_rob(12))
    } else {
        warning("theme_rob not found in plotting environment.")
    }
    
}

###################
# Module Specific #
###################

# ------------------------------------------------------------------------------
# RNA-Seq 
# ------------------------------------------------------------------------------

#' @export
load_rna_seq_env <- function(){
    
    message("→ Loading rna seq env...")
    
    # Count tables
    counts_rna <<- list(brain = "./results/rna_seq/brain/alignment_SalmonTE/EXPR.csv",
                       skin = "./results/rna_seq/skinII/alignment_SalmonTE/EXPR.csv",
                       blood = "./results/rna_seq/blood/alignment_SalmonTE/EXPR.csv")
    
    # Directories & Files
    rna_seq_results_dir <<- 'results/rna_seq/'
    rna_seq_deseq_dir <<- paste0(rna_seq_results_dir, 'deseq2/')
    
    deseq_dds_te <<- "dds_TE_instances_salmonTE.Rdata"
    deseq_results_te <<- "deseq_TE_instances_salmonTE.Rdata"
    deseq_results_te_csv <<- "02_deseq_results_te_instances.csv"
    
    deseq_dds_gene <<- "dds_genes_salmonTE.Rdata"
    deseq_results_gene <<- "deseq_genes_salmonTE.Rdata"
    deseq_results_gene_csv <<- "02_deseq_results_genes.csv"
    
    base_path <- paste0(getwd(), "/env/R/data")
    
    rna_files <- c("load_deseq_tes.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    
}


