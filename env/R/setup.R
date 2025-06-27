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
    library(data.table)
    library(genomation)
    library(Gviz)
    library(GenomicRanges)
    library(plyranges)
    library(BRGenomics)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
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
              "results/te_island",
              "data/shared")
    for (d in dirs) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    }
}

# ------------------------------------------------------------------------------
# biological settings
# ------------------------------------------------------------------------------

tissues = c("brain", "skin", "blood")
female_tissues = c("brain", "skin")
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
                  blood = '#685299',
                  liver = '#adc178',
                  'Gastrocnemius muscle' = '#a98467',
                  'White adipose tissue' = '#f0ead2',
                  male = "red")

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












