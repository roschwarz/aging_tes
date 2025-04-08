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
table_dir = 'results/tables/'
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
#tissue.color <- c(brain = '#264653', skin = '#2A9D8F',  blood = '#E9C46A')
tissue.color <- c(background = 'black',brain = '#58B2AA', skin = '#EFA081',  blood = '#685299')
tissue.text.color <- c( brain = "#ffffff", skin = "#000000", blood = "#ffffff")

#order.color = c(DNA = "#798E87", LINE = "#C27D38", LTR = "#CCC591", SINE = "#29211F")
order.color = c(DNA = "#E49400", LINE = "#831335", LTR = "#14273F", SINE = "#8493AE")
#direction.color = c(up = '#d1495b', down = '#66a182', equal = 'grey' )
direction.color = c(up = '#E63846', down = '#467B9D', Up = '#E63846', Down = '#467B9D', equal = 'grey80', "Not Sig" = 'grey80' )

###################
# Module Specific #
###################

# ------------------------------------------------------------------------------
# RNA-Seq 
# ------------------------------------------------------------------------------

counts_rna <- list(brain = "./results/rna_seq/brain/alignment_SalmonTE/EXPR.csv",
                   skin = "./results/rna_seq/skinII/alignment_SalmonTE/EXPR.csv",
                   blood = "./results/rna_seq/blood/alignment_SalmonTE/EXPR.csv")

rna_seq_results_dir = 'results/rna_seq/'
rna_seq_deseq_dir = paste0(rna_seq_results_dir, 'deseq2/')

deseq_dds_te = "dds_TE_instances_salmonTE.Rdata"
deseq_results_te = "deseq_TE_instances_salmonTE.Rdata"

deseq_dds_gene = "dds_genes_salmonTE.Rdata"
deseq_results_gene = "deseq_genes_salmonTE.Rdata"



