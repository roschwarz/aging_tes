# --------------------------------- Notes --------------------------------------
# 
# This script is used for the CAGE-Seq quantification
#
# Output data (../data/cage; ../tables):
#

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    
    source('./01_load_environment.R')
    
} else if ( !ENVIRONMENT_LOADED ) {
    
    source('./01_load_environment.R')
    
}


# ------------------------------- All ------------------------------------------

if (!exists("deseq.cage.merged")) {
    
    deseq.cage.merged <- read.csv(paste0(table_dir, '07_deseq_cage_all.csv'))
    
}

deseq.cage.merged$tissue <- factor(deseq.cage.merged$tissue, levels = c('brain', 'skin', 'blood'))

volcano.cage <- volcanoPlot(deseq.cage.merged, FDR = 0.05, "tissue") +
    theme_rob(12, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.cage <- color_strips(volcano.cage, 
                             bg_cols = tissue.color,
                             text_cols = tissue.text.color)

grid.draw(volcano.cage)

# ------------------------------ Gene ------------------------------------------

if (!exists("deseq.gene.cage.merged")) {
    
    deseq.gene.cage.merged <- read.csv(paste0(table_dir, '07_deseq_cage_gene.csv'))
    
}

deseq.gene.cage.merged$tissue <- factor(deseq.gene.cage.merged$tissue, levels = c('brain', 'skin', 'blood'))

volcano.gene.cage <- volcanoPlot(deseq.gene.cage.merged, FDR = 0.05, "tissue") +
    theme_rob(12, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.gene.cage <- color_strips(volcano.gene.cage, 
                                  bg_cols = tissue.color,
                                  text_cols = tissue.text.color)

grid.draw(volcano.gene.cage)

# ------------------------------- TEs ------------------------------------------

if (!exists("deseq.te.cage.merged")) {
    
    deseq.te.cage.merged <- read.csv(paste0(table_dir, '07_deseq_cage_te.csv'))
    
}

deseq.te.cage.merged$tissue <- factor(deseq.te.cage.merged$tissue,
                                      levels = c('brain', 'skin', 'blood'))

volcano.te.cage <- volcanoPlot(deseq.te.cage.merged, FDR = 0.05, "tissue") +
    theme_rob(12, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.te.cage <- color_strips(volcano.te.cage, 
                                bg_cols = tissue.color,
                                text_cols = tissue.text.color)


grid.draw(volcano.te.cage)
