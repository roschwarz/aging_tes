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

if (!exists("deseq.quant.merged")) {
    
    deseq.quant.merged <- read.csv(paste0(table_dir, '28_deseq_quant_all.csv'))
    
}

deseq.quant.merged$tissue <- factor(deseq.quant.merged$tissue, levels = c('brain', 'skin', 'blood'))

volcano.quant <- volcanoPlot(deseq.quant.merged, FDR = 0.05, "tissue") +
    theme_rob(12, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.quant <- color_strips(volcano.quant, 
                             bg_cols = tissue.color,
                             text_cols = tissue.text.color)

grid.draw(volcano.quant)


# ------------------------------ Gene ------------------------------------------

if (!exists("deseq.gene.quant.merged")) {
    
    deseq.gene.quant.merged <- read.csv(paste0(table_dir, '28_deseq_quant_gene.csv'))
    
}

deseq.gene.quant.merged$tissue <- factor(deseq.gene.quant.merged$tissue, levels = c('brain', 'skin', 'blood'))

volcano.gene.quant <- volcanoPlot(deseq.gene.quant.merged, FDR = 0.05, "tissue") +
    theme_rob(12, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.gene.quant <- color_strips(volcano.gene.quant, 
                                  bg_cols = tissue.color,
                                  text_cols = tissue.text.color)

grid.draw(volcano.gene.quant)

# ------------------------------- TEs ------------------------------------------

if (!exists("deseq.te.quant.merged")) {
    
    deseq.te.quant.merged <- read.csv(paste0(table_dir, '28_deseq_quant_te.csv'))
    
}

deseq.te.quant.merged$tissue <- factor(deseq.te.quant.merged$tissue,
                                      levels = c('brain', 'skin', 'blood'))

volcano.te.quant <- volcanoPlot(deseq.te.quant.merged, FDR = 0.05, "tissue") +
    theme_rob(12, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.te.quant <- color_strips(volcano.te.quant, 
                                bg_cols = tissue.color,
                                text_cols = tissue.text.color)


grid.draw(volcano.te.quant)
