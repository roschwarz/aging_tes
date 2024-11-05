# --------------------------------- Notes --------------------------------------
# 
# This script creates Volcano Plots for the RNA-Seq data
#
# Output data (../figures):
#
# 03_volcano_te. - DESeq2 results for TEs based on SalmonTE counts
#

# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}
# -------- Volcano plot of TEs during aging in different tissues ---------------

if (!exists("deseq.te.merged")) {

    deseq.te.merged <- read.csv(paste0(table_dir, '02_deseq_results_te_instances.csv'))
    
}

# change order of tissues
deseq.te.merged$tissue <- factor(deseq.te.merged$tissue, 
                                 levels = c('brain', 'skin', 'blood'))

volcano.te <- volcanoPlot(cutPvalue(deseq.te.merged), FDR = 0.05, "tissue") +
    theme_rob(base_size = 10, base_family = 'arial') +
    theme(legend.position = 'None')


volcano.te <- color_strips(volcano.te, 
                           bg_cols = tissue.color[2:4], 
                           text_cols = c( "#ffffff", "#000000","#ffffff"))



ggsave(plot = volcano.te,
       filename = paste0(figure_dir, '03_te_instances_volcano_9.5x5.5_300.pdf'),
       device = cairo_pdf,
       width = 9.5,
       height = 5.5,
       units = "cm",
       dpi = 300
)

pl <- grid.draw(volcano.te)

show(pl)


#########
# Resis #
#########

ggsave(plot = volcano.te,
       filename = 'manuscripts/nature_aging/resis/1B_Volcano.pdf',
       device = cairo_pdf,
       width = 9.5,
       height = 5.5,
       units = "cm",
       dpi = 300
)


write.csv(deseq.te.merged, file = 'manuscripts/nature_aging/resis/1B_volcano.csv')

##### Next step should be to add a submission directory to your 01_load_environment.R file to store
### the stuff which is relevant for a submission.

ggsave(
    filename = paste0(figure_dir, '03_TE_Volcano_RNA.pdf'),
    plot = last_plot(),
    width = 20,
    height = 30,
    units = "cm",
    dpi = 300)

write.table(cutPvalue(deseq.te.merged), 
            file = '../../submission/Resis/figure1/figure_1b.csv', 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')


# -------- Single Volcanos for each tissue -------------------------------------


pl <- sapply(tissues, simplify = FALSE, function(x){
    
    df <- deseq.te.merged %>% 
        filter(tissue == x)

    volcanoPlot(cutPvalue(df), FDR = 0.05) +
        theme_rob(base_size = 20, base_family = 'Arial') +
        theme(legend.position = 'None')

    ggsave(
        filename = paste0(figure_dir, '03_', x, '_TE_Volcano_RNA.pdf'),
        plot = last_plot(),
        width = 20,
        height = 30,
        units = "cm",
        dpi = 300)
    
})


# -------- Volcano plot of gene during aging in different tissues --------------

if(!exists("deseq.gene.merged")){

    deseq.gene.merged <- read.csv(paste0(table_dir, '02_deseq_results_gene.csv'))
    
}

deseq.gene.merged$tissue <- factor(deseq.gene.merged$tissue, levels = c('brain', 'skin', 'blood'))

volcano.gene <- volcanoPlot(deseq.gene.merged, FDR = 0.05, "tissue") +
    theme_rob(10, base_family = 'arial') +
    theme(legend.position = 'None')

volcano.gene <- color_strips(volcano.gene, 
                             bg_cols = tissue.color[2:4], 
                             text_cols = c("#ffffff","#000000","#ffffff"))

grid.draw(volcano.gene)

