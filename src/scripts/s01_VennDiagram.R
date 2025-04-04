
# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}

deseq.te.merged <- read.csv(paste0(table_dir, '02_deseq_results_te_instances.csv'))


expressed.TEs <- sapply(c('brain', 'skin', 'blood'), function(x){
    
    deseq.te.merged %>% 
        filter(tissue == x, !is.na(padj)) %>% 
        pull(te_id) %>% 
        unique()
})


venn <- venn.diagram(
    x = expressed.TEs,
    category.names = names(expressed.TEs),
    
    
    # Circles
    lwd = 2,  
    fill = tissue.color[2:4], #c('#264653', '#2A9D8F',  '#E9C46A'),
    alpha = c(0.7, 0.7, 0.7),
    
    
    # Number
    cex = 1, # font size
    fontface = "bold",
    fontfamily = "arial",
    # 
    # # Set names
    cat.cex = 1.5,
    cat.default.pos = "outer",
    cat.fontface = "bold",
    cat.fontfamily = "arial",
    # cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055, 0.055),
    # main = header,
    scaled = F,
    print.mode = c("raw", "percent"),
    
    # Output
    filename = NULL, #paste0(figures, 'Panel_2C_VennDiagram.svg'),
    imagetype = "svg",
    output = FALSE,
    width = 200,
    height = 500,
    resolution = 300,
    disable.logging = TRUE
    
)


grid::grid.draw(venn)



ggvenn(expressed.TEs, 
       stroke_size = 0.5, 
       set_name_size = 12/.pt, 
       digits = 1,
       label_sep = ",",
       text_size = 8/.pt,
       fill_color = unname(tissue.color[2:4]))


### the stuff which is relevant for a submission.

ggsave(
    filename = paste0(figure_dir, 's01_vennDiagram.pdf'),
    plot = last_plot(),
    width = 10,
    height = 10,
    units = "cm",
    dpi = 300)

ggsave(
    filename = paste0(figure_dir, 's01_vennDiagram.png'),
    plot = last_plot(),
    width = 10,
    height = 10,
    units = "cm",
    dpi = 300)

#########
# RESIS #
#########


ggsave(
    filename = './manuscripts/nature_aging/resis/figures/supplemental_figure_1/S1_vennDiagram.png',
    plot = last_plot(),
    width = 10,
    height = 10,
    units = "cm",
    dpi = 300)
