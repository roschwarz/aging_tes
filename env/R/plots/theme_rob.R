
# from https://benjaminlouis-stat.fr/en/blog/2020-05-21-astuces-ggplot-rmarkdown/
# @export
theme_rob <- function(base_size = 14, base_family = 'sans'){
    
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            
            plot.title = element_text(size = rel(1), margin = margin(0,0,5,0), hjust = 0.5),
            
            panel.grid.minor = element_blank(),
            
            # axes
            axis.title = element_text(size = rel(0.85), face = "bold"),
            axis.text = element_text(size = rel(0.70)),
            
            # Legend
            legend.title = element_text(size = rel(0.85), face = "bold"),
            legend.text = element_text(size = rel(0.70)),
            legend.background = element_rect(colour = 'lightgrey'),
            
            # Facett stuff
            panel.spacing.x = unit(0, "lines"),
            strip.background = element_rect(fill = 'white'),
            strip.text = element_text(size = rel(0.85), face = "bold")
        )
    
    
}
