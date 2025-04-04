#---------------- Own themes -----------

# from https://benjaminlouis-stat.fr/en/blog/2020-05-21-astuces-ggplot-rmarkdown/
theme_rob <- function(base_size = 14, base_family = 'arial'){
    
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            
            plot.title = element_text(size = rel(1), margin = margin(0,0,5,0), hjust = 0.5),
            
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            
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


dropBars <- function(pl){
    
    pl + scale_y_continuous(expand = expansion(mult = c(0, .07)))
    
}

# set the theme of interest
theme_set(theme_rob(12))

# color codes for specific scenarios
#tissue.color <- c(brain = '#264653', skin = '#2A9D8F',  blood = '#E9C46A')
tissue.color <- c(background = 'black',brain = '#58B2AA', skin = '#EFA081',  blood = '#685299')
tissue.text.color <- c( brain = "#ffffff", skin = "#000000", blood = "#ffffff")

#order.color = c(DNA = "#798E87", LINE = "#C27D38", LTR = "#CCC591", SINE = "#29211F")
order.color = c(DNA = "#E49400", LINE = "#831335", LTR = "#14273F", SINE = "#8493AE")
#direction.color = c(up = '#d1495b', down = '#66a182', equal = 'grey' )
direction.color = c(up = '#E63846', down = '#467B9D', equal = 'grey80' )