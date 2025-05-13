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
                    "plot_kimura.R",
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
