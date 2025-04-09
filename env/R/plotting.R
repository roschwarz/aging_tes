#' @export
load_plots<- function(){
    message("→ Loading plotting settings...") 
    
    base_path <- paste0(getwd(), "/env/R/plots")
    # 
    plot_files <- c("theme_rob.R", "plot_volcano.R")
    # 
    # 
    for (f in plot_files){
        source(file.path(base_path, f))
    }
    # 
    # 
}

