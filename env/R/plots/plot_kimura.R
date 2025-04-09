#' ------------------------------------------------------------------------------
#' Histogram of Kimura distances colored by TE classes
#' ------------------------------------------------------------------------------
#' @param results A DESeq2 result object
#' @return A ggplot2 object
#' @export
kimHisto <- function(df, title = ""){
    ggplot(df, aes(Kimura, fill = order)) +
        geom_histogram(position = 'stack', color = 'black', binwidth = 1) +
        scale_fill_manual(values = order.color) +
        scale_x_continuous(breaks = seq(0,50,5)) +
        scale_y_continuous(expand = c(0,0)) +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = c(0.6,0.9),
              legend.direction = 'vertical',
              legend.background = element_blank(),
              legend.text = element_text(size = 7))
    
    
}
