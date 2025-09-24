


fig_index <- function(plot, outdir, meta, index_file = 'index.tsv', width = 6, height = 4, dpi = 300, format = 'pdf') {
    
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
    
  name <- if (!is.null(meta$name)) meta$name else paste0("unamed_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  filename <- paste0(name, ".", format)
  image_path <- file.path(outdir, filename)
  
  # Check if the figure already exists
  file_exists <- file.exists(image_path)
  if (file_exists) {
      cat("WARNING: File", image_path, "allready exists.\n")
      response <- readline("Do you want to overwrite it? (y/n): ")
      if (tolower(response) != "y") {
          cat("Stop: File will not be overwritten.\n")
          return(invisible(NULL))
      } else {
          cat("Overwriting file...\n")
      }
  }
  
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(plot = plot, 
                    filename = image_path,
                    device = cairo_pdf,
                    units = "cm",
                    dpi = dpi,
                    width = width,
                    height = height)
  } else if (inherits(plot, "grob")) {
      ggplot2::ggsave(plot = plot, 
                      filename = image_path,
                      device = cairo_pdf,
                      units = "cm",
                      dpi = dpi,
                      width = width,
                      height = height)
  } else {
    stop("Unsupported plot type. Please provide a ggplot or grob object.")
  }
  
  
  # fill meta object with standard fields
  meta_filled <- list(
      filename = filename,
      date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      description = if (!is.null(meta$description)) meta$description else "",
      tags = if (!is.null(meta$tags)) paste(meta$tags, collapse = ",") else "",
      parameters = if (!is.null(meta$parameters)) toString(meta$parameters) else ""
  )
  # take any additional fields from meta
  for (k in names(meta)) {
      if (is.null(meta_filled[[k]])) meta_filled[[k]] <- meta[[k]]
  }
  
  
  index_path <- file.path(outdir, index_file)
  
  if (file.exists(index_path)) {
    index_data <- read.table(index_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Check for duplicate names
    existing_row <- which(index_data$filename == filename)
    if (length(existing_row) > 0) {
        index_data[existing_row, ] <- as.data.frame(meta_filled, stringsAsFactors = FALSE)
    }else{
        index_data <- rbind(index_data, as.data.frame(meta_filled, stringsAsFactors = FALSE))
    }
  } else {
    index_data <- as.data.frame(meta_filled, stringAsFactors = FALSE)
}
  
  write.table(index_data, file = index_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Figure saved to:", image_path, "\n")
  cat("Index updated at:", index_path, "\n")
  
  invisible(image_path)
  
  
}
