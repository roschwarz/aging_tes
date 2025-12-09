# ------------------------------------------------------------------------------
# Project specific naming 
# ------------------------------------------------------------------------------

getConditions <- function(file_names) {
    # Project specific. It uses the naming convention of the fastq file to
    # determine the tissue and age of each sample.
    
    
    samples <- data.frame(
        SampleID = file_names,
        condition = file_names)
    
    # define the conditions for each sample
    if (grepl("skin", file_names[1], ignore.case = T)) {
        samples <- samples %>%
            separate(condition,
                     c("ID", "age"),
                     sep = "[_]") %>%
            mutate(age = str_replace(age, "[0-9]*Skin", "")) %>%
            mutate(condition = as.factor(case_when(
                age == "MO" ~ "old",
                age == "MY" ~ "young"
            ))) %>%
            dplyr::select(SampleID, condition) %>% 
            column_to_rownames("SampleID")
        
        return(samples)
    }
    
    samples <- samples %>%
        separate(condition,
                 c("ID", "org", "age", "tissue", "number"),
                 sep = "[_]") %>%
        mutate(condition = as.factor(case_when(age == "o" ~ "old",
                                               age == "y" ~ "young"))) %>%
        dplyr::select(SampleID, condition) %>% 
        column_to_rownames("SampleID")
    
    
    return(samples)
    
    
}

getCondition_cage <- function(header){
    
    # Returns a df with sample names (finally the row names of the df) and the
    # condition. The header needs to be separated by an underscore and the
    # condition needs to be at the second position. Usable for DESeq2.
    
    df <- data.frame(SampleID = header)
    df$condition <- sapply(df$SampleID, function(x) strsplit(x, "_")[[1]][2])
    
    rownames(df) <- df$SampleID
    df$SampleID <- NULL
    
    df$condition <- as.factor(df$condition)
    return(df)
}

updateHeader <- function(count.matrix){
    
    oldHeader = names(count.matrix)
    
    newHeaders = c()
    
    for(file in oldHeader){
        newHeader <- paste(getTissue(file), getAge(file), getId(file), sep="_")
        newHeaders = c(newHeaders, newHeader)
    }
    
    names(count.matrix) <- newHeaders
    return(count.matrix)
}



#======= Get functions ======

getTissue <- function(fileName){
    
    if(grepl("[bB]rain", fileName)){
        return('brain')
    }else if(grepl("[Bb]lood", fileName)){
        return('blood')
    }else if(grepl('_Skin_', fileName)){
        return('skinI')
    }else if(grepl('[0-9]Skin_', fileName)){
        return('skinII')
    }else{
        return(NA_character_)
    }
}


getAge <- function(fileName){
    
    if(grepl("_[oO]_|_MO|_Mo_", fileName)){
        return('old')
    }else if(grepl("_[yY]_|_MY|_My_", fileName)){
        return('young')
    }else{
        return(NA_character_)
    }
}

getId <- function(fileName){
    
    if(grepl("no[0-9]*", fileName)){
        
        for(part in strsplit(fileName, '[.]')[[1]]){
            if(grepl("^no[0-9]*", part)){
                return(part)
            }
        }
        
    }else if(grepl("Pool", fileName)){
        for(part in strsplit(fileName , '[_]')[[1]]){
            if(grepl("^[0-9]", part)){
                return(paste0('no00',part))
            }
        }
    }else{
        return(NA_character_)
    }
    
    
}



cutPvalue <- function(df, FDR_CAP = 0.00001){
    # This function allows to set the p-value to a certain max value if the adjuste
    # p-value is below. 
    
    df <- df %>% 
        mutate(padj = case_when(padj <= FDR_CAP ~ FDR_CAP,
                                TRUE ~ padj)) 
    
    return(df)
}

#================== Coloring  strips =============================
# To-Do, make it possible to submit the colors, then I guess you can
# add it to your blackRCloud
color_strips <- function(pl, 
                         bg_cols = c('#264653', '#2A9D8F', "#E9C46A", "#af7ac5", "#F2CC8F", "#3D405B", '#81B29A'),
                         text_cols = c('#ffffff', '#ffffff', "#000000", "#000000", "#000000", "#ffffff", '#000000')){
    
    g <- ggplot_gtable(ggplot_build(pl))
    
    strip_both <- which(grepl('strip-', g$layout$name))
    
    k <- 1
    
    if (length(strip_both) != length(bg_cols)) {
        print('Sorry the number of delivired colours is different compared to the number of facetts.')
        return(g)
    }
    
    for (i in seq_along(strip_both)) {
        
        j <- which(grepl('rect', g$grobs[[strip_both[i]]]$grobs[[1]]$childrenOrder))
        l <- which(grepl('titleGrob', g$grobs[[strip_both[i]]]$grobs[[1]]$childrenOrder))
        
        g$grobs[[strip_both[i]]]$grobs[[1]]$children[[j]]$gp$fill <- bg_cols[i]
        g$grobs[[strip_both[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- text_cols[i]
        k <- k + 1
    }
    
    return(g)
} 


# ------------------------------------------------------------------------------
# TE specific functions...
# ------------------------------------------------------------------------------

order.TEs <- function(df, col, decreasing = TRUE){
    # orders the data frame for a certain column by the occurrences of instances per group 
    # in the column. Is useful for geom_bar where you can order the data frame by decreasing
    # counts.
    # Example:
    # ggplot(order.TEs(TE.cluster.last.instances,'super_family'), aes(super_family)) +
    #     geom_bar() 
    
    df <- within(df,
                 ordered <- factor(df[[col]], levels = names(sort(table(df[[col]]),
                                                                  decreasing = decreasing ))))
    
    df[[col]] <- df$ordered
    
    return(df)
    
}


# ------------------------------------------------------------------------------
# featureCount loading...
# ------------------------------------------------------------------------------

getMeta.featureCounts <- function(featureCounts.table){
    
    require(tidyverse)
    return(featureCounts.table %>% dplyr::select(Chr, Start, End, Geneid, Length, Strand)) 
    
}

getLength.featureCounts <- function(featureCounts.table){
    require(tidyverse)
    return(featureCounts.table %>% dplyr::select(Geneid, Length)) 
    
}

getCountMatrix.featureCounts <- function(featureCounts.table){
    require(tidyverse)
    
    countMatrix <- column_to_rownames(featureCounts.table, var = 'Geneid') %>% 
        dplyr::select(!(Chr:Length))
    
    return(countMatrix)
    
}

getFeatureCountTab <- function(file, filter = TRUE, threshold = 10){
    
    featureCountTable <- readFeatureCountsMat(file)
    
    length <- getLength.featureCounts(featureCountTable)
    
    meta <- getMeta.featureCounts(featureCountTable)
    
    count.matrix <- getCountMatrix.featureCounts(featureCountTable)
    
    if(filter){
        
        count.matrix <- count.matrix[,grepl('unique.sorted', names(count.matrix))]
        count.matrix <- count.matrix[ rowSums(count.matrix) > threshold, ]
        count.tpm <- normalizeCountMatrix(count.matrix, length)
        
    }else{
        
        count.tpm <- normalizeCountMatrix(count.matrix, length)
        
    }
    
    return(list(counts = count.matrix, tpms = count.tpm, meta = meta))  
    
}

getCountTable <- function(file, filter = TRUE, threshold = 10){
    
    featureCountTable <- readFeatureCountsMat(file)
    
    length <- getLength.featureCounts(featureCountTable)
    
    count.matrix <- getCountMatrix.featureCounts(featureCountTable)
    
    if(filter){
        
        count.matrix <- count.matrix[,grepl('unique.sorted', names(count.matrix))]
        count.matrix <- count.matrix[ rowSums(count.matrix) > threshold, ]
        count.tpm <- normalizeCountMatrix(count.matrix, length)
        
    }else{
        
        count.tpm <- normalizeCountMatrix(count.matrix, length)
        
    }
    
    return(list(counts = count.matrix, tpms = count.tpm))  
}


# Skip the first row, as the featureCount call is stored there
readFeatureCountsMat <- function(file){
    read.csv(file, skip =1, sep = "\t", header = T)
}

getCoordFeatureCountsMat <- function(file){
    
    df <- readFeatureCountsMat(file)
    
    return(df[1:5])
}


ask_yes_no <- function(prompt, default = "n") {
    norm <- function(x) tolower(trimws(x))
    def <- if (is.null(default)) NULL else norm(default)
    repeat {
        ans <- if (interactive()) {
            readline(prompt)
        } else {
            cat(prompt)
            line <- readLines(con = stdin(), n = 1, warn = FALSE)
            if (length(line) == 0) "" else line
        }
        ans <- norm(ans)
        if (ans == "" && !is.null(def)) return(def %in% c("y", "yes", "j", "ja"))
        if (ans %in% c("y", "yes", "j", "ja")) return(TRUE)
        if (ans %in% c("n", "no", "nein")) return(FALSE)
        cat("Please provide 'y' or 'n'.\n")
    }
}

# ------------------------------------------------------------------------------
# Save figures with an index table
# ------------------------------------------------------------------------------

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
        response <- ask_yes_no("Do you want to overwrite it? (y/n):") #readline("Do you want to overwrite it? (y/n): ")
        if (!response) {
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
    } else if (inherits(plot, "HeatmapList") || inherits(plot, "Heatmap")) {
        
        # convert cm into inch
        width <- width / 2.54
        height <- height / 2.54
        
        pdf(file = image_path, width = width, height = height)
            draw(plot, heatmap_legend_side = "bottom")
        dev.off()
    } else if (inherits(plot, "upset")) {
        
        # convert cm into inch
        width <- width / 2.54
        height <- height / 2.54
        
        pdf(file = image_path, width = width, height = height)
            show(plot)
        dev.off()
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


#' Generic function to find overlaps between two sets of genomic ranges
#' and summarize the results.
#' @param query_ranges GRanges object representing the query features (e.g., gene TSSs)
#' @param subject_ranges GRanges object representing the subject features (e.g., TEs)
#' @param id_col Optional string specifying the column name in mcols(query_ranges) that uniquely identifies the entity (e.g., gene ID). If NULL, all distinct query features are considered.
#' @param extend_bp Integer specifying the number of base pairs to extend the query ranges on both sides (default is 0, meaning no extension).
#' @param verbose Logical indicating whether to print summary statistics (default is TRUE).
#' @return A list containing:
#'   - query_with_hit: GRanges of query features that have overlaps
#'   - subject_hits: GRanges of subject features that overlap with the query
#'   - stats: A list with total number of query features, number of unique entities with hits, and percentage with hits.
#' @examples
#' # Example usage:
#' results <- process_overlapping(query_ranges = gene_tss_gr, subject_ranges = teRanges, id_col = "ensembl_gene_id", extend_bp = 200)
#' @note
#' # This function requires the GenomicRanges package.
#' 
#' subject_ranges <- GRanges(seqnames = rep("chr1", 3),
#'                           ranges = IRanges(start = c(1,11,21), end = c(10,20, 30)),
#'                           strand = c("+", "+", "+"),
#'                           gene_id = paste0("gene", 1:3),
#'                           score = 1:3)
#' query_ranges <- GRanges(seqnames = rep("chr1", 3),
#'                         ranges = IRanges(1:3, width=1),
#'                         strand = c("+", "+", "+"),
#'                         gene_id = paste0("gene", 1:3),
#'                         score = 1:3)
#'
#' process_overlapping(query_ranges, subject_ranges) 
#' This returns that 100 % of the query features overlap with the subject features. Because all three TSS (
#' query_ranges) overlap with the subject ranges.
#' 
#' #' process_overlapping(subject_ranges, query_ranges)
#' This returns that 33.33 % of the query features overlap with the subject features. Because only one of the 
#' three subject ranges (subject_ranges) overlap with the query ranges. 
#' @export
process_overlapping <- function(
        query_ranges,         # Ranges to search for overlaps (e.g., gene TSSs)
        subject_ranges,       # Ranges to overlap against (e.g., TEs)
        id_col = NULL,        # Name of the column in mcols(query_ranges) that uniquely identifies the entity (optional)
        extend_bp = 0,        # Window size to extend query_ranges (in bp, optional)
        col_query_to_subj = NULL,  # Metadata columns to transfer from subject (character vector)
        verbose = TRUE        # If TRUE, prints results
) {
    # Optionally extend query ranges
    if (extend_bp > 0) {
        query_ranges <- resize(query_ranges, width = extend_bp * 2 + 1, fix = "center")
    }
    # Perform overlap
    hits <- findOverlaps(query_ranges, subject_ranges)
    
    if (length(hits) == 0) {
        warning("No overlaps found.")
        return(NULL)
    }
    
    
    overlapped_query <- query_ranges[queryHits(hits)]
    overlapped_subject <- subject_ranges[subjectHits(hits)]
    
    if (!is.null(col_query_to_subj)) {
        mcols(overlapped_subject)[col_query_to_subj] <- mcols(overlapped_query)[col_query_to_subj]
    }
    # Identifier logic: Use id_col if provided
    if (!is.null(id_col)) {
        # Check existence
        if (!(id_col %in% names(mcols(query_ranges)))) {
            stop(sprintf("Identifier column '%s' not found in query_ranges.", id_col))
        }
        unique_entities <- unique(mcols(overlapped_query)[[id_col]])
    } else {
        # Fall back to all distinct queries
        unique_entities <- unique(queryHits(hits))
    }
    n_query_total <- length(query_ranges)
    n_entities_with_hit <- length(unique_entities)
    pct_with_hit <- 100 * n_entities_with_hit / n_query_total
    
    # Print summary if verbose
    if (verbose) {
        cat(sprintf("Total query features: %d\n", n_query_total))
        cat(sprintf("Unique query entities overlapping: %d (%.2f%%)\n",
                    n_entities_with_hit, pct_with_hit))
        cat(sprintf("Window size: ±%d bp\n", extend_bp))
    }
    
    # Output
    results <- list(
        query_with_hit = overlapped_query,
        subject_hits = overlapped_subject,
        stats = list(
            n_query_total = n_query_total,
            n_entities_with_hit = n_entities_with_hit,
            pct_with_hit = pct_with_hit
        )
    )
    return(results)
}


bootStrapping <- function(n.intronic.tes.target,
                          n.intronic.tes.background,
                          moves = 1000){
    
    random.target <- sample_n(n.intronic.tes.target,
                              moves,
                              replace = T)[['n.intronicTEs']] 
    
    random.background <- sample_n(n.intronic.tes.background,
                                  moves,
                                  replace = T)[['n.intronicTEs']]
    
    
    r <- (random.target + 0.1)/(random.background + 0.1)
    
    
    result <- data.frame(n.target = random.target, 
                         n.background = random.background, 
                         ratio = r)
    
    return(result)
}
