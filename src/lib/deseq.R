#' DESeq2 Utilities for Modular Analysis of Gene and TE Expression
#' Author: Robert Schwarz
#' Last updated: 2025-04-04
#' 

library(DESeq2)
library(tidyverse)

# Helper: Timestamped logging
logmsg <- function(msg) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", msg, "\n")
}


# Filter function placeholder (you can replace with your own)
filterCounts <- function(count.matrix, type) {
    
    if (type == "te") {
        count.matrix <- count.matrix %>% 
            rownames_to_column(var = 'id') %>% 
            filter(grepl('^chr', id)) %>% 
            column_to_rownames(var = 'id')
        return(count.matrix)
    } else if (type == "gene") {
        count.matrix <- count.matrix %>% 
            rownames_to_column(var = 'id') %>% 
            filter(grepl('^ENS', id)) %>% 
            column_to_rownames(var = 'id')
        return(count.matrix)
    } else {
        return(count.matrix)
    }
}


# Condition data frame is created out of the header name. This is highly
# project specific. Peaks with less than 11 reads in sum across all 
# samples are removed. Possibly, you can set that threshold a bit higher.
# Check your notes for that, I recently read a paper where a threshold
# was written. The target attribute allows to select for which features DESeq2
# will be run - all (default): TEs & Genes together; te: te only; gene: gene only. 
# Main DESeq2 Wrapper
doDEseq <- function(count.matrix, 
                    col_data, 
                    design_formula=formula(~condition), 
                    paral=FALSE, 
                    reference = NULL, 
                    target = 'all',
                    count_threshold = 10,
                    relevel_condition = "condition"){
    
    logmsg(paste("Starting DESeq2 for target:", target))
    
    if (target != 'all') {
        count.matrix <- filterCounts(count.matrix, type = target)
    }
    
    # order of columns in count matrix needs to be equal to the order of the
    # col.data table.
    if (!all(rownames(col_data) %in% colnames(count.matrix))) {
        stop("Not all sample names in col_data are found in count matrix column names.")
    }
        
    if (!all(rownames(col_data) == colnames(count.matrix))) {
        print("Reordering count matrix columns to match col_data row order.")
        count.matrix <- count.matrix[, rownames(col_data)]
        
        if (all(rownames(col_data) == colnames(count.matrix))) {
            print("Column names successfully matched.")
        }}
     
    dds <- DESeqDataSetFromMatrix(countData = round(count.matrix),
                                  colData = col_data,
                                  design = design_formula)
    
    logmsg(paste("Filtering features with total counts <", count_threshold))
    
    dds <- dds[rowSums(counts(dds)) >= count_threshold, ]
    
    if (!is.null(reference)) {
        logmsg(paste("Setting reference level to:", reference))
        #dds$condition <- relevel(dds$condition, ref = reference)
        dds[[relevel_condition]] <- relevel(dds[[relevel_condition]], ref = reference)
        
    }
    
    logmsg("Running DESeq...")
    dds <- DESeq(dds, parallel = paral)
    logmsg("DESeq complete.")
    
    return(dds)
}

# Extract the results from a dds object and returns the result table filtered
# for instances with an adjusted p-value. An lfc shrinkage is done by default
# but can turned of with lfcShrink = FALSE. The FDR filter removes all features
# without an adjusted p-value

# Extract and optionally shrink LFCs
getDEseqResults <- function(dds, 
                            lfcShrink = TRUE,
                            coefficient = NULL,
                            parallel = FALSE,
                            FDR.filter = TRUE){
    
    # Implement a trap if something is missing
    # if(lfcShrink & is.null(coefficient)){
    #     stop()
    # }
    
    if (is.null(coefficient)) {
        
        coefficient <- resultsNames(dds)[2]
        logmsg(paste("ℹ️ Automatically using coefficient:", coefficient))
    }
    
    
    
    if (!lfcShrink) {
        
        logmsg("Extracting raw DESeq2 results (no shrinkage)...")
        deseq.res <- results(dds, name = coefficient) 
    }
        
    
    if (lfcShrink) {
        
        logmsg("Applying LFC shrinkage with apeglm...")
        deseq.res <- DESeq2::lfcShrink(dds, 
                                       coef = coefficient,
                                       parallel = parallel,
                                       res = DESeq2::results(dds),
                                       type = 'apeglm')
    }    
    
    
    if (FDR.filter) {
        deseq.res <- as.data.frame(deseq.res) %>% 
            filter(!is.na(padj)) 
    }
    
    return(as.data.frame(deseq.res))
}




### DESeq specific plots

# PCA
# This function takes the dds object from deseq, the name of the column that was used for grouping the data, 
# the number of genes that will be considered, a named list of principle component comparisons that you want, 
# the title of the plot and an additional data frame if you want another label. The additional data frame for
# the labeling needs two columns: 
#  - name (usually the sample name [header of count table])
#  - new_label which contains the string that is added to the plot.
my_pca <- function(deseq_dds,
                   categories = 'condition',
                   ntop = 500,
                   pc_combinations = list(pca_1_2 = c(1,2),
                                          pca_1_3 = c(1,3),
                                          pca_1_4 = c(1,4),
                                          pca_2_3 = c(2,3),
                                          pca_2_4 = c(2,4),
                                          pca_3_4 = c(3,4)),
                   title = "Principal Component Analysis",
                   label.df = NULL
){
    
    vsd <- vst(deseq_dds, blind = FALSE)
    
    rv <- rowVars(assay(vsd))
    select <- order(rv, decreasing = T)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(vsd)[select,]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    
    if (!all(categories %in% names(colData(vsd)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    categories.df <-
        as.data.frame(colData(vsd)[, categories, drop = FALSE])
    
    if (length(categories) > 1) {
        group <- factor(apply(categories.df, 1, paste, collapse = " : "))
    } else {
        group <- colData(vsd)[[categories]]
    }
    
    sapply(pc_combinations, simplify = FALSE, function(pc.comb){
        
        if (length(pc.comb) != 2 | !is.numeric(pc.comb)) {
            stop(
                "the argument 'PCsToPlot' needs to be a numeric vector of exactly two integers like 'c(1,2)'."
            )
        }
        if (max(pc.comb) > dim(pca$x)[1] | min(pc.comb) < 1) {
            cat(
                paste0(
                    "at least one of the requested PCs is out of the dimension range which is 1:",
                    dim(pca$x)[1] ,
                    " for this data set."
                )
            )
            
            next
        }
        
        PCsToPlotText <- paste0("PC", pc.comb)
        
        d <- data.frame(
            PCx = pca$x[, pc.comb[1]],
            PCy = pca$x[, pc.comb[2]],
            group = group,
            categories.df,
            name = colnames(vsd)
        )
        
        if (!is.null(label.df)) {
            d <- merge(d, label.df, by = c('name'))
            names(d)[2] <- PCsToPlotText[1]
            names(d)[3] <- PCsToPlotText[2]
        }else{
            
            names(d)[1] <- PCsToPlotText[1]
            names(d)[2] <- PCsToPlotText[2]
            
        }
        
        # cat(paste0('creation of PC', pc.comb[1], " - PC ", pc.comb[2], ' plot\n'))
        pcaplot_name <-
            paste0("/PC", pc.comb[1], "vsPC", pc.comb[2], ".png")
        
        pl <-
            ggplot(data = d,
                   aes(x = !!sym(PCsToPlotText[1]),
                       y = !!sym(PCsToPlotText[2]),
                       color = group)) +
            geom_point(size = 3) +
            xlab(paste0(PCsToPlotText[1], ": ", round(percentVar[pc.comb[1]] *  100), "% variance")) +
            ylab(paste0(PCsToPlotText[2], ": ", round(percentVar[pc.comb[2]] * 100), "% variance"))  +
            ggtitle(title,
                    subtitle = paste(PCsToPlotText[1], PCsToPlotText[2], sep = " vs. ")) +
            theme_rob() +
            coord_fixed() +
            scale_color_manual(values = age_colors,
                               name = paste("Colored by:\n",
                                            gsub(", ",":",toString(categories)))) +
            geom_text_repel(
                data = d,
                aes(label = !!sym(case_when(
                    !is.null(label.df) ~ "new_label",
                    is.null(label.df) ~ "name"
                    ))),
                size = 3,
                color = "grey",
                direction = "both",
                #min.segment.length = 0.25,
                point.padding = 0.5
            )})
        
  
    
}

