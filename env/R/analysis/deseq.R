#' DESeq2 Utilities for Modular Analysis of Gene and TE Expression
#' Author: Robert Schwarz
#' Last updated: 2025-04-04

library(DESeq2)
library(tidyverse)

# Filter function
filterCounts <- function(count.matrix, type) {
    
    if (type == "te") {
        
        logmsg("Filter count matrix for transposable elements (ID begins with chr)")
        count.matrix <- count.matrix %>% 
            rownames_to_column(var = 'id') %>% 
            filter(grepl('^chr', id)) %>% 
            column_to_rownames(var = 'id')
        return(count.matrix)
    } else if (type == "gene") {
        logmsg("Filter count matrix for genes (ID begins with ENS)")
        count.matrix <- count.matrix %>% 
            rownames_to_column(var = 'id') %>% 
            filter(grepl('^ENS', id)) %>% 
            column_to_rownames(var = 'id')
        return(count.matrix)
    } else {
        logmsg(paste("Count matrix can not be filtered for", type))
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
                    count_threshold = 10){
    
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
        dds$condition <- relevel(dds$condition, ref = reference)
    }
    
    logmsg("Running DESeq...")
    dds <- DESeq(dds, parallel = paral)
    logmsg("✅ DESeq complete.")
    
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
        logmsg(paste("Automatically using coefficient:", coefficient))
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
