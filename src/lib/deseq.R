
# Condition data frame is created out of the header name. This is highly
# project specific. Peaks with less than 11 reads in sum across all 
# samples are removed. Possibly, you can set that threshold a bit higher.
# Check your notes for that, I recently read a paper where a threshold
# was written. The target attribute allows to select for which features DESeq2
# will be run - all (default): TEs & Genes together; te: te only; gene: gene only. 
doDEseq <- function(count.matrix, col_data, design_formula=formula(~condition), paral=FALSE, reference = NULL, target = 'all'){
    
    
    require(DESeq2)
    
    if (target != 'all') {
        count.matrix <- filterSalmonTEcounts(count.matrix, type = target)
    }
    
    # order of columns in count matrix needs to be equal to the order of the
    # col.data table.
    if (!all(rownames(col_data) %in% colnames(count.matrix))) {
        print("Check the file names in the col_data.csv and in your column names")
    }
    
    if (!all(rownames(col_data) == colnames(count.matrix))) {
        print("Change the order of the data file names")
        count.matrix <- count.matrix[, rownames(col_data)]
        
        if (all(rownames(col_data) == colnames(count.matrix))) {
            print("names are successfully rearranged")
    }}
    
    dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                                  colData = col_data,
                                  design = design_formula)
    
    
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    
    if (!is.null(reference)) {
        
        dds$condition <- relevel(dds$condition, ref = reference)
    }
    
    dds <- DESeq(dds, parallel = paral)
    
    return(dds)
}

# Extract the results from a dds object and returns the result table filtered
# for instances with an adjusted p-value. An lfc shrinkage is done by default
# but can turned of with lfcShrink = FALSE. The FDR filter removes all features
# without an adjusted p-value
getDEseqResults <- function(dds, lfcShrink = TRUE, coefficient = NULL, FDR.filter = TRUE){
    
    # Implement a trap if something is missing
    # if(lfcShrink & is.null(coefficient)){
    #     stop()
    # }
    
    require(DESeq2) 
    require(tidyverse)
    
    if (!(lfcShrink)) {
        
        deseq.res <- as.data.frame(results(dds)) %>% 
            filter(!is.na(padj)) 
        
        return(as.data.frame(deseq.res))
        
    } 
    
    deseq.res <- DESeq2::lfcShrink(dds, 
                                   coef = coefficient,
                                   parallel = F,
                                   res = DESeq2::results(dds),
                                   type = 'apeglm')
    
    if (FDR.filter) {
        deseq.res <- as.data.frame(deseq.res) %>% filter(!is.na(padj)) 
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

