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


