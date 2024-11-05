library(data.table)

#======= load functions =================

readInputCage <- function(path, name){
    
    age=strsplit(name, "_")[[1]][2]
    tissue=strsplit(name, '_')[[1]][1]
    strand=strsplit(name, '[._]')[[1]][3]
    
    df <- fread(paste0(path, name), sep = '\t') 
    
    df$age <- age
    df$tissue <- tissue
    df$m.strand <- strand
    
    return(df)
}

readIntersection <- function(path, name, feature='gene', type='peak'){
    
    if(!tools::file_ext(name) == 'bed'){
        print(name)
        return(NULL) 
    }
    
    
    tissue=strsplit(name, '_')[[1]][1]
    age=strsplit(name, "_")[[1]][2]
    strand=strsplit(name, '[._]')[[1]][3]
    genomicFeature=strsplit(name, '[.]')[[1]][3]     
    intersection.type=strsplit(name, '[.]')[[1]][2]
    
    if( type != intersection.type | feature != tolower(genomicFeature)){
        
        return(NULL)
        
    }else{
        
        df <- as.data.frame(fread(paste0(path, name), sep = '\t'))
        
        if(tolower(genomicFeature) == 'te'){

            df['V5'] <- NULL
            names(df) <- c('V1', 'V2', 'V3', 'V4', 'V5')

        }

        df$age <- age
        df$tissue <- tissue
        df$m.strand <- strand
        df$feature <- genomicFeature

        return(df) 
    }
    
}

readInputQuant <- function(path, name){
    
    if(tools::file_ext(name) == 'bed'){
        
        age=strsplit(name, "_")[[1]][2]
        tissue=strsplit(name, '[.]')[[1]][1]
        strand=strsplit(name, '[._]')[[1]][3]
        feature=gsub('CAGEintersect', '', strsplit(name, '[._]')[[1]][4])
        
        df <- as.data.frame(fread(paste0(path, name), sep = '\t')) 
        
        if(feature == 'TE'){
            
            df['V5'] <- NULL
            names(df) <- c('V1', 'V2', 'V3', 'V4', 'V5')
            
        }
        
        
        df$age <- age
        df$tissue <- tissue
        df$m.strand <- strand
        df$feature <- feature
        
        return(df)
        
    }else{
        return(NULL)
    }
}

# Skip the first row, as the featureCount call is stored there
readFeatureCountsMat <- function(file){
    read.csv(file, skip =1, sep = "\t", header = T)
}

getCoordFeatureCountsMat <- function(file){
    
    df <- readFeatureCountsMat(file)
    
    return(df[1:5])
}


readRippchenStats <- function(stat.dir){
    
    files <- list.files(stat.dir, 
                        pattern = '.stat', 
                        full.names = T)
    
    df <- do.call(rbind, lapply(files, function(x) read.csv(x, sep = '\t', header = F)))
    
    names(df) <-  c('sample', 'type', 'count')
    
    df$type <- ifelse(
        df$type == 'polyntclipped reads', 'rawCounts', 
        ifelse(
            df$type == 'mapped reads', 'mappedCounts', 'uniqueCounts'
        ) 
    )
    
    df <- spread(df, type, count)
    
    
    return(df)
}

# simplify is needed to get the df per tissue.
loadPeakAnnotations <- function(list.file.names){
    
    peak.annotation <- sapply(names(list.file.names), simplify = F, function(tissue){
        df <- read.csv(list.file.names[[tissue]], sep = '\t', header = F, stringsAsFactors = F)
        names(df) <- c('peakName', 'geneID')
        df} 
    )
    
    return(peak.annotation)
}



loadClosestFile <- function(file){
    # <- .bed-file (output of bedtools closest with -s -d flag) 
    # -> data.frame
    # There are double entries in the closest file. I assume that are double
    # entries in the align file from repeatMasker that I used to parse the
    # TE annotations into a bed format.
    df <- fread(file, header = F, sep = '\t')
    
    names(df) <- c('chrom', 'te.start', 'te.end', 'TE.ID', 'te.score', 'te.strand', 'gene.chrom', 'gene.start', 'gene.end', 'gene.ID', 'gene.score', 'gene.strand', 'distance')
    
    df <- df %>% filter(!duplicated(TE.ID)) 
    df[,c('te.score', 'gene.chrom', 'gene.score', 'gene.strand')] <- NULL
    
    return(as.data.frame(df))
}


#======= data storage ======

saveDESeqCSV <- function(df, resultDir='.', fileName='deseq.result.csv'){

    if(any(is.na(df$padj))){

        print('Features without an adjusted p-value are removed from the table')

        df <- df %>% filter(!is.na(padj))

    }

    if(names(df)[1] == 'baseMean'){
        df <- df %>% rownames_to_column(var = 'id')    
    }

    if(!grepl('.csv$', fileName)){
            fileName = paste0(fileName, '.csv')
    }

    write.csv(df, paste0(resultDir, fileName))


}

#======= table manipulation ======

prepData <- function(df){
    
    if(!grepl('^ENS', df$V4[1])){
        df <- splitTEID(df, 'V4')
    }
    
    df$identifier <- paste(df$tissue, df$age, sep = '.')
    
    df <- df %>% mutate(age = case_when(age == 'o' ~ 'old (24 mo)',
                                        age == 'y' ~ 'young (6 mo)',
                                        age == 'young' ~ 'young (6 mo)',
                                        age == 'old' ~ 'old (24 mo)',
                                        ),
                        tissue = case_when(tissue == 'SkinI' ~ 'skin_I',
                                           tissue == 'SkinII' ~ 'skin_II',
                                           tissue =='Brain' ~ 'brain',
                                           tissue == 'Blood' ~ 'blood',
                                           TRUE ~ tissue))
    
    df$age <- factor(df$age, levels = c('young (6 mo)', 'old (24 mo)'))
    
    df <- df %>% group_by(tissue, age) %>% filter(!duplicated(V4))
    
}


addTissue <- function(df){
  # extract tissues from the sample name and store them in a new column
  df$tissue <- ifelse(grepl('rain', df$sample), 'brain', 
                                   ifelse(grepl('lood', df$sample), 'blood', 
                                          ifelse(grepl('_Skin', df$sample), 'skin_I', 'skin_II')))

  return(df)
}

addAge <- function(df){
  # extract age from the sample name and store them in a new column
  df$age_group <- ifelse(grepl('Mo|MO|_o_', df$sample), 'old', 'young')
  return(df)
}


splitGeneID <- function(df){
    # SalmonTE table contains an gene identifier that is combined by ensembl_gene_id
    # external gene name and transcript id
    df %>% 
        rownames_to_column(var='ID') %>% 
        tidyr::separate(ID, into = c('trans', 'version', 'gene','ensembl_gene_id', 
                                     'gene.version', 'gen_symbol', 'external_gene_name'), sep = '[.|:]') %>% 
        dplyr::mutate(ensembl_transcript_id = paste(trans, version, sep = '.')) %>% 
        dplyr::select(-c('trans', 'version', 'gene.version', 'gen_symbol', 'gene')) 
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


getDiffTes <- function(df, fdr = 0.05){
    
    if(class(df) != 'data.frame'){
        df <- as.data.frame(df)
    }
   
    if(!("TE.ID" %in% names(df))){
        df <- df %>% rownames_to_column("TE.ID")
    } 
    
    df <- df %>% 
        filter(padj <= fdr) %>% 
        dplyr::select(TE.ID, log2FoldChange, baseMean, padj)
    
    return(df)
}


collectExpFeatureSet <- function(df){
    # Returns all elements of a DESeq table that are considered as expressed.
    # An element is expressed when the adjusted p-value is not NA
    df <- df[!is.na(df['padj']),]
    
    return(rownames(df))
    
}

# The analysis of a GO-enrichment analysis returns a result table where all
# genes of the target gene set are contained that are associated with the 
# GO-terms of interest. This function extracts these genes and select the TEs
# that are associated with those genes.
#
# te.associated.genes contains all expressed TEs that are associated with 
# expressed genes (see genTEassociation.R).

getGoTEs <- function(go.result, te.associated.genes){
    
    if(nrow(go.result) == 0){
        return(NULL)
    }
    
    gene.symbols <- getGoSymbols(go.result)
    
    te <- te.associated.genes %>% 
        dplyr::select(external_gene_name, TEs) %>% 
        filter(external_gene_name %in% gene.symbols) %>% 
        pull(TEs)
    
    te <- unique(do.call('c', strsplit(te, ';')))
    
    te <- data.frame(TE.ID = te)
    
    te <- splitTEID(te, "TE.ID")
    
    return(te)
}


getGoSymbols <- function(go.result){
    
    if(nrow(go.result) == 0){
        return(NULL)
    }
    gene.symbols <- go.result %>% dplyr::pull(geneID)
    
    
    # split the strings (strsplit) and add them to one vector (do.call) and extract
    # only each gene one time (unique)
    gene.symbols <- unique(do.call('c', strsplit(gene.symbols, '/')))
    
    return(gene.symbols)
}


getGOgenCounts <- function(go.result, con = 'NA'){
    
    sym <- getGoSymbols(go.result)
    
    count <- length(sym)
    
    if(count == 0){
        return(NULL)
    }
    
    df <- data.frame(tissue = con,
                     n = count,
                     label = paste('n of TE asso DEGs:', count)
    )
    return(df)
}
    

#This function takes a data frame with a TE.ID column and split these IDs
# to create a data frame that contains all field of a .bed-file.
getTEbed <- function(df){
    
    if(any(duplicated(df[['TE.ID']]))){
        warning('There are double entries in your data frame. Double entries are removed')
        df <- df %>% dplyr::filter(!duplicated(TE.ID)) 
    }
    
    df <- df %>% splitTEID('TE.ID')
    
    df %>%
        dplyr::select('chromosome', 'start', 'end', 'TE.ID', 'order', 'strand')
    
    
}

mergeTissueRes <- function(seq.object, target = 'te.dres'){
# reads the DESeq result tables and merge the different tissues
# seq.object are objects that are created in my quantification.Rmd and contains
# a lot of tables. You can apply to tables by the target and it will merge the
# tables of the different tissues by adding a tissue column and merging of the
# data frames.
    
   df <- do.call('rbind', sapply(names(seq.object[[target]]), simplify = F, function(x){
    
        df <- seq.object[[target]][[x]]
    
        df <- df %>% filter(!duplicated(TE.ID))
     if(!is.null(df)){
        
            df$tissue <- x
     }
    
        return(df)
    }
    )) 
   
   return(df)
} 
#======== Write functions ==================
# This function expects a data frame that will be stored at the hard drive.
# To-DO:
#   - check for the correct format
writeBedFile <- function(bed.df, file.name=NULL, directory = NULL){
    
    if(is.null(file.name)){
        warning('File name is not set. It is set to dummy.bed')
        file.name = 'dummy.bed'
    }
    
    if(!is.null(directory)){
        file.name=paste0(directory, '/', file.name)
    } 
    
    if(any(duplicated(bed.df[['TE.ID']]))){
        warning('There are double entries in your data frame. Double entries are removed')
        bed.df <- bed.df %>% filter(!duplicated(TE.ID)) 
    }
    
    write.table(bed.df, 
                file = file.name, 
                sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
}

#======== Work with count tables ===========


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


aggregatePeaks <- function(count.matrix, annotation){
    # Multiple peaks can intersect with one gene, those are aggregated by
    # that function.
    require(tidyverse)
    
    
    count.matrix <- rownames_to_column(count.matrix, "peakName")
    
    count.matrix <- merge(count.matrix, annotation, by = 'peakName')
    
    count.matrix$peakName <- NULL
    
    count.matrix <- aggregate(.~ geneID, data = count.matrix, sum)
    
    count.matrix <- column_to_rownames(count.matrix, 'geneID')
    
    return(count.matrix)
    
}


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


#========= Preparation of correlation plots =======

getValues <- function(df){
    
    # returns the columns of interests for the correlation plots
    
    require(tidyverse)
    
    df <- df %>% 
        dplyr::select(baseMean, log2FoldChange, padj) %>% 
        rownames_to_column("geneID")
    
    return(df)
    
}

mergeDESeqRes <- function(results_Seq1, results_Seq2, seq_names){
    
    # Two DESeq results of different sequencing technologies are merged and
    # elements that contain an NA in the adjusted p-value in at least one
    # result table are removed. A sig column is added where it is stored where
    # the element is significantly differential expressed.
    
    Seq1_name <- seq_names[1]
    Seq2_name <- seq_names[2]
    
    results_Seq1 <- getValues(results_Seq1)
    results_Seq2 <- getValues(results_Seq2)
    
    names(results_Seq1) <- paste(Seq1_name, names(results_Seq1), sep = '.')
    names(results_Seq2) <- paste(Seq2_name, names(results_Seq2), sep = '.')
    
    merged.df <- merge(results_Seq1, results_Seq2, by.x = paste0(Seq1_name, '.geneID'), by.y = paste0(Seq2_name, '.geneID'))
    
    merged.df <- merged.df %>% filter(!is.na(merged.df[[paste0(Seq1_name, '.padj')]]), !is.na(merged.df[[paste0(Seq2_name, '.padj')]]))
    
    
    merged.df$sig <- ifelse(merged.df[[paste0(Seq1_name, '.padj')]] <= 0.05 & merged.df[[paste0(Seq2_name, '.padj')]] <= 0.05, paste(Seq1_name, '&', Seq2_name), 
                            ifelse(merged.df[[paste0(Seq1_name, '.padj')]] <= 0.05 & !merged.df[[paste0(Seq2_name, '.padj')]] <= 0.05, Seq1_name,
                                   ifelse(!merged.df[[paste0(Seq1_name, '.padj')]] <= 0.05 & merged.df[[paste0(Seq2_name, '.padj')]] <= 0.05, Seq2_name, 'NONE')))
    
    merged.df$sig <- factor(merged.df$sig, levels = c(paste(Seq1_name, '&', Seq2_name), Seq1_name, Seq2_name, 'NONE')) 
    
    return(merged.df)
    
}


#======== Work with deseq result tables ===========

getBioMart <- function(request = c('entrezgene_id','ensembl_gene_id', 'external_gene_name', 'entrezgene_accession')){
    
    require(biomaRt)
    
    ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 102)
    
    
    gene.annotation <- getBM(request, mart = ensembl)
    
    return(gene.annotation)
}

updateDESeqResult <- function(df, closest = NULL, featureType = 'gene'){
    
    # This function removes all features without an adjusted p-value as they are not
    # considered as expressed.
    # Additionally, the external gene names and entrez ids are added to the table. 
    # Sometimes multiple entrez id's exist per ensembl id. I removed all double 
    # ensembl ids, so that only unique ones are contained in the df. The ensembl
    # ids are back transformed to the row names.
    require(tidyverse) 
    
    if(class(df) != "data.frame"){
        df <- as.data.frame(df)
    } 
    
    df <- df[!is.na(df$padj),] 
    
    if(featureType == 'te'){
        
        df <- rownames_to_column(df, var='TE.ID')
        
        df <- merge(df, closest[c('TE.ID', 'gene.ID', 'distance')], by = 'TE.ID')
        
        df <- df %>% dplyr::rename(ensembl_gene_id = gene.ID)
        
        df <- merge(df, getBioMart(), by = 'ensembl_gene_id')
        
        df <- mutate(df, localization = case_when(distance == 0 ~ 'intragenic',
                                                TRUE ~ 'intergenic'),
                     type = case_when(padj <= 0.05 ~ 'DETE',
                                     TRUE ~ 'eTE'))
        
        #df <- column_to_rownames(df, var='TE.ID')
    }
    
    if(featureType == 'gene'){
        
        df <- rownames_to_column(df, var='ensembl_gene_id')
        
        df <- merge(df, getBioMart(), by = 'ensembl_gene_id')
        
        df <- df %>% filter(!duplicated(ensembl_gene_id))
        
        df <- column_to_rownames(df, var='ensembl_gene_id')
    }
    
    return(df)    
}


cutPvalue <- function(df, FDR_CAP = 0.00001){
# This function allows to set the p-value to a certain max value if the adjuste
# p-value is below. 
    
    df <- df %>% 
        mutate(padj = case_when(padj <= FDR_CAP ~ FDR_CAP,
                                TRUE ~ padj)) 
    
    return(df)
}


