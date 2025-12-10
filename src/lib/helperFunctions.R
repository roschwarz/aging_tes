
#----------------------------- general helper ----------------------------------

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

order_tissue <- function(df, levels = c('brain', 'skin', 'blood')){
    
    df$tissue <- factor(df$tissue, levels = levels)
    
    return(df)
}

# Calculates the z-score for a vector
# The max z-score can be capped
calc_z_score <- function(x, capped = NULL){
    
    zscore = (x - mean(x)) / sd(x)
    
    if (!is.null(capped)) {
        
        zscore <- sapply(zscore,USE.NAMES = F, function(x){
            
            if (abs(x) < capped) {
                return(x) 
            }
            
            if (x < 0 ) {
                return(-1*capped)
            } 
            return(capped)
        })
    }
    return(zscore)
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

# --------------------------- Fence plot helper --------------------------------

getFlankTSSCoordinates <- function(df, type = NULL){
    # type is up or down
    names(df) <- c('chrom', 
                   'flank.start', 
                   'flank.end', 
                   'TE.region.ID', 
                   'flank.value',
                   'flank.strand')
    
    df <- df %>% dplyr::select(TE.region.ID, flank.start, flank.end)
    
    names(df) <- c('TE.region.ID',
                   paste0(type, '.tss.start'),
                   paste0(type, '.tss.end'))
    
    return(df)
    
    
}

calcRelativePosition <- function(regions.of.interest){
    # Dose NA in flank mean that there is no intersection? 
    
    df.down <- getFlankTSSCoordinates(regions.of.interest[['down.stream']], 'down')
    df.up <- getFlankTSSCoordinates(regions.of.interest[['up.stream']], 'up')
    
    df <- unique(merge(regions.of.interest[['te.regions']], df.up, by='TE.region.ID', all.x = T))
    df <- unique(merge(df, df.down, by='TE.region.ID', all.x = T))
    
    
    df <- df %>% 
        mutate(up.tss.pos =  ifelse(TE.region.strand == '-', up.flank.end - up.tss.start, up.tss.start - up.flank.start),
               down.tss.pos =  ifelse(TE.region.strand == '-', down.flank.end - down.tss.start, down.tss.start - down.flank.start),
               up.relative.position = (up.tss.pos*100/up.window.size) - 101,
               down.relative.position =(down.tss.pos*100/down.window.size) + 101,
        )
    
    
    return(df)
    
    
}


#--------------------------------- Annotation stuff ----------------------------


create_TE_Instances <- function(teIslandRanges, teRanges, file = 'data/shared/mm10_TE_island_instances_extended.Rdata') {
    # This function intersects TE islands with TEs to get the information about the composition of TE islands. 
    # In addition, the relative position of each TE within the TE island is calculated and the ID of the TEs is
    # split to get meta information about the respective TE.
    #
    #  - add TE position within the TE Island? Could be useful for an enrichment analysis.
    
    if (file.exists(file)) {
        print("An annotation is already available, which is loaded.")
        teIslandInstances <- loadRdata(file)
    } else{
        
        teIslandInstances <-
            intersectGranger(teIslandRanges, teRanges, tab = 'all')
        
        # TEs can intersect with multiple genes, so that they occur multiple times in
        # the teRange objects. Therefore an unique is applied at the end of the pipeline
        teIslandInstances <- teIslandInstances %>%
            dplyr::select(
                "query.start" ,
                "query.end",
                "query.strand",
                "query.width",
                "query.names",
                "subject.start",
                "subject.end",
                "subject.strand",
                "subject.te.id",
                "subject.position"
            ) %>%
            filter(!duplicated(subject.te.id)) %>%
            unique()
        
        
        names(teIslandInstances) <-
            c(
                "te_island_start",
                "te_island_end",
                "te_island_strand",
                "te_island_width",
                "te_island_id",
                "te_start",
                "te_end",
                "te_strand",
                "te_id",
                "te_position"
            )
        
        # Set the start and end of TEs to the borders of TE Islands. Is this really needed for this approach?
        teIslandInstances <- teIslandInstances %>%
            mutate(
                te_start = case_when(te_start < te_island_start ~ te_island_start, .default = te_start),
                te_end = case_when(te_end > te_island_end ~ te_island_end, .default = te_end)
            )
        
        teIslandInstances <- teIslandInstances %>%
            group_by(te_island_id) %>%
            mutate(
                te_rel_start = case_when(
                    te_island_strand == "+" ~ ((te_start - te_island_start) * 100 / te_island_width),
                    .default = (te_island_end - te_end) *
                        100 / te_island_width
                ),
                te_rel_end = case_when(
                    te_island_strand == "+" ~ ((te_end - te_island_start) * 100 / te_island_width),
                    .default = (te_island_end - te_start) *
                        100 / te_island_width
                )
            ) %>%
            ungroup() %>%
            blackRcloud::splitTEID('te_id')
        
        save(teIslandInstances, file = file)
        
    }
    
    return(teIslandInstances)
    
}

#---------------------------- prepare df for plots -----------------------------


prepare_Composition <- function(df){
    # - Collect all TEs that are not part of TEs of interest in Other.
    
    order_count <- df %>%
        filter(order %in% orders.of.interest) %>% 
        dplyr::count(order)
    
    order_count$per <- round(order_count$n / sum(order_count$n), 2)
    
    order_count$labels <-
        paste(round(order_count$n / sum(order_count$n) * 100, 2), "%")
    
    return(order_count)
}

prepare_super_composition <- function(df) {
    df %>%
        dplyr::select(order, super_family) %>%
        dplyr::count(super_family, order) %>%
        mutate(per = n / sum(n) * 100) %>%
        group_by(order) %>%
        mutate(super_family = reorder(super_family, per))
    
}

prepare_position_composition <- function(df) {
    # Is this function still needed?
    # I think a more generalized version could make sense, because the te island definition is changed at the
    # end of the document. There the CAGE-Peak and the last Quant-peak define the ends. Therefore, the overlap
    # of Tes can be different.
    
    blackRcloud::fullIntersectionResult(te_island_gRange,
                                        teRanges) %>% 
        #filter(subject.position == 'intergenic') %>% 
        splitTEID('subject.te.id') %>% 
        
        filter(order %in% order_of_interest)
    
    density_df <-
        TE_Island_Overlap %>% dplyr::select(
            "query.start" ,
            "query.end",
            "query.strand",
            "query.width",
            "query.teisland_id",
            "subject.start",
            "subject.end",
            "subject.strand",
            "subject.te.id"
        )
    
    names(density_df) <-
        c(
            "tei_start",
            "tei_end",
            "tei_strand",
            "tei_width",
            "tei_id",
            "te_start",
            "te_end",
            "te_strand",
            "te_id"
        )
    
    
    density_df <- density_df %>%
        mutate(
            te_start = case_when(te_start < tei_start ~ tei_start, .default = te_start),
            te_end = case_when(te_end > tei_end ~ tei_end, .default = te_end)
        )
    
    density_df <- density_df %>%
        group_by(tei_id) %>%
        mutate(
            te_rel_start = case_when(
                tei_strand == "+" ~ ((te_start - tei_start) * 100 / tei_width),
                .default = (tei_end - te_end) *
                    100 / tei_width
            ),
            te_rel_end = case_when(
                tei_strand == "+" ~ ((te_end - tei_start) * 100 / tei_width),
                .default = (tei_end - te_start) *
                    100 / tei_width
            )
        ) %>%
        ungroup() %>%
        blackRcloud::splitTEID('te_id')
    
    
    
}


#--------------------------------- Stuff to Sort -------------------------------
# getCondition <- function(header){
#     
#     # Returns a df with sample names (finally the row names of the df) and the 
#     # condition. The header needs to be separated by an underscore and the 
#     # condition needs to be at the second position. Usable for DESeq2.
#     
#     df <- data.frame(SampleID = header)
#     df$condition <- sapply(df$SampleID, function(x) strsplit(x, "_")[[1]][2])
#     
#     rownames(df) <- df$SampleID
#     df$SampleID <- NULL
#     
#     df$condition <- as.factor(df$condition)
#     return(df)
# }


# takes the count matrixs and filters for all elements of the list where the 
# TE.ID starts with a chr, which usually identifies the TE id. Gene ids start
# with ENS
filterForTEs <- function(count.matrix){
    
    require(tidyverse)
    
    count.matrix %>% 
        rownames_to_column(var = 'TE.ID') %>% 
        filter(grepl('^chr', TE.ID)) %>% 
        column_to_rownames(var = 'TE.ID')
}

# takes the count matrixs and filters for all elements of the list where the 
# TE.ID starts with a chr, which usually identifies the TE id. Gene ids start
# with ENS
filterSalmonTEcounts <- function(count.matrix, type = 'te'){
    
    require(tidyverse)
    
    if(type == 'te'){
        count.matrix <- count.matrix %>% 
            rownames_to_column(var = 'TE.ID') %>% 
            filter(grepl('^chr', TE.ID)) %>% 
            column_to_rownames(var = 'TE.ID')
        
    }else if(type == 'gene'){
        count.matrix <- count.matrix %>% 
            rownames_to_column(var = 'TE.ID') %>% 
            filter(grepl('^ENS', TE.ID)) %>% 
            column_to_rownames(var = 'TE.ID')
        
    }else{
        stop("selected type not known")
    }
}
# doDEseq <- function(count.matrix, condition, paral=TRUE, reference = NULL, TEspecific = FALSE){
#     
#     # Condition data frame is created out of the header name. This is highly
#     # project specific. Peaks with less than 11 reads in sum across all 
#     # samples are removed. Possibly, you can set that threshold a bit higher.
#     # Check your notes for that, I recently read a paper where a threshold
#     # was written.
#     
#     require(DESeq2)
#     
#     # condition <- getCondition(names(count.matrix))
#     
#     if(TEspecific){
#         count.matrix <- filterForTEs(count.matrix)
#     }
#     
#     dds <- DESeqDataSetFromMatrix(countData = count.matrix,
#                                   colData = condition,
#                                   design = ~condition)
#     
#     
#     dds <- dds[rowSums(counts(dds)) >= 10, ]
#     
#     if(!is.null(reference)){
#         
#         dds$condition <- relevel(dds$condition, ref = reference)
#     }
#     
#     dds <- DESeq(dds, parallel = paral)
#     
#     return(dds)
# }

# Extract the results from a dds object and returns the result table filtered
# for instances with an adjusted p-value. An lfc shrinkage is done by default
# but can turned of with lfcShrink = FALSE.
# getDEseqResults <- function(dds, lfcShrink = TRUE, coefficient = NULL){
#     
#     # Implement a trap if something is missing
#     # if(lfcShrink & is.null(coefficient)){
#     #     stop()
#     # }
#     
#     require(DESeq2) 
#     require(tidyverse)
#     
#     if(!(lfcShrink)){
#         
#         deseq.res <- as.data.frame(results(dds)) %>% 
#             filter(!is.na(padj)) 
#         
#         return(deseq.res)
#         
#     } 
#     
#     deseq.res <- DESeq2::lfcShrink(dds, 
#                                    coef = coefficient,
#                                    parallel = F,
#                                    res = DESeq2::results(dds),
#                                    type = 'apeglm')
#     
#     deseq.res <- as.data.frame(deseq.res) %>% filter(!is.na(padj)) 
#     
#     return(deseq.res)
# }


# deseqTemplate <- function(list.of.files, preFilter = FALSE){
#     # This function takes a named list of count matrices and do DESeq2
#     # It returns a named list with the count.matrices (cm), dds objects (dds), 
#     # and deseq results (dres) per entry in the input list.
#     # To-Do:
#     #   - add arguments to turn on and off the parallel stuff
#     
#     
#     count.matrices <- sapply(names(list.of.files), function(x) getCountTable(list.of.files[[x]], filter = F)[['counts']])
#     
#     # In my RNA-Seq data are more columns than necessary, so that they will
#     # be removed here.
#     count.matrices <- sapply(names(list.of.files), function(x) {
#         
#         counts <- count.matrices[[x]]
#         counts <- counts[names(counts[,grepl("*unique.sorted.bam", names(counts))])]
#         counts
#     })
#     
#     count.matrices <- sapply(names(count.matrices), function(x) updateHeader(count.matrices[[x]]))
#     
#     if(preFilter){
#         # TE data can contain a bunch of lines without any read which can kill 
#         # the memory, therefore a filter step could be useful
#         count.matrices <- sapply(names(count.matrices), 
#                                  function(x) {(count.matrices[[x]][rowSums(count.matrices[[x]]) >= 10,])}
#                                  )
#     }
#     
#     dds.objects <- sapply(names(count.matrices), 
#                           function(x) doDEseq(count.matrices[[x]], 
#                                               paral = F))
#     
#     #deseq.results <- sapply(names(dds.objects), function(x) DESeq2::results(dds.objects[[x]]))
#     deseq.results <- sapply(names(dds.objects), simplify = F, function(x){
#         
#         df <- DESeq2::lfcShrink(dds.objects[[x]],
#                           coef = "condition_old_vs_young",
#                           parallel = F,
#                           res = DESeq2::results(dds.objects[[x]]),
#                           type = 'apeglm')
#         df <- as.data.frame(df) %>% 
#             filter(!is.na(padj))
#         
#         return(df)
#         }
#     )
#     
#     return(
#         list(
#             cm = count.matrices,
#             dds = dds.objects,
#             dres = deseq.results
#         )
#     )
# }
# 
# deseqAgg <- function(deseq.temp.object, peak.annotation, feature.type = 'gene'){
#     
#     # Multiple-Peaks can intersect one gene/te. Here the peaks are aggregated per
#     # gene/te and a gene/te-wise DESeq analysis is done. The results are added
#     # to the list of the incoming object. All elements that have no adjusted p-value
#     # are removed from the result table.
#     
#     conditions <- names(deseq.temp.object[['cm']])
#     
#     agg.cm <- sapply(conditions, 
#                       function(x) aggregatePeaks(deseq.temp.object[['cm']][[x]], 
#                                                  peak.annotation[[x]]))
#     
#     agg.dds <- sapply(conditions, function(x) doDEseq(agg.cm[[x]]))
#     
#     agg.dres <- sapply(conditions, function(x) DESeq2::results(agg.dds[[x]]))
#     
#     
#     
#    deseq.temp.object[[paste0(feature.type, '.cm')]] <- agg.cm
#    deseq.temp.object[[paste0(feature.type, '.dds')]] <- agg.dds
#    deseq.temp.object[[paste0(feature.type, '.dres')]] <- agg.dres
#    
#     
#     return(deseq.temp.object) 
# }


