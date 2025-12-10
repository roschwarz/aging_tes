if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

source("src/lib/helperFunctions.R")

aging_tes::load_te_island_env()
aging_tes::load_annotations()
aging_tes::load_rna_seq_env()
aging_tes::load_cage_seq_env()
aging_tes::load_quant_seq_env()

load_te_island_annotation()
load_te_ranges()
load_cage_peak_annotation()
load_quant_peak_annotation()

rna_deseq_results <- blackRcloud::loadRdata(paste0(rna_seq_deseq_dir, deseq_results_te))

te_island_instances <- create_TE_Instances(te_islandRanges, teRanges)

expr_te_islands <- sapply(names(rna_deseq_results), simplify = F, function(x){
    
    expr_te_rna <- rownames(rna_deseq_results[[x]] %>% filter(!is.na(padj)))
    
    expr_island <- te_island_instances %>%
        filter(te_id %in% expr_te_rna) %>%
        distinct(te_island_id)
    
    return(expr_island$te_island_id)
})

tss_carrying_TEIsland <- sapply(names(cageRanges), simplify = F, USE.NAMES = TRUE, function(x){
    
    intersectGranger(cageRanges[[x]], te_island_5primeRanges, 'all')
    
})



# Identify TE islands that carry at least one express TE and overlap with CTSS peaks
tss_carrying_expr_TEIslands <- sapply(names(cageRanges), simplify = F, USE.NAMES = TRUE, function(x){
    
   intersectGranger(cageRanges[[x]], te_island_5primeRanges, 'all') %>% 
        filter(subject.names %in% expr_te_islands[[x]]) %>%
        filter(!duplicated(subject.names))
    
})

tts_carrying_TEIsland <- sapply(names(quantRanges), simplify = F, USE.NAMES = TRUE, function(x){
    
    intersectGranger(quantRanges[[x]], te_island_3primeRanges, 'all') 
    
})

# Identify TE islands that carry at least one express TE and overlap with QTTS peaks
tts_carrying_expr_TEIslands <- sapply(names(quantRanges), simplify = F, USE.NAMES = T, function(x){
    
    df <- intersectGranger(quantRanges[[x]], te_island_3primeRanges, 'all') %>% 
        filter(subject.names %in% expr_te_islands[[x]]) %>%
        filter(!duplicated(subject.names))
    
    bed_file <- df %>% dplyr::select(subject.seqnames,
                                     subject.start,
                                     subject.end,
                                     subject.names,
                                     subject.width,
                                     subject.strand)
    
    file_name <- paste0('data/shared/', x, '_tss_carrying_expressed_TE_Island.bed')
    
    write.table(bed_file, file = file_name, sep = '\t', col.names = F, row.names = F, quote = F)
    return(df)
    
    
})

# The definition of TE island transcripts is the intersection of a TE island with a CAGE-, RNA-, and Quant-
# signal. Therefore, I just need to filter the TSS carrying TE islands for TE island ids that are contained
# in the set of TE islands that contain a Quant-signal and are expressed (sto_expr).
teIslands <- sapply(tissues, 
                    simplify = F, 
                    USE.NAMES = T, 
                    function(tissue){
                        
                        tss_carrying_expr_TEIslands[[tissue]] %>% 
                            filter(subject.names %in% tts_carrying_expr_TEIslands[[tissue]]$subject.names) %>% 
                            filter(!duplicated(subject.names))
                        
                    })


# Next, I want to annotate canonical TE Island transcripts, which means I want to set the start coordinate to 
# the start coordinate of the left most CAGE peak and end the stop coordinate to the end coordinate of the 
# right most Quant peak.
# 
# Filter for TE Island transcripts with a width > 0, as when the Quant-peak is located in front of the CAGE-
# peak you will get a negative length and this is not a real TE island transcript.
cannonical_teIslands <- do.call('rbind', sapply(tissues, 
                                                simplify = F, 
                                                USE.NAMES = T, 
                                                function(tissue){
                                                    
                                                    # collect the first cage peak in a TE island
                                                    cage_tmp <- tss_carrying_TEIsland[[tissue]] %>% 
                                                        group_by(subject.names) %>% 
                                                        filter(case_when(subject.strand == "-" ~ query.end == (max(query.end)),
                                                                         subject.strand == "+" ~ query.start == (min(query.start))
                                                        )
                                                        )
                                                    
                                                    names(cage_tmp) <- str_replace_all(names(cage_tmp), 
                                                                                       c('query' = 'cage', 'subject' = 'te_region'))
                                                    
                                                    # collect the last Quant peak in a TE island
                                                    quant_tmp <- tts_carrying_TEIsland[[tissue]] %>% 
                                                        group_by(subject.names) %>% 
                                                        filter(case_when(subject.strand == "-" ~ query.end == (min(query.end)),
                                                                         subject.strand == "+" ~ query.start == (max(query.start))
                                                        )
                                                        ) 
                                                    
                                                    
                                                    names(quant_tmp) <- str_replace_all(names(quant_tmp), 
                                                                                        c('query' = 'quant', 'subject' = 'te_region'))
                                                    
                                                    quant_tmp <- quant_tmp %>% dplyr::select(-c("te_region.seqnames", 
                                                                                                "te_region.start",
                                                                                                "te_region.end", 
                                                                                                "te_region.width",
                                                                                                "te_region.strand"))
                                                    
                                                    df <- merge(cage_tmp, quant_tmp, by = 'te_region.names')
                                                    
                                                    df_x <- df %>% 
                                                        filter(te_region.names %in% expr_te_islands[[tissue]]) %>% 
                                                        mutate(coord = case_when(te_region.strand == "+" ~ paste0(te_region.seqnames,
                                                                                                                  ":",
                                                                                                                  cage.start,
                                                                                                                  "-",
                                                                                                                  quant.end),
                                                                                 te_region.strand == "-" ~ paste0(te_region.seqnames,
                                                                                                                  ":",
                                                                                                                  quant.start,
                                                                                                                  "-",
                                                                                                                  cage.end)),
                                                               te_island.start = case_when(te_region.strand == "+" ~ cage.start,
                                                                                           te_region.strand == "-" ~ quant.start),
                                                               te_island.end = case_when(te_region.strand == "+" ~ quant.end,
                                                                                         te_region.strand == "-" ~ cage.end),
                                                               te_island.width = te_island.end - te_island.start,
                                                               tissue = tissue
                                                        ) %>% 
                                                        filter(te_island.width > 0)
                                                    
                                                }))


write.csv(cannonical_teIslands,
          file = 'results/TEItx/canonical_te_island_transcripts.csv',
          row.names = F,
          col.names = T,
          quote = F)
