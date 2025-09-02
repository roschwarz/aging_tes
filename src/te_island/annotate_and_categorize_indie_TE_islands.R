# ==============================================================================
# Annotation and characterization of individually expressed TE islands
# ==============================================================================
#
# This script takes the by distance annotated TE islands (5'prime extended; 
# see) and overlaps it with CAGE-peaks and RNA-Seq signals to identify individually
# expressed TE islands. Individually expressed TE islands are defined as an TE island
# that overlaps with at least one CAGE-peak (position doesn't matter) and at least
# one expressed TE (adjusted p-value is not NA).
#
# In addition, the TE island is characterized by the number of TEs overlapping the
# TE island. TE islands are subdivided into three categories:
#   - single - only one TE
#   - double - two TEs
#   - multi - more the two TEs
#
# ==============================================================================

# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_rna_seq_env()
aging_tes::load_cage_peak_annotation()
aging_tes::load_annotations()
aging_tes::load_te_island_env()

load_te_ranges()
load_te_island_annotation()
load_gene_ranges()

# Collect TE instances of TE-regions to get information about the composition of the TE
# regions later in the analysis.
logmsg("Overlap 5'-extended TE islands with te instances")
te_island_instances <- blackRcloud::intersectGranger(te_island_5primeRanges, teRanges, tab = 'all')

# TEs can intersect with multiple genes, so that they occur multiple times in
# the teRange object. Therefore a unique is applied at the end of the pipeline
te_island_instances <- te_island_instances %>% 
    dplyr::select(query.names, subject.te.id, subject.position) %>% 
    filter(!duplicated(subject.te.id)) %>% 
    unique()

names(te_island_instances) <- c("te_island_id", "te_id", "te_position")


#####################################
# individually expressed TE islands #
#####################################

# Usage of CAGE- and RNA-seq to get individually expressed TE islands. 
#  1. Intersect the extended TE islands (te_island_5primeRanges) with CAGE-peaks to get potential 
#     individually expressed TE islands.
#  2. Get individually expressed TE islands that also contain at least one expressed TE instance that 
#     is detected via RNA-Seq. An TE instance is considered as expressed when an adjusted p-value after
#     the DESeq analysis is available.

# Cage intersection
logmsg("Determine individually TE islands.")
indie_te_islands <- sapply(names(cageRanges), 
                          simplify = F, 
                          USE.NAMES = T, 
                          function(x){
                              
                              blackRcloud::intersectGranger(cageRanges[[x]], te_island_5primeRanges, tab = 'all')
                              
                          })

# RNA intersection
load_deseq_results() # loads rna_deseq_res_te
logmsg("Determine individually expressed TE islands.")
expressed_te_islands <- sapply(names(rna_deseq_res_te), simplify = F, function(x){
    
    expressed_TEs_rna <- rownames(rna_deseq_res_te[[x]]  %>% filter(!is.na(padj)))
    
    exp_islands <- te_island_instances %>% 
        filter(te_id %in% expressed_TEs_rna) %>% 
        dplyr::pull('te_island_id')
    
    # regions Ids occur multiple times because of the composition of the islands
    # However, only one id is needed to get the expressed TE island 
    return(unique(exp_islands))
})

# filter the cage intersecting TE islands for those that intersect with 
# expressed TEs detected via RNA-Seq
indie_te_islands <- sapply(names(indie_te_islands), simplify = F, function(x){
    
    te_island <- indie_te_islands[[x]]
    exp_te <- expressed_te_islands[[x]]
    
    te_island %>% 
        filter(subject.names %in% exp_te) %>% 
        filter(!duplicated(subject.names))
    
})


## Write bed-files of individually expressed TE islands
logmsg('Write .bed-files of individually expressed TE islands.')
sapply(names(indie_te_islands), function(x){
    
    df <- indie_te_islands[[x]]
    
    bed_file <- GRanges(seqnames = df$subject.seqnames,
                        ranges = IRanges(start = df$subject.start,
                                         end = df$subject.end,
                                         names = df$subject.names),
                        strand = df$subject.strand,
                        te_island_id = df$subject.names,
                        value = 1)
    
    bed_file <- transGrange(bed_file)
    
    bed_file <- bed_file[c('seqnames', 'start', 'end', 'te_island_id',  'value', 'strand')]
    
    filename = indie_te_island_bed[[x]]
    logmsg(filename)
    write.table(bed_file, file = filename, sep = '\t', col.names = F, row.names = F, quote = F)
    
})

# ==============================================================================
# Characterization of TE islands
# ==============================================================================
#
# Composition of TE islands #
# Take the TE region instance table and filter for all TE regions that are 
# individually expressed. Determine the position of the instance within the TE
# island and categorize TE islands into single (consists of one instance),
# double (consists of two instances) and multiple (consists of more than 2 instances).
#
# ==============================================================================
logmsg('Determine the composition of the individually expressed TE islands')
indie_te_islands_composition <- sapply(names(indie_te_islands), simplify = F, function(x){
    
    indie_te_islands_instances <- te_island_instances %>% 
        filter(te_island_id %in% indie_te_islands[[x]][['subject.names']])
    
    
    names(indie_te_islands_instances) <- c('te_island_id', 'te_id', 'te_genomic_pos')
    
    indie_te_islands_instances <- splitTEID(indie_te_islands_instances, 'te_id')
    
    indie_te_islands_instances <- indie_te_islands_instances %>% 
        group_by(te_island_id) %>% 
        mutate(member = n(),
               position = ifelse(strand == '-', rev(1:n()), 1:n()),
               island_type = ifelse(member == 1, 'single', 
                                    ifelse(member == 2, 'double', 'multiple')))
    
    
    return(indie_te_islands_instances)
})


indie_te_island_super_family_composition <- sapply(names(indie_te_islands_composition), simplify = F, function(x){
    
    df <- indie_te_islands_composition[[x]]
    
    df <- df %>% 
        filter(member > 2) %>% 
        mutate(super_family = case_when(super_family == "Alu" ~ "B1", .default = super_family)) %>% 
        group_by(te_island_id) %>% 
        mutate(te_position = ifelse(position == 1, 'first', ifelse(position == max(position), 'last', 'body'))) %>% 
        ungroup() %>% 
        dplyr::count(super_family, te_position) %>% 
        group_by(te_position) %>% 
        mutate(percent = n/sum(n))
    
    df$te_position <- factor(df$te_position, levels = c('first', 'body', 'last'))
    
    
    df <- df %>% 
        dplyr::select(-n) %>% 
        filter(super_family != "hAT-Charlie") %>% 
        spread(key = te_position, value = percent, fill = 0) %>% 
        column_to_rownames(var = 'super_family') %>% 
        filter(rowSums(.) > 0.05)
    
    return(as.matrix(df))
    
})


indie_te_island_categorized <- list("instance" = indie_te_islands_composition,
                                 "super_fam" = indie_te_island_super_family_composition)

logmsg(paste0("Store the categorized individually expressed TE islands (", tables_and_co,"indie_te_island_categorized.Rdata)"))
save(indie_te_island_categorized,
     file = paste0(tables_and_co, "indie_te_island_categorized.Rdata"))


# ======================================================================
# TE island gene association 
# ======================================================================

logmsg('Overlap gene annotation with 5 prime extended te islands.')

te_islands_genes_df <- intersectGranger(geneRanges, te_island_5primeRanges, 'all') %>% 
    dplyr::select(query.ensembl_gene_id, query.external_gene_name, subject.names) %>% 
    dplyr::rename(ensembl_gene_id = query.ensembl_gene_id, 
                  external_gene_name = query.external_gene_name, 
                  te_island_id = subject.names) 

logmsg(paste0("load DESeq table for te island:", table_dir, "02_deseq_results_te_island.csv"))
deseq_results <- data.table::fread(paste0(table_dir, "02_deseq_results_te_island.csv"))

logmsg("Categorization of TE islands with respect to gene localization.")
deseq_results_extended <- merge(deseq_results, te_islands_genes_df, by = 'te_island_id', all.x = TRUE) %>% 
    mutate(gene_association = case_when(is.na(ensembl_gene_id) ~ 'intergenic', .default = 'intragenic'))

logmsg("Categorization of TE islands with to individual expression.")
# Combine all indie_te tables in one DataFrame, with tissue-Info
df_all_indie <- imap_dfr(indie_te_islands, ~ mutate(.x, tissue = .y)) %>% 
    dplyr::select(subject.names, tissue) %>% 
    dplyr::rename(te_island_id = subject.names) %>% 
    mutate(individually = TRUE)

deseq_results_extended <- deseq_results_extended %>% 
    left_join(df_all_indie, by = c("te_island_id", 'tissue')) %>% 
    mutate(individually = replace_na(individually, FALSE))

logmsg(paste0("load DESeq table for genes:", table_dir, "02_deseq_results_genes.csv"))
deseq_results_genes <- data.table::fread(paste0(table_dir, "02_deseq_results_genes.csv")) %>% 
    #distinct(tissue, ensembl_gene_id, baseMean, log2FoldChange, padj, tissue, gene_biotype ) %>% 
    dplyr::select(tissue, ensembl_gene_id, baseMean, log2FoldChange, padj, tissue, gene_biotype ) %>% 
    dplyr::rename(gene_baseMean = baseMean,
                  gene_log2FoldChange = log2FoldChange,
                  gene_padj = padj)

deseq_results_extended <- deseq_results_extended %>% 
    left_join(deseq_results_genes, by = c("ensembl_gene_id", 'tissue'))

write.csv(deseq_results_extended,
          file = paste0(table_dir, "02_deseq_results_te_island_extended.csv"))
