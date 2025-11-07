# The script calculates z-scores of vst of features of a certain region within the protocadherin beta cluster.
# Additionally, the log2 fold changes and FDR values of the respective features in the regions
# are extracted from the DESeq2-result tables. The data is stored in two csv-file.
#
# INPUT:
# - dds object of te island and genes (te_island_rna_seq_dds.Rdata)
# - deseq result tables of te island and genes (02_deseq_results_genes.csv, 02_deseq_results_te_island_extended.csv)
# - annotations (teIslandRanges and geneRanges)
#
# OUTPUT:
#   ./tables/pcdhb15_locus_heatmap_zscore_vst.csv
#   ./tabels/pcdhb15_locus_heatmap_log2fc_fdr.csv

# Author: Robert Schwarz
# Created: 2025-11-07

# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()
aging_tes::load_rna_seq_env()
aging_tes::load_annotations()

load_gene_ranges()
load_te_island_annotation()

# ========================= Get variance stabilized counts for TE islands and genes ==========================

# -------------------------------- TE island -----------------------------------------------------------------

counts_TE_island <- loadRdata(paste0(rna_seq_deseq_dir,
                                     "te_island_rna_seq_dds.Rdata"))

vst_counts_te_island_brain <- data.frame(getVarianceStabilizedData(counts_TE_island$brain))

vst_counts_te_island_brain <- rownames_to_column(vst_counts_te_island_brain, var = "te_island_id")

# -------------------------------- Genes ---------------------------------------

counts_gene <- loadRdata(paste0(rna_seq_deseq_dir, deseq_dds_gene))

vst_counts_gene_brain <- data.frame(getVarianceStabilizedData(counts_gene$brain))

vst_counts_gene_brain <- rownames_to_column(vst_counts_gene_brain,
                                             var = "ensembl_gene_id")

# ====================== Get TE regions and genes of interest ==================

# Handpicked in the genome browser for figure 5 in the publication
te_islands_of_interest <- c("TE_Cluster_728951",
                         "TE_Cluster_728962",
                         "TE_Cluster_728983",
                         "TE_Cluster_728986",
                         "TE_Cluster_728991")

Pcdhb1_coordinates <- transGrange(geneRanges) %>%
    filter(external_gene_name == "Pcdhb1")

Pcdhb22_coordinates <- transGrange(geneRanges) %>%
    filter(external_gene_name == "Pcdhb22")

genes <- transGrange(geneRanges) %>%
    filter(seqnames == "chr18",
           start >= Pcdhb1_coordinates$start,
           end <= Pcdhb22_coordinates$end,
           strand == "+")

# ==================== Filter count tables and calc z scores ===================

# --------------------------- TE Islands ----------------------------------------

pro_locus_vst_counts_te_island_brain <-  vst_counts_te_island_brain %>%
    filter(te_island_id %in% te_islands_of_interest) %>%
    gather(key = "sample",
           value = "vst",
           names(vst_counts_te_island_brain)[2:length(vst_counts_te_island_brain)])


pro_locus_vst_counts_te_island_brain <- pro_locus_vst_counts_te_island_brain %>%
    separate("sample",
             c("age.name", "age", "sample.id", "tissue"),
             sep = "[_]",
             remove = FALSE) %>%
    dplyr::select(-age.name)


pro_locus_vst_counts_te_island_brain <- pro_locus_vst_counts_te_island_brain %>%
    dplyr::select(te_island_id, sample, vst) %>%
    spread(sample, vst) %>%
    column_to_rownames(var = "te_island_id")

# --------------------------- Genes ----------------------------------------

pro_locus_vst_counts_gene_brain <- vst_counts_gene_brain %>%
    filter(ensembl_gene_id %in% genes$ensembl_gene_id)  %>%
    gather(key = "sample",
           value = "vst",
           names(vst_counts_gene_brain)[2:length(vst_counts_gene_brain)])

pro_locus_vst_counts_gene_brain <- pro_locus_vst_counts_gene_brain %>%
    separate("sample",
             c("age.name", "age", "sample.id", "tissue"),
             sep = "[_]",
             remove = FALSE) %>%
    dplyr::select(-age.name)

pro_locus_vst_counts_gene_brain <- pro_locus_vst_counts_gene_brain %>%
    dplyr::select(ensembl_gene_id, sample, vst) %>%
    spread(sample, vst)

pro_locus_vst_counts_gene_brain <- merge(pro_locus_vst_counts_gene_brain,
                          transGrange(geneRanges) %>%
                              dplyr::select(ensembl_gene_id,
                                            external_gene_name),
                          by = "ensembl_gene_id") %>%
    dplyr::select(-ensembl_gene_id) %>%
    column_to_rownames("external_gene_name")


# --------------------------- Merge and sort tables ----------------------------

vst_counts_pro_locus <- rbind(pro_locus_vst_counts_te_island_brain, pro_locus_vst_counts_gene_brain)

order_of_elements <- rbind(
    
    genes %>%
        dplyr::select(-ensembl_gene_id, -gene_biotype) %>%
        dplyr::rename(names = external_gene_name),
    
    transGrange(te_island_5primeRanges) %>%
        filter(names %in% te_islands_of_interest)
    
)

vst_counts_pro_locus <- vst_counts_pro_locus[order_of_elements[order(order_of_elements$start),
                                             "names"], ]


vst_counts_pro_locus <- t(apply(as.matrix(vst_counts_pro_locus),
                       1,
                       function(x) calc_z_score(x, capped = 2)))

write.table(vst_counts_pro_locus,
            file = paste0(table_dir, "pcdhb15_locus_heatmap_zscore_vst.csv"),
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = ","
)

# ======================== Get quantification numbers ========================================================

gene_deseq <- read.csv(paste0(table_dir, "02_deseq_results_genes.csv")) %>%
    filter(external_gene_name %in% genes$external_gene_name,
           tissue == "brain") %>%
    dplyr::select(external_gene_name, log2FoldChange, padj) %>% 
    dplyr::rename(id = external_gene_name)

te_island_deseq <- read.csv(paste0(table_dir, "02_deseq_results_te_island_extended.csv")) %>%
    filter(te_island_id %in% te_islands_of_interest, tissue == "brain") %>%
    dplyr::select(te_island_id, log2FoldChange, padj) %>%
    dplyr::rename(id = te_island_id)


df <- rbind(gene_deseq, te_island_deseq)

df <- df[match(rownames(vst_counts_pro_locus), df$id),]


write.table(df,
            file = paste0(table_dir, "pcdhb15_locus_heatmap_log2fc_fdr.csv"),
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = ","
)




