# Description: Find overlaps between TE islands and biotypes of genes
# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_annotations()
aging_tes::load_cage_seq_env()
load_cage_peak_annotation()
load_te_ranges()

# === Transcription start sites of genes using biomaRt ===

# 1) Mouse-TSS from ensembl 
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "mmusculus_gene_ensembl",
                      version = 102)  # Ensembl Mouse v102

tss_df <- getBM(
    mart = ensembl,
    attributes = c("ensembl_gene_id", "mgi_symbol", "ensembl_transcript_id",
                   "chromosome_name", "strand", "transcription_start_site", "external_gene_name")
)


# only keep standard chromosomes
# mitochondrial chromsome ignored
std_chr <- c(as.character(1:19), "X", "Y")
tss_df <- tss_df %>% filter(chromosome_name %in% std_chr)

# UCSC-Style Seqnames ('chr' Prefix; MT -> chrM)
tss_df <- tss_df %>%
    mutate(seqname = ifelse(chromosome_name == "MT", "chrM", paste0("chr", chromosome_name)))

# all TSS, multiple TSS per gene are possible
tss_tx_gr <- GRanges(
    seqnames = tss_df$seqname,
    ranges   = IRanges(start = tss_df$transcription_start_site,
                       end   = tss_df$transcription_start_site),
    strand   = ifelse(tss_df$strand == 1, "+", "-"),
    ensembl_gene_id       = tss_df$ensembl_gene_id,
    mgi_symbol            = tss_df$mgi_symbol,
    ensembl_transcript_id = tss_df$ensembl_transcript_id
)

# === Find overlaps between TE and TSS of genes ===
# ~4 Million TEs, 5707 do overlap with at least one gene Tss
# 6395 gene TSS overlap with at least one TE
overlap_results <- process_overlapping(query_ranges = teRanges, 
                                       subject_ranges = tss_tx_gr,
                                       col_query_to_subj = 'te.id') 

df_te <- data.frame(te_id = overlap_results$query_with_hit$te.id,
                    ensembl_gene_id = mcols(overlap_results$query_with_hit)$ensembl_gene_id) %>% 
    blackRcloud::splitTEID("te_id") %>%
    filter(order %in% orders.of.interest)

te_stat_results <- analyze_te_statistics(df_te, colors = order.color)

te_stat_results$plots$super_family_barplot
te_stat_results$plots$kimura_distribution_histogram


# === Look for cage detected TSSs within the TE-derived gene TSSs ===

# Brain

brain_cage <- cageRanges$brain

# 6395 gene TSS overlap with at least one TE
# 444 of these TE derived gene TSS overlap with at least one CAGE peak in brain
overlap_results_brain <- process_overlapping(overlap_results$subject_hits, brain_cage)

df_te_brain <- data.frame(te_id = overlap_results_brain$query_with_hit$te.id,
                    ensembl_gene_id = mcols(overlap_results_brain$query_with_hit)$ensembl_gene_id) %>% 
    blackRcloud::splitTEID("te_id") %>%
    filter(order %in% orders.of.interest)

te_stat_brain <- analyze_te_statistics(df_te_brain, colors = order.color)
te_stat_brain$plots$kimura_distribution_histogram

# Skin

skin_cage <- cageRanges$skin

overlap_results_skin <- process_overlapping(overlap_results$subject_hits, skin_cage)

df_te_skin <- data.frame(te_id = overlap_results_skin$query_with_hit$te.id,
                    ensembl_gene_id = mcols(overlap_results_skin$query_with_hit)$ensembl_gene_id) %>% 
    blackRcloud::splitTEID("te_id") %>%
    filter(order %in% orders.of.interest)

te_stat_skin <- analyze_te_statistics(df_te_skin, colors = order.color)

# Blood

blood_cage <- cageRanges$blood

overlap_results_blood <- process_overlapping(overlap_results$subject_hits, blood_cage)

df_te_blood <- data.frame(te_id = overlap_results_blood$query_with_hit$te.id,
                    ensembl_gene_id = mcols(overlap_results_blood$query_with_hit)$ensembl_gene_id) %>% 
    blackRcloud::splitTEID("te_id") %>%
    filter(order %in% orders.of.interest)

te_stat_blood <- analyze_te_statistics(df_te_blood, colors = order.color)


#### Combine figures

super_fam_derived_tss <- ggarrange(te_stat_results$plot$super_family_barplot +
              ggtitle("TE-derived gene TSS (background)") +
              theme(plot.title = element_text(hjust = 0.5)),
          te_stat_brain$plot$super_family_barplot +
              ggtitle("TE-derived gene TSS (brain CAGE peaks)") +
              theme(plot.title = element_text(hjust = 0.5)),
          te_stat_skin$plot$super_family_barplot +
              ggtitle("TE-derived gene TSS (skin CAGE peaks)") +
              theme(plot.title = element_text(hjust = 0.5)),
          te_stat_blood$plot$super_family_barplot +
              ggtitle("TE-derived gene TSS (blood CAGE peaks)") +
              theme(plot.title = element_text(hjust = 0.5)),
          ncol = 2, nrow = 2,
          labels = c("A", "B", "C", "D"))


meta <- list(name = 'super_fam_derived_tss_for_genes',
             description = 'Barplots showing the super family composition of TEs that overlap with gene TSS (A) and those that overlap with CAGE peaks in brain (B), skin (C) and blood (D). Only TE orders of interest are shown.',
             tags = c('CAGE-Seq', 'TEs', 'TSS'),
             parameters = list(tissues = c('brain', 'skin', 'blood')),
             script = 'TE_derived_gene_TSSs.R'
)

fig_index(plot = super_fam_derived_tss,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 25,
          height = 17,
          dpi = 300,
          format = 'pdf')



# === Find overlaps between TE and gene TSS (most extrem TSS per gene) ===
# Rule: per gene the most extrem TSS with respect to their strand (Plus: min; Minus: max)
tx_tbl <- as.data.frame(tss_tx_gr) %>%
    group_by(ensembl_gene_id, mgi_symbol, seqnames, strand) %>%
    summarise(gene_tss = ifelse(first(strand) == "+",
                                min(start),
                                max(start)),
              .groups = "drop")


gene_tss_gr <- GRanges(
    seqnames = tx_tbl$seqnames,
    ranges   = IRanges(start = tx_tbl$gene_tss, end = tx_tbl$gene_tss),
    strand   = tx_tbl$strand,
    ensembl_gene_id = tx_tbl$ensembl_gene_id,
    mgi_symbol      = tx_tbl$mgi_symbol
)

gene_tss_canonical <- process_overlapping(teRanges, gene_tss_gr)

df_te_canonical <- data.frame(te_id = gene_tss_canonical$te_hits$te.id,
                    ensembl_gene_id = mcols(gene_tss_canonical$te_hits)$ensembl_gene_id) %>% 
    blackRcloud::splitTEID("te_id") %>%
    filter(order %in% orders.of.interest)

te_stat_results_canonical <- analyze_te_statistics(df_te_canonical, colors = order.color)

# === Find overlaps between TE and gene TSS (±200 bp) ===

gene_tss_win <- resize(gene_tss_gr, width = 201, fix = "center")

overlap_results_win <- process_overlapping(teRanges, gene_tss_win)

df_te_win <- data.frame(te_id = overlap_results_win$te_hits$te.id,
                    ensembl_gene_id = mcols(overlap_results_win$te_hits)$ensembl_gene_id) %>% 
    blackRcloud::splitTEID("te_id")

te_stat_results_win <- analyze_te_statistics(df_te_win, colors = order.color)

