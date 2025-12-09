# ================================== Setup ================================================

# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

load_cage_peak_annotation()
load_te_ranges()
load_gene_ranges()

tss_tx_gr <- load_gene_tss_gr()

mcols(teRanges)$te_class <- data.frame(teRanges) %>% blackRcloud::splitTEID("te.id") %>% pull(order)

# -----------------------------------------------------------------------------------------------------------
# Pie Charts
# -----------------------------------------------------------------------------------------------------------


tss_covered_by_cage <- do.call('rbind', sapply(tissues, simplify = FALSE, function(t){
    
    tss_covered_by_cage <- blackRcloud::intersectGranger(tss_tx_gr, cageRanges[[t]], tab = 'all')
    n_tss_covered_by_cage <- length(unique(tss_covered_by_cage$query.ensembl_transcript_id))
    
    # Nenner: alle CAGE-Peaks
    n_total <- length(tss_tx_gr)
    n_not   <- n_total - n_tss_covered_by_cage
    
    df <- data.frame(
        status = c("TSS covered by CAGE", "TSS wo CAGE"),
        n      = c(n_tss_covered_by_cage, n_not),
        tissue = t
    )
    
    return(df)
    
}))

tss_covered_by_cage <- tss_covered_by_cage %>% 
    group_by(tissue) %>% 
    mutate(pct = n/sum(n),
           label = paste0(n, " (", label_percent(accuracy = 0.1)(pct), ")" )
    )

tss_covered_by_cage$tissue <- factor(tss_covered_by_cage$tissue, levels = tissues)

tss_cov_pl <- ggplot(tss_covered_by_cage, aes(x = "", y = n, fill = status)) +
    geom_col(width = 1, color = "white") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 6/.pt) + 
    coord_polar(theta = "y", clip = 'off') +
    scale_fill_manual(values = c("TSS covered by CAGE" = "steelblue", "TSS wo CAGE" = "lightgray")) +
    #labs(title = "Known TSS covered by CAGE-peaks", x = NULL, y = NULL) +
    facet_grid(tissue~.) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.position = 'bottom',
          aspect.ratio = 1,
          plot.margin = margin(0,0,0,0, 'pt'),
          panel.spacing = unit(0, 'pt'),
          legend.key.size = unit(6, 'pt'),
          legend.text = element_text(size = 6))


prepare_complex_pi <- function(cage_peaks_gr, exonRanges, geneRanges, teRanges) {
    # prepare data frame
    # Rules: TE > Exon > Intron > Intergen
    cat <- rep("intergenic", length(cage_peaks_gr))
    
    idx_te <- which(overlapsAny(cage_peaks_gr, teRanges, ignore.strand = FALSE))
    cat[idx_te] <- "TE"
    
    idx_remaining <- setdiff(seq_along(cage_peaks_gr), idx_te)
    idx_exon <- which(overlapsAny(cage_peaks_gr[idx_remaining], gene.exonRanges, ignore.strand = FALSE))
    cat[idx_exon] <- "exon"
    
    idx_remaining <- setdiff(seq_along(idx_remaining), idx_exon)
    idx_intron <- which(overlapsAny(cage_peaks_gr[idx_remaining], geneRanges, ignore.strand = FALSE))
    cat[idx_intron] <- "intron"
    
    # count and percent labeling
    df <- as.data.frame(table(category = factor(cat, levels = c("TE","intron","exon","intergenic"))))
    
    
    names(df) <- c("category","n")
    df$pct <- df$n / sum(df$n)
    df$label <- paste0(df$category, ": ", df$n, " (", label_percent(accuracy = 0.1)(df$pct), ")")
    return(df)
}

complex_pi_df <- do.call('rbind', sapply(tissues, simplify = FALSE, function(t) {
        peaks <- cageRanges[[t]]
        tss_peaks <- which(overlapsAny(peaks, tss_tx_gr, ignore.strand = FALSE))
        non_tss_peaks <- setdiff(seq_along(peaks), tss_peaks)
    
        cage_peaks <- peaks[non_tss_peaks]
        
        n_peaks <- length(cageRanges[[t]])
        
        df <- prepare_complex_pi(cage_peaks,
                                 gene.exonRanges, 
                                 geneRanges,
                                 teRanges)
        df$tissue <- t
        
        print(paste0("n_cage: ", n_peaks))
        print(paste0("n_cage filtered: ", n_peaks))
        
        return(df)
    
    }))

complex_pi_df$tissue <- factor(complex_pi_df$tissue, levels = tissues)

# 4) Torten-Plot mit Prozenten
complex_pi_pl <- ggplot(complex_pi_df, aes(x = 1, y = n, fill = category)) +
    geom_col(width = 1, color = "white") +
    #geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 8/.pt) +
    geom_label_repel(aes(label = label),
                     #position = position_stack(vjust = 0.5),
                     size = 6/.pt,
                     show.legend = FALSE,
                     nudge_x = 0.7,
                     segment.size = 0.2,
                     color = 'black',
                     segment.color = 'grey50',
                     max.overlaps = Inf,
                     box.padding = 0.3,
                     force = 5,
    ) +
    coord_polar(theta = "y", clip = 'off') +
    scale_fill_manual(values = c("TE" = "#8b8c89", "intron" = "#d5bdaf", "exon" = "#6c584c", "intergenic" = "lightgray")) +
    facet_grid(rows = vars(tissue), scales = 'free') +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.position = 'bottom',
          aspect.ratio = 1,
          plot.margin = margin(0,0,0,0, 'pt'),
          panel.spacing = unit(0, 'pt'),
          legend.key.size = unit(6, 'pt'),
          legend.text = element_text(size = 6),
          panel.spacing.y = unit(0, 'pt')
          )

pl <- ggarrange(tss_cov_pl, complex_pi_pl, 
          col = 2,
          padding = 0)

# Save figure with metadata to index
meta <- list(name = 'tss_cage_coverage_pie_charts',
             description = 'Pie charts showing the coverage of known TSS by CAGE peaks and the genomic context of CAGE peaks not overlapping known TSS.',
             tags = c('tss', 'cage-seq', 'genomic context', 'pie chart'),
             parameters = list(tissues = c('brain', 'skin', 'blood'), sex = c('male')),
             script = 'characterizing_TSS.R'
)

fig_index(plot = pl,
          outdir = figure_dir,
          meta = meta,
          index_file = 'figure_index.tsv',
          width = 16,
          height = 10,
          dpi = 300,
          format = 'pdf')

# ------------------------------------------------------------------------------------------------------------
# Sankey Plot
# ------------------------------------------------------------------------------------------------------------
# Prepare data for Sankey plot

library(ggalluvial)
prepare_sankey_table <- function(cageRanges, exonRanges, geneRanges, teRanges){
    
    # 1. Assign exon overlaps first
    overlaps_exon <- logical(length(cageRanges))
    exon_hits <- findOverlaps(cageRanges, exonRanges)
    overlaps_exon[unique(queryHits(exon_hits))] <- TRUE
    
    # 2. Assign gene overlaps
    overlaps_gene <- logical(length(cageRanges))
    gene_hits <- findOverlaps(cageRanges, geneRanges)
    overlaps_gene[unique(queryHits(gene_hits))] <- TRUE
    
    # 3. Determine annotation:
    # exon if overlaps exon
    # intron if overlaps gene but NOT exon
    # intergenic otherwise
    annotation <- rep("intergenic", length(cageRanges))
    annotation[overlaps_gene] <- "intron"
    annotation[overlaps_exon] <- "exon"
    
    # 4. Find overlaps with TEs and TE classes
    te_hits <- findOverlaps(cageRanges, teRanges)
    has_te <- rep(FALSE, length(cageRanges))
    te_class <- rep(NA_character_, length(cageRanges))
    has_te[queryHits(te_hits)] <- TRUE
    te_class[queryHits(te_hits)] <- mcols(teRanges)$te_class[subjectHits(te_hits)]
    
    # 5. Set target to TE if overlapping TE, else annotation
    target <- ifelse(has_te, "TE", annotation)
    
    # 6. Summarize counts for Sankey plot
    sankey_df <- data.frame(
        source = annotation,
        target = target,
        te_class = ifelse(has_te, te_class, NA_character_)) %>% 
    mutate(te_class = str_replace(te_class, "\\?", "")) %>% 
        group_by(source, target, te_class) %>%
        summarize(count = dplyr::n(), .groups = 'drop') %>% 
        arrange(factor(source, levels = c('exon', 'intron', 'intergenic')), target, te_class)
    
    return(sankey_df)
}

prepare_sankey_stop_flows <- function(sankey_df){
    
    # sankey_df contains source, target, te_class, count
    
    # Non-TE flows: Create a third axis same as target (to stop flows)
    non_te <- sankey_df %>%
        filter(target != "TE") %>%
        mutate(axis3 = target)
    
    # TE flows: third axis = te_class
    te_flows <- sankey_df %>%
        filter(target == "TE") %>%
        mutate(axis3 = te_class)
    
    # Combine flows
    combined <- bind_rows(non_te, te_flows) %>%
        arrange(factor(source, levels = c("exon", "intron", "intergenic")),
                target, axis3)
    
    return(combined)
    
    
} 

df_sankey <- prepare_sankey_table(cageRanges$brain, gene.exonRanges, geneRanges, teRanges)

# Plotting

ggplot(df_sankey, aes(axis1 = source, axis2 = target, y = count)) +
    geom_alluvium(aes(fill = te_class), width = 1/18, alpha = 0.7) +
    geom_stratum(aes(fill = te_class), width = 1/18, color = "grey") +
    scale_fill_manual(values = order.color) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Genomic Context", "TE Overlap")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = NULL,
         y = "Number of CAGE Peaks")

