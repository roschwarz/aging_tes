

if (!exists("ENVIRONMENT_LOADED")) {
    source("./01_load_environment.R")
} else if (!ENVIRONMENT_LOADED) {
    source("./01_load_environment.R")
}

# ========================= Calc TPMs of TE regions and genes ==================

# -------------------------------- TE  region ----------------------------------

#counts_TE_region <- loadRdata("../data/rna/dds.TE.region.sum.SalmonTE.Rdata")
counts_TE_region <- loadRdata(paste0(data_dir, "rna_RS_dds_TE_region_SalmonTE.Rdata"))


brain_counts <- data.frame(counts(counts_TE_region$brain))

te_region_length <- transGrange(teregionRanges) %>%
    dplyr::select(names, width) %>%
    dplyr::rename(te_region_id = names) %>%
    filter(te_region_id %in% rownames(brain_counts))

counts_normalized <- normalizeCountMatrix(brain_counts, te_region_length)

names(counts_normalized) <- str_replace(names(counts_normalized),
                                        "TPM.",
                                        "age_")

counts_normalized <- rownames_to_column(counts_normalized,
                                        var = "te_region_id")

# -------------------------------- Genes ---------------------------------------

counts_gene <- loadRdata(paste0(data_dir, "rna_RS_dds_gene_SalmonTE.Rdata"))

counts_gene_brain <- data.frame(counts(counts_gene$brain))

gene_length <- transGrange(geneRanges) %>%
    dplyr::select(ensembl_gene_id, width) %>%
    filter(ensembl_gene_id %in% rownames(counts_gene_brain),
           !duplicated(ensembl_gene_id))

counts_gene_brain <- counts_gene_brain %>%
    rownames_to_column("gen.id") %>%
    filter(gen.id %in% gene_length$ensembl_gene_id) %>%
    column_to_rownames("gen.id")


counts_gene_normalized <- normalizeCountMatrix(counts_gene_brain,
                                               gene_length)

names(counts_gene_normalized) <- str_replace(names(counts_gene_normalized),
                                             "TPM.",
                                             "age_")

counts_gene_normalized <- rownames_to_column(counts_gene_normalized,
                                             var = "ensembl_gene_id")


# ====================== Get TE regions and genes of interest ==================

# handpicked in the genome browser for figure 5 in the publication
cluster.of.interest <- c("TE_Cluster_728951",
                         "TE_Cluster_728962",
                         "TE_Cluster_728983",
                         "TE_Cluster_728986",
                         "TE_Cluster_728991")

Pcdhb1.coordinates <- transGrange(geneRanges) %>%
    filter(external_gene_name == "Pcdhb1")

Pcdhb22.coordinates <- transGrange(geneRanges) %>%
    filter(external_gene_name == "Pcdhb22")

genes <- transGrange(geneRanges) %>%
    filter(seqnames == "chr18",
           start >= Pcdhb1.coordinates$start,
           end <= Pcdhb22.coordinates$end,
           strand == "+")


# ==================== Filter count tables and calc z scores ===================

# --------------------------- TE region ----------------------------------------

counts_normalized_pro <-  counts_normalized %>%
    filter(te_region_id %in% cluster.of.interest) %>%
    gather(key = "sample",
           value = "TPM",
           names(counts_normalized)[2:length(counts_normalized)])


counts_normalized_pro <- counts_normalized_pro %>%
    separate("sample",
             c("age.name", "age", "sample.id", "tissue"),
             sep = "[_]",
             remove = FALSE) %>%
    dplyr::select(-age.name)


counts_normalized_pro <- counts_normalized_pro %>%
    dplyr::select(te_region_id, sample, TPM) %>%
    spread(sample, TPM) %>%
    column_to_rownames(var = "te_region_id")

# --------------------------- Genes ----------------------------------------

counts_gene_norm <- counts_gene_normalized %>%
    filter(ensembl_gene_id %in% genes$ensembl_gene_id)  %>%
    gather(key = "sample",
           value = "TPM",
           names(counts_gene_normalized)[2:length(counts_gene_normalized)])

counts_gene_norm <- counts_gene_norm %>%
    separate("sample",
             c("age.name", "age", "sample.id", "tissue"),
             sep = "[_]",
             remove = FALSE) %>%
    dplyr::select(-age.name)

counts_gene_norm <- counts_gene_norm %>%
    dplyr::select(ensembl_gene_id, sample, TPM) %>%
    spread(sample, TPM) #%>%
#column_to_rownames(var ='ensembl_gene_id')

counts_gene_norm <- merge(counts_gene_norm,
                          transGrange(geneRanges) %>%
                              dplyr::select(ensembl_gene_id,
                                            external_gene_name),
                          by = "ensembl_gene_id") %>%
    dplyr::select(-ensembl_gene_id) %>%
    column_to_rownames("external_gene_name")


# --------------------------- Merge and sort tables ----------------------------

norm.counts <- rbind(counts_normalized_pro, counts_gene_norm)

order.of.elements <- rbind(
    
    genes %>%
        dplyr::select(-ensembl_gene_id) %>%
        dplyr::rename(names = external_gene_name),
    
    transGrange(teregionRanges) %>%
        filter(names %in% cluster.of.interest)
    
)

norm.counts <- norm.counts[order.of.elements[order(order.of.elements$start),
                                             "names"], ]


norm.counts <- t(apply(as.matrix(norm.counts),
                       1,
                       function(x) calc_z_score(x, capped = 2)))

norm.counts.resis <- as.data.frame(norm.counts) %>%
    rownames_to_column("feature_id")


write.table(norm.counts.resis,
            file = "./manuscripts/nature_aging/resis/figures/figure5/5B_pcdhb_heatmap.csv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = ",")

# ======================== get quantification numbers ##########################

gene.deseq <- read.csv(paste0(table_dir, "02_deseq_results_gene.csv")) %>%
    filter(external_gene_name %in% genes$external_gene_name,
           tissue == "brain") %>%
    dplyr::select(external_gene_name, log2FoldChange, padj) %>% 
    dplyr::rename(id = external_gene_name)

#gene.deseq %>% filter(external_gene_name == "Pcdhb1")

te.region.deseq <- read.csv(paste0(table_dir, "02_deseq_results_te_region.csv")) %>%
    filter(te_region_id %in% cluster.of.interest, tissue == "brain") %>%
    dplyr::select(te_region_id, log2FoldChange, padj) %>%
    dplyr::rename(id = te_region_id)


df <- rbind(gene.deseq, te.region.deseq)

df <- df[match(rownames(norm.counts), df$id),]

#col_fun_l2fc = colorRamp2(c(-2,0,2), c('orange','white','darkgreen'))

col_fun_l2fc = colorRamp2(c(-2,-1,0,1,2), c('#264653','#2a9d8f', 'white','#f4a261', '#e76f51'))
col_fun_fdr = colorRamp2(c(0.000001,0.00001,0.0001, 0.001, 0.01, 0.1, 1), 
                         c('#AC0A0A','#B82C2C','#C44F4F','#D17272','#DD9595', '#E9B8B8', '#F6DBDB'))

col_fun_fdr = colorRamp2(c(0.00000001,0.0000001,0.000001,
                           0.00001,0.0001,0.001, 
                           0.05, 0.1, 0.5,
                           1), 
                         c("#ea698b","#d55d92","#c05299","#ac46a1","#973aa8","#822faf","#6d23b6","#6411ad","#571089","#47126b"))

row_legend = list(direction = "horizontal",
                  legend_width = unit(2, 'cm'),
                  grid_height = unit(0.1, 'cm'),
                  title_gp = gpar(fontsize = 6),
                  labels_gp = gpar(fontsize = 6),
                  title_position = "topcenter")

row_ha = rowAnnotation(L2FC = df$log2FoldChange,
                       FDR = df$padj,
                       col = list(L2FC = col_fun_l2fc,
                                  FDR = col_fun_fdr),
                       simple_anno_size = unit(0.3, 'cm'),
                       annotation_name_gp = gpar(fontsize = c(6,6)),
                       annotation_name_side = c('top', 'top'),
                       annotation_legend_param = list(L2FC = row_legend, 
                                                      FDR = row_legend))



#colnames(norm.counts) <- NULL

ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points") 
col_fun = colorRamp2(c(-1.5, 0, 1.5), c('#457b9d','white','#e63946'))

row_text_size = 6

pcdhb_heatmap <- Heatmap(norm.counts,
        col = col_fun ,
        column_title = "Protocadherin beta cluster",
        column_title_gp = gpar(color = 'white', border = 'black', fontsize = 8),
        border = TRUE,
        column_gap = unit(0, 'mm'),
        row_names_gp = gpar(fontsize = row_text_size),
        column_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = TRUE,
        row_names_side = "left",
        cluster_columns = FALSE,
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_title_gp = gpar(fontsize = 6),
        top_annotation = HeatmapAnnotation(age = anno_block(gp = gpar(fill = c('white', 'white')),
                                                            labels = c('young', 'old'),
                                                            labels_gp = gpar(col = 'black', fontsize = 8))),
        right_annotation = row_ha,
        column_km = 2,
        cluster_rows = FALSE,
        column_labels = rep('', 10),
        heatmap_legend_param = list(title = "z-score of TPM", 
                                    direction = 'horizontal',
                                    legend_width = unit(2, 'cm'),
                                    grid_height = unit(0.1, 'cm'),
                                    title_gp = gpar(fontsize = 6),
                                    labels_gp = gpar(fontsize = 6),
                                    title_position = "topcenter"))


pdf(file = paste0(figure_dir, '20_heatmap_pcdhb.pdf'), width = 3.2, height = 4.1)
draw(pcdhb_heatmap,  heatmap_legend_side = "bottom")
dev.off()



pl <- pheatmap(norm.counts,
               fontsize = 8,
               #fontfamily = "Arial",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               #breaks = seq(-2, 2, length.out = 100),
               color = colorRampPalette(c("#457b9d", "white", "#e63946"))(100),
               #cutree_rows = 2,
               #cutree_cols = 2,
               treeheight_row = 0, # removes dendrogram
               treeheight_col = 0
)


ggsave(
    filename = paste0(figure_dir, "20_heatmap_protocadherins.svg"),
    plot = pl,
    width = 10,
    height = 12,
    units = "cm",
    dpi = 300)




