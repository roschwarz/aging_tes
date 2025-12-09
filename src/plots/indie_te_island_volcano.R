# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_plotting_env()
aging_tes::load_annotations()

load_gene_ranges()
load_te_island_annotation()


gene_te_island_association_merged = data.table::fread(file = paste0(table_dir, '02_deseq_results_te_island_extended.csv'))

data_brain <- gene_te_island_association_merged %>% filter(tissue == 'brain', 
                                                           individually
                                                           )

data_brain$host_gene <- ifelse(is.na(data_brain$ensembl_gene_id), FALSE, TRUE)
data_brain$diff <- ifelse(is.na(data_brain$diff), FALSE, data_brain$diff)

brain.x.max = max(abs(data_brain$log2FoldChange)) + 0.1

df_brain_labels <- data_brain %>% 
    filter(padj <= FDR, !duplicated(te_island_id))


brain_no_raster <-  data_brain %>% filter(padj <= FDR)
brain_raster <-  data_brain %>% filter(padj > FDR)

data_brain$tissue <- 'brain'

brain_labels <- c("Fam32a", "Mast4", "Kcnh7", "Dop10", "Nrxn3",
                  "Gm37013", "Rasgrf1", "9330111N05Rik", "Lncpint",
                  "AU020206", "Gm12104")

df_brain_labels <- data_brain %>% 
    filter(padj <= FDR, !duplicated(te_island_id), 
           external_gene_name %in% brain_labels)

vol_brain <- ggplot(data_brain, aes(log2FoldChange, -log10(padj))) +
    geom_hline(yintercept = -log10(FDR), linetype = 2) +
    geom_point_rast(data = brain_raster, aes(color = host_gene, shape = diff), alpha = 0.7) + 
    geom_point(data = brain_no_raster, aes(color = host_gene, shape = diff, alpha = 0.7)) +
    xlim(c(-brain.x.max, brain.x.max)) +
    facet_grid(.~tissue) + 
    labs(y = expression(-log[10](FDR)),
         x = expression(paste(log[2], "(fold ", change[TE_region], ")"))) + 
    geom_text_repel(data = df_brain_labels,
                    aes(label = external_gene_name,
                        segment.size = 0.5,
                        segment.alpha = 0.7), 
                    size = 6/.pt,
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    family = 'Arial',
                    min.segment.length = 0.3) +
    scale_color_manual(values = c("TRUE" = '#192E37', "FALSE" = 'grey80')) +
    scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 16 )) +
    theme(legend.position = 'None',
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          strip.background = element_rect(fill = tissue.color['brain']),
          strip.text = element_text(size = 10, color = 'white'))