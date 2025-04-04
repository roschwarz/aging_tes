# ---------------------------- Prepare data ------------------------------------
#
# An autonomous TE region is a region that contains at least one CAGE-Peak and
# TE instance that is considered as expressed and was already defined in 
# TE.region.R and stored as a bed file. Currently, the coordinates of these 
# autonomous TE region are stored in ../data/densityPlot/, but that will changed 
# when the scripts are updated.
#
# Load coordinates of autonomously expressed TE regions as a Grange object
# These are all TE regions with at least one CAGE-Peak and RNA-Seq signal
# Subsequently, the autonomous TEs are intersected with genes to identify the
# host genes of autonomous TEs that are located in genes.
#
# Afterwards, the deseq results of CAGE and RNA for TE regions are loaded and
# filtered for the set of autonomous TEs that are located in genes. The same
# is done for the RNA deseq results of genes, whereas this table is filtered 
# for the host genes.
#
# Finally, the filtered tables are merged and can be used for the down stream
# analysis. The merged table is stored as 10_auto_TE_Host_correlation.csv in
# the table directory. The rna_quadrant column is the information in which 
# quadrant the data point is located considering only the correlation between
# the rna signal of the TE region and the host gene.

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    
    source('./01_load_environment.R')
    
} else if ( !ENVIRONMENT_LOADED ) {
    
    source('./01_load_environment.R')
    
}


####################
# Plot Preparation #
####################

###########
# Volcano #
###########

# ======================= Load autonomous TE regions ===========================

auto.TE.cluster.files <- list(brain = "data/shared/brain_downsampled_independent_TE_regions.bed",
                              blood = "data/shared//blood_independent_TE_regions.bed",
                              skin = "data/shared/skinII_independent_TE_regions.bed"
)

auto.TE.regions <- sapply(tissues, simplify = F, function(x){
    df <- read.csv(auto.TE.cluster.files[[x]], 
                   sep = '\t', 
                   header = F) %>% 
        dplyr::pull(V4)
})



# ========================== Filter Count matrix ==============================

dds.TE.regions <- loadRdata("./data/processed/rna_RS_dds_TE_region_SalmonTE.Rdata")

counts.TE.regions <- sapply(tissues, simplify = F, function(x){
   
    counts <- as.data.frame(counts(dds.TE.regions[[x]])) %>% 
        rownames_to_column('te_region_id') %>% 
        filter(te_region_id %in% auto.TE.regions[[x]]) %>% 
        column_to_rownames('te_region_id')
    
})

# ========================== get conditions ====================================

# assign the sample names to a condition (age)
conditions <- sapply(tissues, simplify = F, function(x){
    
    samples <- data.frame(SampleID = names(counts.TE.regions[[x]]),
                          condition = names(counts.TE.regions[[x]])) 
    
    # define the conditions for each sample 
    if (x == 'skin') {
        samples <- samples %>% 
            separate(condition, 
                     c('ID', 'age'), 
                     sep = '[_]') %>%
            mutate(age = str_replace(age, '[0-9]*Skin', '')) %>% 
            mutate(condition = as.factor(case_when(age == 'MO' ~ 'old',
                                                   age == 'MY' ~ 'young'))) %>% 
            dplyr::select(SampleID, condition) %>% 
            column_to_rownames('SampleID')
        
        return(samples)
    }
    
    samples <- samples %>% 
        separate(condition, 
                 c('ID', 'org', 'age', 'tissue', 'number'), 
                 sep = '[_]') %>%
        mutate(condition = as.factor(case_when(age == 'o' ~ 'old',
                                               age == 'y' ~ 'young'))) %>% 
        dplyr::select(SampleID, condition) %>% 
        column_to_rownames('SampleID')
    
    return(samples)
    
})


# ================================ Run DESeq2 ==================================
# runDESeq
dds.te.regions <- sapply(tissues, simplify = F,  function(x){
    
    doDEseq(count.matrix = round(counts.TE.regions[[x]]),
            col_data = conditions[[x]],
            reference = 'young',
            target = 'all'
    )
})

# ============================ Get DESeq results ===============================

deseq.te.deseq <- sapply(tissues, simplify = F, function(x){
    
    getDEseqResults(dds.te.regions[[x]],
                    coefficient = 'condition_old_vs_young',
                    FDR.filter = F)
    
})


# ================================ DESeq Genes =================================

auto_TE_regions_Hosts <- sapply(tissues, simplify = F, function(x){
    
    auto.df <- readGeneric(auto.TE.cluster.files[[x]], strand = 6, meta.cols = list(names = 4))
    
    auto.gene.df <- intersectGranger(geneRanges, auto.df, 'all') %>% 
        dplyr::select(query.ensembl_gene_id, query.external_gene_name, subject.names) %>% 
        dplyr::rename(ensembl_gene_id = query.ensembl_gene_id, 
                      external_gene_name = query.external_gene_name, 
                      te_region_id = subject.names) 
    
    return(auto.gene.df)
    
})

# ========================== Filter Count matrix ==============================

dds.genes <- loadRdata("./data/processed/rna_RS_dds_gene_SalmonTE.Rdata")

counts.host_genes <- sapply(tissues, simplify = F, function(x){
   
    counts <- as.data.frame(counts(dds.genes[[x]])) %>% 
        rownames_to_column('ensembl_gene_id') %>% 
        filter(ensembl_gene_id %in% auto_TE_regions_Hosts[[x]][['ensembl_gene_id']]) %>% 
        column_to_rownames('ensembl_gene_id')
    
})


dds.host_genes <- sapply(tissues, simplify = F,  function(x){
    
    doDEseq(count.matrix = round(counts.host_genes[[x]]),
            col_data = conditions[[x]],
            reference = 'young',
            target = 'all'
    )
})

deseq.host_genes <- sapply(tissues, simplify = F, function(x){
    
    getDEseqResults(dds.host_genes[[x]],
                    coefficient = 'condition_old_vs_young',
                    FDR.filter = F)
    
})

# ============================= Write tables ===================================

deseq.res.merged <- do.call('rbind', sapply(tissues, simplify = F, function(x){
    
    df <- deseq.te.deseq[[x]] %>% 
        rownames_to_column('te_region_id')
    
    df <- merge(df, 
                auto_TE_regions_Hosts[[x]],
                by = 'te_region_id', all.x = T)
    
    df$tissue <- x
    
    return(df)
    
}))

write.table(deseq.res.merged, 
            file = paste0(table_dir, '14_deseq_results_auto_te_regions.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

deseq.res.host.genes.merged <- do.call('rbind', sapply(tissues, simplify = F, function(x){
    
    df <- deseq.host_genes[[x]] %>% rownames_to_column('ensembl_gene_id')
    
    df <- df %>% mutate(sig = ifelse(padj <= FDR, TRUE, FALSE),
                        direction = case_when(log2FoldChange > 0 ~ 'up',
                                              log2FoldChange < 0 ~ 'down',
                                              TRUE ~ 'no'))
    
    df <- merge(df, auto_TE_regions_Hosts[[x]], by = 'ensembl_gene_id', all.x = T) 
    
    df$tissue <- x
    return(df)
}))

write.table(deseq.res.host.genes.merged, 
            file = paste0(table_dir, '14_deseq_results_host_genes.csv'), 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

host_gene_expression_info <- sapply(tissues, simplify = F, function(x){
   
    deseq.res.host.genes.merged %>% 
        filter(tissue == x) %>% 
        dplyr::select(ensembl_gene_id, sig, direction, external_gene_name) %>% 
        filter(!duplicated(ensembl_gene_id))
    
})


###############
# Correlation #
###############

# auto.TE.cluster.files <- list(brain = "data/shared/brain_downsampled_independent_TE_regions.bed",
#                               blood = "data/shared/blood_independent_TE_regions.bed",
#                               skin = "data/shared/skinII_independent_TE_regions.bed"
# )

auto_TE_regions <- sapply(names(auto.TE.cluster.files), simplify = F, function(x){
    
    auto.df <- readGeneric(auto.TE.cluster.files[[x]], strand = 6, 
                           meta.cols = list(names = 4))
    
    auto.gene.df <- intersectGranger(geneRanges, auto.df, 'all') %>% 
        dplyr::select(query.ensembl_gene_id, query.external_gene_name, subject.names) %>% 
        dplyr::rename(ensembl_gene_id = query.ensembl_gene_id, 
                      external_gene_name = query.external_gene_name, 
                      te_region_id = subject.names) 
    
    return(auto.gene.df)
    
})

# ================================= Load TE region rna =========================

te.region.rna <- read.csv(paste0(table_dir, "02_deseq_results_te_region.csv"))

auto_TE_regions_rna <- sapply(tissues, simplify = F, function(x){
    
    df <- te.region.rna %>% filter(tissue == x,
                                   te_region_id %in% auto_TE_regions[[x]][['te_region_id']]) %>% 
        dplyr::select(te_region_id, log2FoldChange, padj, baseMean, pvalue)
    
    names(df) <- c('te_region_id', 'rna_L2FC', 'rna_padj', 'te_region_baseMean', 'te_region_pvalue')
    
    return(df)
})

# ================================= Load TE region cage ========================

te.region.cage <- read.csv(paste0(table_dir, "07_deseq_cage_te_region.csv"))

auto_TE_regions_cage <- sapply(tissues, simplify = F, function(x){
    
    df <- te.region.cage %>% 
        filter(tissue == x, te_region.names %in% auto_TE_regions[[x]][['te_region_id']]) %>% 
        dplyr::select(peak_id, te_region.names, ensembl_gene_id, log2FoldChange, padj, baseMean, pvalue)
    
    names(df) <- c('peak_id', 'te_region_id', 'ensembl_gene_id', 'cage_L2FC', 'cage_padj', 'cage_baseMean', 'cage_pvalue')
    
    return(df)
})

# ================================= Load Gene data =============================


gene.rna <- read.csv(paste0(table_dir, "02_deseq_results_gene.csv"))

genes_with_auto_TE_regions <- sapply(tissues, simplify = FALSE, function(x){
    
    df <- gene.rna %>% 
        filter(tissue == x, ensembl_gene_id %in% auto_TE_regions[[x]][['ensembl_gene_id']]) %>% 
        dplyr::select(ensembl_gene_id, log2FoldChange, padj, external_gene_name, baseMean, pvalue)
    
    names(df) <- c('ensembl_gene_id', 'gene_L2FC', 'gene_padj', 'external_gene_name', 'gene_baseMean', 'gene_pvalue')
    
    return(df)
})

# ================================= Merge Tables ===============================

combined.df <- sapply(tissues, simplify = F, function(x){
    
    df <- merge(auto_TE_regions_cage[[x]], auto_TE_regions_rna[[x]], by = 'te_region_id')
    df <- merge(df, genes_with_auto_TE_regions[[x]], by = 'ensembl_gene_id')
    
    return(df)
    
})

combined.df.merged <- do.call('rbind', sapply(tissues, simplify = F, function(x){
    combined.df[[x]] %>% 
        mutate(tissue = x) %>% 
        mutate(rna_quadrant = ifelse(gene_L2FC > 0 & rna_L2FC > 0, 1,
                                     ifelse(gene_L2FC > 0 & rna_L2FC < 0, 2,
                                            ifelse(gene_L2FC < 0 & rna_L2FC < 0, 3, 4))))
}))

write.table(combined.df.merged,
            file = paste0(table_dir, '14_auto_TE_Host_correlation.csv'),
            col.names = T, 
            row.names = F,
            quote = F,
            sep = ',')

# combined.resis <- combined.df.merged %>% 
#     group_by(tissue) %>% 
#     filter(!duplicated(te_region_id))
# 
# 
# write.table(combined.resis,
#             file = '../../submission/Resis/figure4/figure_4b.csv',
#             col.names = T, 
#             row.names = F,
#             quote = F,
#             sep = ',')



#################
# Plots Panel 4 #
#################

# ---------------------------- Volcano plots -----------------------------------

# ================================ Brain =======================================

data_brain <- deseq.te.deseq$brain %>% 
    rownames_to_column('te_region_id')

data_brain <- merge(data_brain, 
                    auto_TE_regions_Hosts$brain, by = 'te_region_id', all.x = T)

data_brain <- merge(data_brain, host_gene_expression_info$brain, by = c('ensembl_gene_id', 'external_gene_name'), all.x = T)

data_brain$host_gene <- ifelse(is.na(data_brain$ensembl_gene_id), FALSE, TRUE)
data_brain$sig <- ifelse(is.na(data_brain$sig), FALSE, data_brain$sig)

brain.x.max = max(abs(data_brain$log2FoldChange)) + 0.1

df_brain_labels <- data_brain %>% 
    filter(padj <= FDR, host_gene, !duplicated(te_region_id))


brain_no_raster <-  data_brain %>% filter(padj <= FDR)
brain_raster <-  data_brain %>% filter(padj > FDR)

data_brain$tissue <- 'brain'

brain_labels <- c("Fam32a", "Mast4", "Kcnh7", "Dop10", "Nrxn3",
                 "Gm37013", "Rasgrf1", "9330111N05Rik", "Lncpint",
                 "AU020206", "Gm12104")

df_brain_labels <- data_brain %>% 
    filter(padj <= FDR, host_gene, 
           !duplicated(te_region_id), 
           external_gene_name %in% brain_labels)

vol_brain <- ggplot(data_brain, aes(log2FoldChange, -log10(padj))) +
    geom_hline(yintercept = -log10(FDR), linetype = 2) +
    geom_point_rast(data = brain_raster, aes(color = host_gene, shape = sig), alpha = 0.7) + 
    geom_point(data = brain_no_raster, aes(color = host_gene, shape = sig, alpha = 0.7)) +
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


# ================================ skin =======================================

data_skin <- deseq.te.deseq$skin %>% rownames_to_column('te_region_id')

data_skin <- merge(data_skin, auto_TE_regions_Hosts$skin, by = 'te_region_id', all.x = T)

data_skin <- merge(data_skin, host_gene_expression_info$skin, by = c('ensembl_gene_id', 'external_gene_name'), all.x = T)
data_skin$sig <- ifelse(is.na(data_skin$sig), FALSE, data_skin$sig)

data_skin$host_gene <- ifelse(is.na(data_skin$ensembl_gene_id), FALSE, TRUE)

skin.x.max = max(abs(data_skin$log2FoldChange)) + 0.1

df_skin_labels <- data_skin %>% filter(padj <= FDR, host_gene, !duplicated(te_region_id))

skin_no_raster <-  data_skin %>% filter(padj <= FDR)
skin_raster <-  data_skin %>% filter(padj > FDR)

data_skin$tissue <- 'skin'

skin_labels_vol <- c("Skint5", "Skint11", "Rian", "AI506816", "Gja1",
                     "Rnf180", "Chac1", "Gm14703", "Krt27", "Krt35", "Krt28",
                     "Padi3", "Mt4", "Farp2", "Il31ra", "Plaat3", "Krt71", "Krt25",
                     "Krt26", "Calhm4", 'Abi1')

vol_skin <- ggplot(data_skin, aes(log2FoldChange, -log10(padj))) +
    geom_hline(yintercept = -log10(FDR), linetype = 2) +
    geom_point_rast(data = skin_raster, aes(color = host_gene, shape = sig), alpha = 0.7) + 
    geom_point(data = skin_no_raster, aes(color = host_gene, shape = sig, alpha = 0.7)) +
    xlim(c(-skin.x.max, skin.x.max)) +
    labs(y = expression(-log[10](FDR)),
         x = expression(paste(log[2], "(fold ", change[TE_region], ")"))) + 
    facet_grid(.~tissue) + 
    geom_text_repel(data = df_skin_labels %>% filter(external_gene_name %in% skin_labels_vol),
                    aes(label = external_gene_name,
                        segment.size = 0.5,
                        segment.alpha = 0.7), 
                    size = 6/.pt,
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    family = 'Arial',
                    min.segment.length = 0.3,
                    ) +
    scale_color_manual(values = c("FALSE" = 'grey80', "TRUE" = '#192E37')) +
    scale_shape_manual(values = c('FALSE' = 16, 'TRUE' = 8)) +
    theme(legend.position = 'None',
          strip.background = element_rect(fill = tissue.color['skin']),
          strip.text = element_text(size = 10, color = 'black'),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8))



# ================================ blood =======================================

data_blood <- deseq.te.deseq$blood %>% rownames_to_column('te_region_id')

data_blood <- merge(data_blood, auto_TE_regions_Hosts$blood, by='te_region_id', all.x = T)

data_blood <- merge(data_blood, host_gene_expression_info$blood, by = c('ensembl_gene_id', 'external_gene_name'), all.x = T)
data_blood$sig <- ifelse(is.na(data_blood$sig), FALSE, data_blood$sig)

data_blood$host_gene <- ifelse(is.na(data_blood$ensembl_gene_id), FALSE, TRUE)

blood.x.max = max(abs(data_blood %>% filter(!is.na(padj)) %>% pull(log2FoldChange))) + 0.1

df_blood_labels <- data_blood %>% filter(padj <= FDR, host_gene, !duplicated(te_region_id))

blood_no_raster <-  data_blood %>% filter(padj <= FDR)
blood_raster <-  data_blood %>% filter(padj > FDR)

data_blood$tissue <- 'blood'

vol_blood <- 
    ggplot(data_blood, aes(log2FoldChange, -log10(padj))) +
    geom_hline(yintercept = -log10(FDR), linetype = 2) +
    geom_point_rast(data = blood_raster, aes(color = host_gene, shape = sig), alpha = 0.7) + 
    geom_point(data = blood_no_raster, aes(color = host_gene, shape = sig, alpha = 0.7)) +
    xlim(c(-blood.x.max, blood.x.max)) +
    facet_grid(.~tissue) + 
    labs(y = expression(-log[10](FDR)),
         x = expression(paste(log[2], "(fold ", change[TE_region], ")"))) + 
    geom_text_repel(data = df_blood_labels,
                    aes(label = external_gene_name), 
                    size = 6/.pt,
                    box.padding = 0.3,
                    max.overlaps = 50,
                    family = 'Arial') +
    scale_color_manual(values = c("TRUE" = '#192E37', "FALSE" = 'grey80')) +
    scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 16)) +
    theme(strip.background = element_rect(fill = tissue.color['blood']),
          strip.text = element_text(size = 10, color = 'white'),
          legend.position = "None",
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),)


# ggsave(
#     filename = paste0(figure_dir, '14_blood_volcano_autonomous_TE_regions.svg'),
#     plot = last_plot(),
#     width = 10,
#     height = 13,
#     units = "cm",
#     dpi = 300)+
#     theme(legend.position = 'bottom')

# 
# data_brain$tissue <- 'brain'
# data_skin$tissue <- 'skin'
# data_blood$tissue <- 'blood'
# 
# data.table <- rbind(data_brain, rbind(data_skin, data_blood))
# 
# write.table(data.table, 
#             file = '../../submission/Resis/figure4/figure_4a.csv',
#             col.names = T, 
#             row.names = F,
#             quote = F,
#             sep = ',')


# ------------------------------ CORRELATION PLOTS ------------------ 
#
# The tables may contain multiple entries per TE region, as there is a 
# possibility that more than one CAGE-peak may exists per TE region. 
# Consequently, a filtered table is used for the log2FC correlation plots.
# 
# ================================ Brain =======================================


brain_data <- combined.df$brain %>%
    filter(!duplicated(te_region_id))

pcdhb <- gene.rna %>% 
    filter(tissue == 'brain',
           grepl(pattern = '^Pcdhb', x = external_gene_name))


pcdhb.bed <- transGrange(geneRanges) %>% 
    filter(ensembl_gene_id %in% pcdhb$ensembl_gene_id) %>% 
    dplyr::select(seqnames, start ,end, external_gene_name, width, strand) %>% 
    arrange(start)


write.table(pcdhb.bed,
            file = paste0(table_dir, '10_pcdhb.bed'),
            col.names = F, row.names = F, quote = F, 
            sep = '\t')

# handpicked in the genome browser
cluster.of.interest <- c('TE_Cluster_728951', 
                         'TE_Cluster_728962', 
                         'TE_Cluster_728983', #
                         'TE_Cluster_728986', #
                         'TE_Cluster_728991') #

te.region.bed <- transGrange(teRegionRanges) %>% 
    filter(names %in% cluster.of.interest) %>% 
    dplyr::select(seqnames, start ,end, names, width, strand)


write.table(te.region.bed,
            file = paste0(table_dir, '10_pcdhb.te.regions.bed'),
            col.names = F, row.names = F, quote = F, 
            sep = '\t')


df <- system(paste0('bedtools closest -a ', table_dir, '10_pcdhb.te.regions.bed -b ', 
                    table_dir, '10_pcdhb.bed -s'), wait = TRUE, intern = TRUE)

pcdhb.region.pairs <- read.table(text = df) %>% 
    dplyr::select(V4, V10) %>% 
    dplyr::rename(te_region_id = V4,
                  external_gene_name = V10)


pcdhb.region.pairs <- merge(pcdhb.region.pairs, 
                            te.region.rna %>% filter(tissue == 'brain'), 
                            by.x = 'te_region_id', by.y = 'te_region_id') %>% 
    dplyr::select(te_region_id, external_gene_name, log2FoldChange, padj) %>% 
    dplyr::rename(rna_L2FC = log2FoldChange,
                  rna_padj = padj)


pcdhb.region.pairs <- merge(pcdhb.region.pairs,
                            gene.rna %>% filter(tissue == 'brain'),
                            by = 'external_gene_name') %>% 
    mutate(peak_id = NA,
           cage_L2FC = NA,
           cage_padj = NA) %>% 
    dplyr::select(ensembl_gene_id, te_region_id, 
                  peak_id, cage_L2FC, cage_padj,
                  rna_L2FC, rna_padj, log2FoldChange, padj, external_gene_name) %>% 
    dplyr::rename(gene_L2FC = log2FoldChange,
                  gene_padj = padj)

brain_data.pcdhb <- rbind(brain_data %>% 
                              dplyr::select(ensembl_gene_id, te_region_id, 
                                            peak_id, cage_L2FC, cage_padj,
                                            rna_L2FC, rna_padj, gene_L2FC,
                                            gene_padj, external_gene_name),
                          pcdhb.region.pairs)

faded <- brain_data.pcdhb %>% filter(abs(rna_L2FC) < 0.1)
focus <- brain_data.pcdhb %>% filter(abs(rna_L2FC) > 0.1)


brain_data.pcdhb$tissue <- 'brain'


brain_labels_corr <- c("Kcnh7", "Olfr1033", "Pla2g4e", "Banp", "Papola",
                       "Zkscan16", "B230303A05Rik", "Pcdhb11", "Pcdhb9",
                       "Pcdhb15", "Rasgrf1", "Gm37388")

corr_brain <- ggplot(brain_data.pcdhb, aes(rna_L2FC, gene_L2FC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point_rast(data = faded, 
                    alpha = 0.05, 
                    color = '#192E37',
                    size = 0.3) +
    geom_point(data = focus , 
               color = '#192E37',
               size = 0.3) +
    geom_smooth(method = 'lm', 
                color = '#028763',
                linewidth = 0.3) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'Rho',
             p.accuracy = .Machine$double.xmin,
             size = 6/.pt,
             label.sep = '\n',
             label.x = -2, 
             label.y = 1.5,
             family = 'Arial',) +
    xlim(c(-2,2)) +
    ylim(c(-2,2)) +
    geom_text_repel(data = focus %>% filter(external_gene_name %in% brain_labels_corr),
                    aes(label = external_gene_name,
                        segment.size = 0.5,
                        segment.alpha = 0.7), 
                    family = 'Arial',
                    size = 6/.pt,
                    box.padding = 0.4,
                    min.segment.length = 0.2,
                    max.overlaps = 25) +
    stat_quadrant_counts(xintercept = 0, 
                         yintercept = 0,  
                         size = 6/.pt,
                         family = 'Arial',) +
    facet_grid(.~tissue) +
    labs(x = expression(paste(log[2], "(fold ", change[TE_region], ")")),
         y = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          #text = element_text(family = 'arial'),
          strip.background = element_rect(fill = tissue.color['brain']),
          strip.text = element_text(size = 10, color = 'white'))

# =============================== skin =========================================


skin_data <- combined.df$skin %>% 
    filter(!duplicated(te_region_id))

skin_faded <- skin_data %>% filter(abs(rna_L2FC) < 0.3)
skin_focus <- skin_data %>% filter(abs(rna_L2FC) > 0.3)

skin_data$tissue <- 'skin'


skin_labels_corr <- c("Skint5", "Skint11", "Plaat3", "Il31ra", "Chit1",
                      "Pvalb", "C7", "Itih5", "Krt27", "Krt35", "Padi3",
                      "Tpd52l1", "Krt28", "Tenm2", "Arl15", "Mical2", 
                      "Rnf180", "Lap3", "Gja1", "Mt4", "Samd5", "Calhm4")

corr_skin <- ggplot(skin_data, aes(rna_L2FC, gene_L2FC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point_rast(data = skin_faded, 
                    alpha = 0.05,
                    size = 0.3,
                    color = '#192E37') +
    geom_point(data = skin_focus,
               color = '#192E37',
               size = 0.3) +
    geom_smooth(method = 'lm', 
                color = '#028763',
                linewidth = 0.3) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'Rho',
             size = 6/.pt,
             family = "Arial",
             p.accuracy = .Machine$double.xmin,
             label.sep = '\n',
             label.x = -4, 
             label.y = 3) +
    xlim(c(-4,4)) +
    ylim(c(-4,4)) +
    geom_text_repel(data = skin_focus %>% filter(external_gene_name %in% skin_labels_corr), # filter(abs(rna_L2FC) > 1.5),
                    aes(label = external_gene_name,
                        segment.size = 0.5,
                        segment.alpha = 0.7), 
                    size = 6/.pt,
                    family = 'Arial', 
                    box.padding = 0.3,
                    min.segment.length = 0.3,
                    max.overlaps = 50) +
    stat_quadrant_counts(xintercept = 0,
                         yintercept = 0,
                         size = 6/.pt,
                         family = 'Arial') +
    labs(x = expression(paste(log[2], "(fold ", change[TE_region], ")")),
         y = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    facet_grid(.~tissue) +
    #theme_rob() +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          strip.background = element_rect(fill = tissue.color['skin']),
          strip.text = element_text(size = 10, color = 'black'))

# =============================== blood =========================================

blood_data <- combined.df$blood %>% 
    filter(!duplicated(te_region_id), !is.na(rna_padj))

blood_faded <- blood_data %>% filter(abs(rna_L2FC) < 0.2)
blood_focus <- blood_data %>% filter(abs(rna_L2FC) > 0.2)

blood_data$tissue <- 'blood'

corr_blood <- ggplot(blood_data, aes(rna_L2FC, gene_L2FC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point_rast(data = blood_faded,
                    alpha = 0.05, 
                    color = '#192E37',
                    size = 0.3) +
    geom_point(data = blood_focus,
               size = 0.3,
               color = '#192E37') +
    geom_smooth(method = 'lm',
                color = '#028763',
                linewidth = 0.3) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'Rho',
             size = 6/.pt,
             family = 'Arial',
             p.accuracy = .Machine$double.xmin,
             #p.digits = 300,
             label.sep = '\n',
             label.x = -3, 
             label.y = 2) +
    xlim(c(-3,3)) +
    ylim(c(-3,3)) +
    geom_text_repel(data = blood_focus %>% filter(abs(rna_L2FC) > 0.75),
                    aes(label = external_gene_name,
                        segment.size = 0.5,
                        segment.alpha = 0.7), 
                    size = 6/.pt,
                    family = 'Arial',
                    box.padding = 0.3,
                    min.segment.length = 0.3,
                    max.overlaps = 50) +
    facet_grid(.~tissue) +
    stat_quadrant_counts(xintercept = 0,
                         yintercept = 0,
                         size = 6/.pt,
                         family = 'Arial') +
    labs(x = expression(paste(log[2], "(fold ", change[TE_region], ")")),
         y = expression(paste(log[2], "(fold ", change[gene], ")"))) +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          strip.background = element_rect(fill = tissue.color['blood']),
          strip.text = element_text(size = 10, color = 'white'))



######################
# Panel Construction #
######################

panel_4 <- ggarrange(vol_brain, vol_skin, vol_blood,
                     corr_brain, corr_skin, corr_blood, 
                     nrow = 2,
                     ncol = 3)

ggsave(filename = paste0(figure_dir, '14_te_island_panel.pdf'),
       device = cairo_pdf,
       plot = panel_4,
       width = 20,
       height = 16,
       units = 'cm',
       dpi = 300)


