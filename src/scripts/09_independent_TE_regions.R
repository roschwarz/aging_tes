# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}

rna.deseq.res.file <- 'data/processed/rna_RS_deseq_TE_instances_SalmonTE.Rdata'
#rna.deseq.res.file <- '/misc/paras/data/rschwarz/workspace/phd/projects/mCQuaRna/packages/publication_ready/data/rna/deseq.TE.SalmonTE.Rdata'


# cage.files <- list(
#     brain = '../../results/cage_RS/brain_downsampled_gclipped/counts/brain_downsampled_peak_counts.csv',
#     blood = '../../results/cage_RS/blood_gclipped/counts/blood_peak_counts.csv',
#     skin = '../../results/cage_RS/skinII_gclipped/counts/skinII_peak_counts.csv')
# 
# quant.files <- list(
#     brain = 'results/quant_RS/brain/counts/brain_peak_counts.csv',
#     blood = 'results/quant_RS/blood/counts/blood_peak_counts.csv',
#     skin = 'results/quant_RS/skinII/counts/skinII_peak_counts.csv')
# 
# tissues = c('brain', 'skin', 'blood')
# 
# teRanges <- getRanges("mmusculus", 102, "te", 1)

# featureCounts was used to count reads that mapped to the called cage peaks (by peakachu)
# The matrix from featureCounts does contain the coordinates of the called peaks. These
# coordinates are used to create a GRanges object of the called CAGE-peaks
cageRanges <- sapply(tissues, simplify = F, function(tissue){
    
    directory <- strsplit(cage.files[[tissue]], "/")[[1]][3]
    
    file <- paste0('results/cage_RS/',
                   directory,
                   '/raw_peaks/',
                   tissue,
                   '.cage_peaks.bed')
    
    return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
    
})


quantRanges <- sapply(tissues, simplify = F, function(tissue){
    
    directory <- strsplit(quant.files[[tissue]], "/")[[1]][3]
    
    file <- paste0('results/quant_RS/',
                   directory,
                   '/raw_peaks/',
                   tissue,
                   '.quant_peaks.bed')
    
    return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
    
})


#create and load TE region annotation
#te.region.file = '../../data/shared/mm10_TE_region_extended.bed'
# te.region.file = '../../data/shared/mm10_TE_region_5prime_extended.bed'
# 
# TE.region <- readGeneric(te.region.file, strand = 6, meta.cols = list(names = 4)) 

# Collect TE instances of TE-regions to get information about the composition of the TE
# regions later in the analysis.
#TE.region.instances <- intersectGranger(TE.region, teRanges, tab = 'all')
TE.region.instances <- intersectGranger(TEregionRanges_5prime, teRanges, tab = 'all')

# TEs can intersect with multiple genes, so that they occur multiple times in
# the teRange objects. Therefore a unique is applied at the end of the pipeline
TE.region.instances <- TE.region.instances %>% 
    dplyr::select(query.names, subject.te.id, subject.position) %>% 
    filter(!duplicated(subject.te.id)) %>% 
    unique()

names(TE.region.instances) <- c("te_region_id", "te_id", "te_position")

#########################
# autonomous TE regions #
#########################
# Usage of CAGE- and RNA-seq to get autonomously expressed TE regions. 
#  1. Intersect the extended TE regions with CAGE-peaks to get potential autonomous TE region.
#  2. Get autonomous TE regions that also contain at least one expressed TE that 
#       is detected via RNA-Seq. An element is considered as expressed when a 
#       adjusted p-value after the DESeq analysis is available.


# Cage intersection
auto.TE.region <- sapply(names(cageRanges), 
                         simplify = F, 
                         USE.NAMES = T, 
                         function(x){
                             
                             intersectGranger(cageRanges[[x]], TEregionRanges_5prime, tab = 'all')
                             
                         })


# RNA intersection
rna.deseq.res <- loadRdata(rna.deseq.res.file)

expressed.TE.regions <- sapply(names(rna.deseq.res), simplify = F, function(x){
    
    expressed.TEs.rna <- rownames(rna.deseq.res[[x]]  %>% filter(!is.na(padj)))
    
    exp.regions <- TE.region.instances %>% 
        filter(te_id %in% expressed.TEs.rna) %>% 
        dplyr::pull('te_region_id')
    
    # regions Ids occur multiple times because of the composition of the regions
    # However, only one id is needed to get the expressed TE regions 
    return(unique(exp.regions))
})



# filter the cage intersecting TE regions for those that intersect with 
# expressed TEs detected via RNA-Seq
auto.TE.region <- sapply(names(auto.TE.region), simplify = F, function(x){
    
    region <- auto.TE.region[[x]]
    exp.te <- expressed.TE.regions[[x]]
    
    region %>% 
        filter(subject.names %in% exp.te) %>% 
        filter(!duplicated(subject.names))
    
})

## Write bed-files of autonomously expressed TE clusters
sapply(names(auto.TE.region), function(x){
    
    df <- auto.TE.region[[x]]
    
    bed.file <- GRanges(seqnames = df$subject.seqnames,
                        ranges = IRanges(start = df$subject.start,
                                         end = df$subject.end,
                                         names = df$subject.names),
                        strand = df$subject.strand,
                        cluster.ID = df$subject.names,
                        value = 1)
    
    bed.file <- transGrange(bed.file)
    
    bed.file <- bed.file[c('seqnames', 'start', 'end', 'cluster.ID',  'value', 'strand')]
    
     
    if (x == 'brain') {
        
        filename <- paste0('data/shared/', x,
                           '_downsampled_independent_TE_regions.bed')
        
    }else if (x == 'skin') {
        
        filename <- paste0('data/shared/', x,
                           'II_independent_TE_regions.bed')
    }else{
        
        filename <- paste0('data/shared/', x,
                           '_independent_TE_regions.bed')
    }
    
    write.table(bed.file, file = filename, sep = '\t', col.names = F, row.names = F, quote = F)
    
})

## Write bedGraph files of the CAGE peaks
sapply(names(cageRanges), function(x){
    
    df <- transGrange(cageRanges[[x]]) %>% dplyr::select(seqnames, start, end)
    
    df$value <- 1
    if (x == 'brain') {
        
        filename <- paste0('results/cage_RS/', x,
                           '_downsampled_gclipped/coverage/',
                           x, '_downsampled.cage_peaks.one.bedgraph')
        
    }else if (x == 'skin') {
        
        filename <- paste0('results/cage_RS/', x,
                           'II_gclipped/coverage/', x,
                           'II.cage_peaks.one.bedgraph')
    }else{
        
        filename <- paste0('results/cage_RS/', x,
                           '_gclipped/coverage/', x, '.cage_peaks.one.bedgraph')
    }
    
    write.table(df,
                file = filename, sep = '\t', col.names = F, row.names = F, quote = F)
})


## Write bedGraph files of the Quant peaks
sapply(names(quantRanges), function(x){
    
    df <- transGrange(quantRanges[[x]]) %>% dplyr::select(seqnames, start, end)
    
    df$value <- 1
    if (x == 'skin') {
        
        filename <- paste0('results/quant_RS/', x,
                           'II/coverage/', x,
                           'II.quant_peaks.one.bedgraph')
    }else{
        
        filename <- paste0('results/quant_RS/', x,
                           '/coverage/', x, '.quant_peaks.one.bedgraph')
    }
    
    write.table(df,
                file = filename, sep = '\t', col.names = F, row.names = F, quote = F)
})


#############################
# Composition of TE regions #
#############################
# Take the TE region instance table and filter for all TE regions that are 
# autonomously expressed. Determine the position of the instance within the TE
# region and categorize the TE regions into single (consists of one instance),
# double (consists of two instances) and multiple (consists of more than 2 instance)

auto.TE.region.composition <- sapply(names(auto.TE.region), simplify = F, function(x){
    
    auto.TE.region.instances <- TE.region.instances %>% 
        filter(te_region_id %in% auto.TE.region[[x]][['subject.names']])
    
    
    names(auto.TE.region.instances) <- c('te_region_id', 'te_id', 'te_genomic_pos')
    
    auto.TE.region.instances <- splitTEID(auto.TE.region.instances, 'te_id')
    
    auto.TE.region.instances <- auto.TE.region.instances %>% 
        group_by(te_region_id) %>% 
        mutate(member = n(),
               position = ifelse(strand == '-', rev(1:n()), 1:n()),
               region.type = ifelse(member == 1, 'single', 
                                    ifelse(member == 2, 'double', 'multiple')))
    
    
    return(auto.TE.region.instances)
})



auto.region.types <- sapply(c('brain', 'skin', 'blood'), simplify = F, function(x){
    
    df <- auto.TE.region.composition[[x]]
    
    df <- df %>% filter(!duplicated(te_region_id))
    
    df$region.type <- factor(df$region.type, levels = c('single', 'double', 'multiple'))
    
    print(paste(x, ":", nrow(df)))
    
    pl <- ggplot(df, aes(region.type)) +
        geom_bar() +
        labs(title = x) +
        theme(axis.text.x = element_text(angle = 90))
    
    return(pl)
    
})

gridExtra::grid.arrange(grobs = auto.region.types, nrow = 1)

auto.single.region <- sapply(names(auto.TE.region.composition), simplify = F, function(x){
    
    df <- auto.TE.region.composition[[x]]
    
    df <- df %>%
        filter(member == 1) %>%
        order.TEs('super_family', decreasing = F)
    
    pl <- ggplot(df, aes(super_family)) +
        geom_bar() +
        coord_flip()
    
    return(pl)
})

gridExtra::grid.arrange(grobs = auto.single.region, nrow = 1)

auto.region.super.composition <- sapply(names(auto.TE.region.composition), simplify = F, function(x){
    
    df <- auto.TE.region.composition[[x]]
    
    df <- df %>% 
        filter(member > 2) %>% 
        mutate(super_family = case_when(super_family == "Alu" ~ "B1", .default = super_family)) %>% 
        group_by(te_region_id) %>% 
        mutate(instance.position = ifelse(position == 1, 'first', ifelse(position == max(position), 'last', 'body'))) %>% 
        ungroup() %>% 
        dplyr::count(super_family, instance.position) %>% 
        group_by(instance.position) %>% 
        mutate(percent = n/sum(n))
    
    df$instance.position <- factor(df$instance.position, levels = c('first', 'body', 'last'))
    
    
    df <- df %>% 
        dplyr::select(-n) %>% 
        filter(super_family != "hAT-Charlie") %>% 
        spread(key = instance.position, value = percent, fill = 0) %>% 
        column_to_rownames(var = 'super_family') %>% 
        filter(rowSums(.) > 0.05)
    
    return(as.matrix(df))
}) 
    
ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points")    # setting to get the titels in the heatmaps centered from a horizontal perspective see: https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html?q=color%20title#heatmap-titles
col_fun = colorRamp2(c(0, 0.4), c('white','#e63946'))

pl_brain <-  Heatmap(auto.region.super.composition$brain,
                     col = col_fun,
                     border = TRUE,
                     column_gap = unit(0, 'mm'),
                     column_title_gp = gpar(fill = tissue.color['brain'], color = 'black', border = 'black', fontsize = 8), 
                     column_title = 'brain',
                     row_names_side = 'left',
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     show_heatmap_legend = FALSE,
                     cluster_columns = FALSE,
                     column_names_rot = 0,
                     column_names_centered = TRUE,
                     cluster_rows = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.2f", auto.region.super.composition$brain[i, j]), x, y, gp = gpar(fontsize = 6))
                     },
                     heatmap_legend_param = list(title = "Proportion of TE-families", 
                                                 direction = 'horizontal',
                                                 legend_width = unit(3, 'cm'))
)

pl_skin <-  Heatmap(auto.region.super.composition$skin,
                     col = col_fun,
                     border = TRUE,
                     column_title = "skin",
                     column_title_gp = gpar(fill = tissue.color['skin'], color = 'black', border = 'black', fontsize = 8), 
                     column_gap = unit(0, 'mm'),
                     row_names_side = 'left',
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     show_heatmap_legend = TRUE,
                     cluster_columns = FALSE,
                     column_names_rot = 0,
                     column_names_centered = TRUE,
                     cluster_rows = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.2f", auto.region.super.composition$skin[i, j]), x, y, gp = gpar(fontsize = 6))
                     },
                     heatmap_legend_param = list(title = "Proportion of TE-families", 
                                                 direction = 'horizontal',
                                                 legend_width = unit(3, 'cm'),
                                                 grid_height = unit(0.2, 'cm'),
                                                 title_gp = gpar(fontsize = 6),
                                                 labels_gp = gpar(fontsize = 6),
                                                 title_position = "topcenter"))

pl_blood <-  Heatmap(auto.region.super.composition$blood,
                     col = col_fun,
                     column_title = "blood",
                     column_title_gp = gpar(fill = tissue.color['blood'], color = 'white', border = 'black', fontsize = 8), 
                     border = TRUE,
                     column_gap = unit(0, 'mm'),
                     row_names_side = 'left',
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     show_heatmap_legend = FALSE,
                     cluster_columns = FALSE,
                     column_names_rot = 0,
                     column_names_centered = TRUE,
                     cluster_rows = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.2f", auto.region.super.composition$blood[i, j]), x, y, gp = gpar(fontsize = 6))
                     },
                     heatmap_legend_param = list(title = "Proportion of TE-families", 
                                                 direction = 'horizontal',
                                                 legend_width = unit(3, 'cm'),
                                                 legend_heigt = unit(0.25, 'cm'))
)


ht_list <- pl_brain + pl_skin + pl_blood

pdf(file = paste0(figure_dir, '09_superFam_composition_te_island_3.7x2.5_300.pdf'), width = 3.78, height = 2.8)
    draw(ht_list, heatmap_legend_side = "bottom")
dev.off()
    


#########
# RESIS #
#########


pdf(file = "./manuscripts/nature_aging/resis/figures/figure3/3B_TE_island_composition.pdf", width = 3.78, height = 2.8)
    draw(ht_list, heatmap_legend_side = "bottom")
dev.off()


TE_composition_df <- do.call('rbind', sapply(tissues, simplify = F, function(tissue){
    
    df <- as.data.frame(auto.region.super.composition[[tissue]])
    
    df <- df %>% rownames_to_column(var = 'super_family')
    
    df$tissue <- tissue
    
    return(df)
})
)

write.csv(TE_composition_df, file = "./manuscripts/nature_aging/resis/figures/figure3/3B_TE_island_composition.csv", row.names = F)



# ========================== UNDER CONSTRUCTION ================================

# Quant intersection
# quant.TE.region <- sapply(names(quantRanges), 
#                          simplify = F, 
#                          USE.NAMES = T, 
#                          function(x){
#                              
#                              intersectGranger(quantRanges[[x]], TE.region, tab = 'all')
#                              
#                          })
# 
# # filter the cage intersecting TE regions for those that intersect with 
# # expressed TEs detected via RNA-Seq
# auto.TE.quant.region <- sapply(names(auto.TE.region), simplify = F, function(x){
#     
#     region <- auto.TE.region[[x]]
#     quant.te <- quant.TE.region[[x]]
#     
#     region %>% 
#         filter(subject.names %in% quant.te[['subject.names']]) %>% 
#         filter(!duplicated(subject.names))
#     
# })
# 
# 
# brain.auto.TE.region <- auto.TE.region$brain
# 
# names(brain.auto.TE.region) <- str_replace_all(names(brain.auto.TE.region), c('query' = 'cage', 'subject' = 'te_region'))
# 
# brain.quant.TE.region <- quant.TE.region$brain
# names(brain.quant.TE.region) <- str_replace_all(names(brain.quant.TE.region), c('query' = 'quant', 'subject' = 'te_region'))
# 
# brain.quant.TE.region <- brain.quant.TE.region %>% dplyr::select(-c("te_region.seqnames", 
#                                            "te_region.start",
#                                            "te_region.end", 
#                                            "te_region.width",
#                                            "te_region.strand"))
# 
# te.islands <- merge(brain.auto.TE.region, brain.quant.TE.region, by = 'te_region.names')
# 
# te.islands.bed <- te.islands %>% 
#     mutate(te_island.length = case_when(te_region.strand == "+" ~ (quant.end - cage.start),
#                                         te_region.strand == "-" ~ (quant.start - cage.end))) %>% 
#     filter(te_island.length > 0) %>% 
#     group_by(te_region.names) %>%
#     filter(te_island.length == max(te_island.length)) %>% 
#     ungroup() %>% 
#     mutate(chr = cage.seqnames,
#            start = case_when(te_region.strand == "+" ~ cage.start,
#                              te_region.strand == "-" ~ quant.start),
#            end = case_when(te_region.strand == "+" ~ quant.end,
#                              te_region.strand == "-" ~ cage.end),
#            te_region_id = te_region.names,
#            score = te_island.length,
#            strand = te_region.strand
#            ) %>% 
#     dplyr::select(chr, start, end, te_region_id, score, strand)
# 
# 
# ggplot(te.islands.bed, aes(score)) +
#     geom_density() +
#     geom_vline(xintercept = mean(te.islands.bed$score), color = "red") +
#     geom_text(x = log10(mean(te.islands.bed$score)), 
#               y = 0.0, 
#               label = paste0("mean = " ,round(mean(te.islands.bed$score),2)), 
#               angle = 90, 
#               vjust = -1, 
#               hjust = 0) +
#     geom_text(x = log10(max(te.islands.bed$score)) - 0.2, 
#               y = 0.75, 
#               label = paste0("n = " ,nrow(te.islands.bed)), 
#               
#               ) +
#     scale_x_continuous(trans = 'log10')
# 
