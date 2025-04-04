# --------------------------------- Notes--- -----------------------------------
# This script is used to identify and caracterize TE Islands based on 
# CAGE-, RNA-, and Quantdata
# 
# A TE-Island is defined as a TE region ("artificial annotation of TE Island"; see buildTE-region.annotation.sh)
# that intersects with at least one CAGE- and Quant-peak and does have at leas on expressed TE inside (!is.na(padj))
#
#
# Output data (../data/rna; ../tables):
#
# Load Environment
if (!exists("ENVIRONMENT_LOADED")) {
    
    source("./01_load_environment.R")
    
} else if (!ENVIRONMENT_LOADED) {
    
    source("./01_load_environment.R")
    
}

rna.deseq.res.file <- 'data/processed/rna_RS_deseq_TE_instances_SalmonTE.Rdata'
#rna.deseq.res.file <- '/misc/paras/data/rschwarz/workspace/phd/projects/mCQuaRna/packages/publication_ready/data/rna/deseq.TE.SalmonTE.Rdata'

rnaDeseq <- loadRdata(rna.deseq.res.file)

teRegionInstances <-
    create_TE_Instances(teRegionRanges = teRegionRanges, 
                        teRanges = teRanges)


teRegionLength <- data.frame(teRegionRanges) %>%
    dplyr::select(names, width) %>%
    mutate(tissue = 'background')

teRegionMemberCount <- teRegionInstances %>%
    dplyr::count(te_region_id)

teRegionLength <- merge(teRegionLength,
                        teRegionMemberCount,
                        by.x = 'names',
                        by.y = 'te_region_id')

order_legend <- legendBuilder(teRegionInstances)

# Collect TE regions with at least one expressed TE. All TEs with an adjusted p-value in the DESeq2 analysis
# of the RNA-Seq data are considered as expressed
expr_te_regions <- sapply(names(rnaDeseq), simplify = F, function(x){
    
    expr_te_rna <- rownames(rnaDeseq[[x]]  %>% filter(!is.na(padj)))
    
    expr_regions <- teRegionInstances %>% 
        filter(te_id %in% expr_te_rna) %>% 
        dplyr::pull('te_region_id')
    
    return(unique(expr_regions)) # unique, as multiple TEs per island can be expressed, 
    # so that TE Islands can be counted multiple times
})


# Get the set of expressed TE regions that also intersect with a CAGE-peak (individually expressed TE regions)

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

idv_teRegions <- sapply(names(cageRanges), 
                       simplify = F, 
                       USE.NAMES = T, 
                       function(x){
                           
                           intersectGranger(cageRanges[[x]], TEregionRanges_5prime, tab = 'all')
                           
                       })


idv_expr_teRegions <- sapply(tissues, 
                             simplify = F, 
                             USE.NAMES = T, 
                             function(tissue){
                                 
                                 idv_teRegions[[tissue]] %>% 
                                     filter(subject.names %in% expr_te_regions[[tissue]]) %>% 
                                 filter(!duplicated(subject.names))
                                 
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



sto_teRegions <- sapply(names(quantRanges), 
                          simplify = F, 
                          USE.NAMES = T, 
                          function(x){
                              
                              intersectGranger(quantRanges[[x]], TEregionRanges_3prime, tab = 'all')
                              
                          })

sto_expr_teRegions <- sapply(tissues, 
                             simplify = F, 
                             USE.NAMES = T, 
                             function(tissue){
                                 
                                 df <- sto_teRegions[[tissue]] %>% 
                                     filter(subject.names %in% expr_te_regions[[tissue]]) %>% 
                                 filter(!duplicated(subject.names))
                                 
                                 
                                 bed_file <- df %>% dplyr::select(subject.seqnames,
                                                                  subject.start,
                                                                  subject.end,
                                                                  subject.names,
                                                                  subject.width,
                                                                  subject.strand)
                                 
                                 filename <- paste0('data/shared/', tissue,
                                                    '_TE_island_w_quant.bed')
                                 
                                 write.table(bed_file, file = filename, sep = '\t', col.names = F, row.names = F, quote = F)
                                 return(df)
                                 
                             })



# The definition of TE island is the intersection of a TE region with a CAGE-, RNA-, and Quant-signal.
# Therefore, I just need to filter the individually TE regions for TE regions ids that are contained
# in the set of TE regions that contain a quant signal and are expressed (sto_expr).
teIslands <- sapply(tissues, 
                    simplify = F, 
                    USE.NAMES = T, 
                    function(tissue){
                                                
                        idv_teRegions[[tissue]] %>% 
                            filter(subject.names %in% sto_expr_teRegions[[tissue]]$subject.names) %>% 
                            filter(!duplicated(subject.names))
                        
                    })



# Next, I want to annotate canonical TE Islands, which means I want to set the start coordinate to the 
# start coordinate of the left most CAGE peak and end the stop coordinate to the end coordinate of the 
# right most Quant peak.
# 
# Filter for TE Islands with a width > 0, as when the Quant-peak is located in front of the CAGE peak
# you will get a negative length and this is not a real TE island.
cannonical_teIslands <- do.call('rbind', sapply(tissues, 
                                                simplify = F, 
                                                USE.NAMES = T, 
                                                function(tissue){
                                                    
                                                    # collect the first cage peak in a TE island
                                                    cage_tmp <- idv_teRegions[[tissue]] %>% 
                                                        group_by(subject.names) %>% 
                                                        filter(case_when(subject.strand == "-" ~ query.end == (max(query.end)),
                                                                         subject.strand == "+" ~ query.start == (min(query.start))
                                                        )
                                                        )
                                                    
                                                    names(cage_tmp) <- str_replace_all(names(cage_tmp), 
                                                                                       c('query' = 'cage', 'subject' = 'te_region'))
                                                    
                                                    # collect the last Quant peak in a TE island
                                                    quant_tmp <- sto_teRegions[[tissue]] %>% 
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
                                                    
                                                    df <- df %>% 
                                                        filter(te_region.names %in% expr_te_regions[[tissue]]) %>% 
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

#### Plotting the data ##########
rpl_density(cannonical_teIslands %>% 
                 dplyr::rename(width = te_island.width)) +
    theme(legend.position = c(0.1,0.8))


n_func <- function(df, y = 1){
    return(data.frame(y = y, 
                      label = paste0("n = ", length(df))))
}

max_func <- function(df, y = 4.3){
    
    return(data.frame(y,
                      label = paste0("max = ", 10^max(df), " bp")))
    
}

min_func <- function(df, y = 1.35){
    
    return(data.frame(y,
                      label = paste0("min = ", 10^min(df), " bp")))
    
}

mean_func <- function(df){
    
    return(data.frame(y = mean(df),
                      label = paste0(round(10^mean(df), 2), " bp")))
    
}

cannonical_teIslands$tissue <- factor(cannonical_teIslands$tissue, levels =  c('brain', 'skin', 'blood') )

cannonical_teIslands %>% 
    filter(te_island.width > 50) %>% 
    group_by(tissue) %>% 
    summarise(max.length = max(te_island.width),
              min.length = min(te_island.width))

# Filtered for TE islands with a length bigger than 50 bp
pl_length_violin <- ggplot(cannonical_teIslands %>% filter(te_island.width > 50), aes(tissue, te_island.width)) +
    geom_violin(trim = FALSE,
                aes(fill = tissue),
                color = NA,
                #scale = "count"
    ) +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6) +
    scale_y_log10(labels = scales::comma,
                  expand = expansion(mult = c(0.05, 0.0))) +
    stat_summary(geom = "text",
                 fun.data = n_func,
                 fun.args = c(y = 1.3),
                 size = 6/.pt) +
    # stat_summary(geom = "text",
    #              fun.data = max_func,
    #              fun.args = c(y = 4.5),
    #              size = 6/.pt) +
    stat_summary(geom = "text",
                 fun.data = mean_func,
                 #fun.args = c(y = 4.5),
                 size = 6/.pt,
                 vjust = 1.5,
                 color = "black",
                 angle = 90) +
    # stat_summary(geom = "text",
    #              fun.data = min_func,
    #              fun.args = c(y = 1.4),
    #              size = 6/.pt) +
    scale_fill_manual(values = tissue.color) +
    labs(y = "Length of TE island in bp") +
    theme(legend.position = "None",
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.title.x = element_blank(),
          panel.grid = element_blank())

pl_length_violin

ggsave(pl_length_violin,
       filename = paste0(figure_dir, '30_te_island_length_dis_2x2_300.pdf'),
       device = cairo_pdf,
       width = 2,
       height = 2,
       units = "in",
       dpi = 300
)






################
# Double Checking if I really update the TE Island instance stuff

cannonical_teIslands_composition <- do.call('rbind',
                                            sapply(tissues, simplify = FALSE, function(x) {
                                                teIslandBed <- cannonical_teIslands %>%
                                                    filter(tissue == x) %>%
                                                    dplyr::select(
                                                        te_region.seqnames,
                                                        te_island.start,
                                                        te_island.end,
                                                        te_region.names,
                                                        te_island.width,
                                                        te_region.strand
                                                    ) %>%
                                                    dplyr::rename(
                                                        chr = te_region.seqnames,
                                                        start = te_island.start,
                                                        end = te_island.end,
                                                        te_island_id = te_region.names,
                                                        score = te_island.width,
                                                        strand = te_region.strand
                                                    )
                                                
                                                teIslandRanges <-
                                                    GRanges(
                                                        seqnames = teIslandBed$chr,
                                                        ranges = IRanges(
                                                            teIslandBed$start,
                                                            end = teIslandBed$end,
                                                            names = teIslandBed$te_island_id
                                                        ),
                                                        strand = teIslandBed$strand,
                                                        te_island_id = teIslandBed$te_island_id
                                                    )
                                                
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
                                                        "query.te_island_id",
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
                                                        "te_region_start",
                                                        "te_region_end",
                                                        "te_region_strand",
                                                        "te_region_width",
                                                        "te_region_id",
                                                        "te_start",
                                                        "te_end",
                                                        "te_strand",
                                                        "te_id",
                                                        "te_position"
                                                    )
                                                
                                                teIslandInstances <- 
                                                    blackRcloud::splitTEID(teIslandInstances, 'te_id') %>% 
                                                    mutate(super_family = case_when(super_family == 'Alu' ~ 'B1',
                                                                                    .default = super_family))
                                                
                                                
                                                teIslandInstances$tissue <- x
                                                
                                                return(teIslandInstances)
                                            })
)



save(cannonical_teIslands_composition, file = paste0(data_dir, "30_TE_island_composition.Rdata"))

### TE density within TEItx

mean_dens_func <- function(df){
    
    return(data.frame(y = mean(df),
                      label = paste0(round(mean(df), 2))))
    
}

max_dens_func <- function(df){
    
    return(data.frame(y = max(df) + 1,
                      label = paste0("n(max) = ", max(df))))
    
}

teitx_density <- cannonical_teIslands_composition %>% 
    filter(te_region_width > 50) %>% 
    group_by(tissue, te_region_id) %>% 
    dplyr::count() %>% 
    order_tissue()


teitx_density %>% group_by(tissue) %>% summarise(n(), mean(n), median(n)) 

teitx_density_pl <- ggplot(teitx_density, aes(tissue, n)) +
    geom_violin(aes(fill = tissue), adjust = 2,
                color = NA,) +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6) +
     # stat_summary(geom = "text",
     #              fun.data = max_dens_func,
     #              size = 6/.pt) +
    stat_summary(geom = "text",
                 fun.data = mean_dens_func,
                 size = 6/.pt,
                 vjust = 1.5,
                 color = "black",
                 angle = 90) +
    scale_fill_manual(values = tissue.color) +
    labs(y = "# TEs in TE island") +
    theme(legend.position = "None",
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.title.x = element_blank(),
          panel.grid = element_blank())


### Diversity of TEItx
teitx_diversity <- cannonical_teIslands_composition %>% 
    filter(te_region_width > 50) %>% 
    group_by(tissue, te_region_id) %>% 
    filter(!duplicated(family)) %>% 
    dplyr::count() %>% 
    order_tissue()

teitx_diversity %>% group_by(tissue) %>% summarise(n(), mean(n), median(n)) 

teitx_diversity_pl <- ggplot(teitx_diversity, aes(tissue, n)) +
    geom_violin(aes(fill = tissue), adjust = 2, color = NA,) +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6) +
    # stat_summary(geom = "text",
    #              fun.data = max_dens_func,
    #              size = 6/.pt) +
    stat_summary(geom = "text",
                 fun.data = mean_dens_func,
                 size = 6/.pt,
                 vjust = 1.5,
                 color = "black",
                 angle = 90) +
    scale_fill_manual(values = tissue.color) +
    labs(y = "# TE families in TE island") +
    theme(legend.position = "None",
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.title.x = element_blank(),
          panel.grid = element_blank())



teitx_char_violines <- ggarrange(pl_length_violin, teitx_density_pl, teitx_diversity_pl, nrow = 1)

ggsave(teitx_char_violines,
       filename = paste0(figure_dir, '30_te_island_char_dis_11x4_300.pdf'),
       device = cairo_pdf,
       width = 11,
       height = 4,
       units = "cm",
       dpi = 300
)

#########
# RESIS #
#########

write.csv(cannonical_teIslands %>% filter(te_island.width > 50),
          file = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_length.csv')

ggsave(pl_length_violin,
       filename = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_length.pdf',
       device = cairo_pdf,
       width = 5,
       height = 5,
       units = "cm",
       dpi = 300 )

write.csv(teitx_density, file = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_density.csv')

ggsave(teitx_density_pl,
       filename = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_density.pdf',
       device = cairo_pdf,
       width = 5,
       height = 5,
       units = "cm",
       dpi = 300 )

write.csv(teitx_diversity, file = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_diversity.csv')

ggsave(teitx_diversity_pl,
       filename = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_diversity.pdf',
       device = cairo_pdf,
       width = 5,
       height = 5,
       units = "cm",
       dpi = 300 )

### Composition of TEItx

cannonical_teIslands_composition_super <- do.call('rbind', sapply(tissues, simplify = FALSE, function(x){
    df <- cannonical_teIslands_composition %>% 
        filter(te_region_width > 50) %>%  
        filter( tissue == x)
    
    order_counts <- prepare_super_composition(df)
    order_counts$tissue <- x
    
    return(order_counts)
})
)


cannonical_teIslands_composition_super$tissue <- factor(cannonical_teIslands_composition_super$tissue,
                                                        levels = c('brain', 'skin', 'blood'))


teIsland_super_bar <- cannonical_teIslands_composition_super %>% 
    filter(per >= 0.5) %>%
    ggplot(aes(x = super_family, y = per, fill = order)) +
    geom_col() + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_flip() +
    scale_fill_manual(values = order.color) +
    facet_grid(.~tissue) +
    labs(y = "Proportion of TE island [%]") +
    theme_rob(6) +
    theme(
        legend.position = 'bottom',
        legend.background = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        #strip.text = element_blank(),
        panel.spacing.x = unit(0.2,'cm')
        #plot.margin = unit(c(0, 0, 0, 0), "mm"),
        
    )

teIsland_super_bar

ggsave(teIsland_super_bar,
       filename = paste0(figure_dir, '30_super_bar_4x3_300.pdf'),
       device = cairo_pdf,
       width = 4.2,
       height = 3,
       units = "in",
       dpi = 300
)

#########
# RESIS #
#########

write.csv(cannonical_teIslands_composition_super %>% filter(per >= 0.5), 
          file = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_composition.csv')



ggsave(teIsland_super_bar,
       filename = './manuscripts/nature_aging/resis/figures/figure7/7D_TE_island_composition.pdf',
       device = cairo_pdf,
       width = 4.2,
       height = 3,
       units = "in",
       dpi = 300
)

############# Merge Figures to a sub-panel #######

legendCol <- get_legend(teIsland_super_bar +
                            theme(legend.margin = margin(0,0,0,0),
                                  legend.title = element_text(size = 6), 
                                  legend.text = element_text(size = 6)
                            )  )

teitx_char_violines <- ggarrange(pl_length_violin, 
                                 teitx_density_pl, 
                                 teitx_diversity_pl,
                                 nrow = 1)

teitx_char <- 
    ggarrange(teitx_char_violines, 
              teIsland_super_bar + theme(legend.position = 'None'), 
              legendCol,
              nrow = 3,
              heights = c(0.3, 0.6, 0.1),
              align = 'v'
    )


teitx_char

# grid.arrange(teitx_char_violines,
#              teIsland_super_bar,
#              nrow = 2,
#              heights )

ggsave(teitx_char,
       filename = paste0(figure_dir, '30_te_island_char_10x12_300.pdf'),
       device = cairo_pdf,
       width = 10,
       height = 12,
       units = "cm",
       dpi = 300
)






c_teIslands_comp_order <- do.call('rbind',
                                  sapply(tissues, simplify = FALSE, function(x) {
                                      
                                      df <- cannonical_teIslands_composition %>% filter(tissue == x)
                                      
                                      order_counts <-
                                          prepare_Composition(df)
                                      order_counts$tissue <-
                                          x
                                      
                                      return(order_counts)
                                  }))


ggplot(c_teIslands_comp_order %>% filter(tissue == 'skin'), aes(x = "", y = per, fill = order)) +
    geom_col(color = "black") +
    coord_polar("y", start = 0) +
    geom_label(aes(label = labels), position = position_stack(vjust = 0.3), color = "white", size = 10/.pt) +
    guides(fill = guide_legend(title = "Order of TEs")) +
    scale_fill_manual(values = order.color) +
    facet_grid(.~tissue) +
    theme_void() +
    theme(legend.position = "bottom",
          plot.margin = unit(c(0, -10, 0, 12), "mm"), #top, right, bottom, left
          strip.text = element_text(size = 12),
          #panel.spacing.x = unit(1, 'cm')
    ) 






c_teIslands_comp_order$tissue <- factor(c_teIslands_comp_order$tissue,
                                                        levels = c('brain', 'skin', 'blood'))


c_teIslands_comp_order$order <- factor(c_teIslands_comp_order$order,
                                       levels = c('LINE', 'SINE', 'DNA',  'LTR' ))

teIsland_pies <- ggplot(c_teIslands_comp_order, aes(x = "", y = per, fill = order)) +
    geom_col(color = "black", position = position_fill()) +
    coord_polar("y", start = 0) +
    geom_label(aes(label = labels), position = position_fill(vjust = 0.3), color = "white", size = 6/.pt) +
    guides(fill = guide_legend(title = "Order of TEs")) +
    scale_fill_manual(values = order.color) +
    facet_grid(.~tissue) +
    theme_void() +
    theme(legend.position = "None",
          plot.margin = unit(c(0, -5, 0, -5), "mm"), #top, right, bottom, left
          
          strip.text = element_text(size = 10),
          #panel.spacing.x = unit(1, 'cm')
    ) 

teIsland_pies

ggsave(teIsland_pies,
       filename = paste0(figure_dir, '30_pies_4x2_300.pdf'),
       device = cairo_pdf,
       width = 4,
       height = 1.5,
       units = "in",
       dpi = 300
)


#####




# Store bedfiles of the TE Islands for deeptools
# Store bedfiles of promoter of canonical TE Islands for Homer
sapply(tissues, function(x){
    
    bed_file <- cannonical_teIslands %>% 
        filter(tissue == x, te_island.width > 50) %>% 
        dplyr::select(cage.seqnames,
                      te_island.start,
                      te_island.end,
                      te_region.names,
                      te_island.width,
                      te_region.strand)
    filename <- paste0(common_data, x,
                       '_can_TE_island.bed')
    
    write.table(bed_file, file = filename, sep = '\t', col.names = F, row.names = F, quote = F)
    
    bed_file_promoter <- cannonical_teIslands %>% 
        filter(tissue == x, te_island.width > 50) %>% 
        dplyr::select(cage.seqnames,
                      te_island.start,
                      te_island.end,
                      te_region.names,
                      te_island.width,
                      te_region.strand) %>% 
        mutate(te_island_promoter.start = case_when(te_region.strand == '+' ~ te_island.start - 500,
                                           .default = te_island.end),
               te_island_promoter.end = case_when(te_region.strand == '-' ~ te_island.end + 500,
                                         .default = te_island.start)) %>% 
        dplyr::select(cage.seqnames,
                      te_island_promoter.start,
                      te_island_promoter.end,
                      te_region.names,
                      te_island.width,
                      te_region.strand)
    
    filename_promoter <- paste0(common_data, x,
                       '_promoter_can_TE_island.bed')
    
    
    write.table(bed_file_promoter, file = filename_promoter, sep = '\t', col.names = F, row.names = F, quote = F)
    
    
    return(cat(paste(filename, "hase been written.\n")))
    
})


# CAGE AND Quant-peaks are enriched in B1 elements. 
# Do we see something special in the composition of the B1 and L1 family?

b1_composition <- cannonical_teIslands_composition %>% 
    filter(super_family == 'B1') %>% 
    group_by(tissue) %>% 
    dplyr::count(family) %>% 
    mutate(per = n/sum(n)*100)


ggplot(b1_composition, aes(family, per, fill = tissue)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = tissue.color) +
    coord_flip()

L1_composition <- cannonical_teIslands_composition %>% 
    filter(super_family == 'L1') %>% 
    group_by(tissue) %>% 
    dplyr::count(family) %>% 
    mutate(per = n/sum(n)*100)


ggplot(L1_composition %>% filter(per > 0.5), aes(family, per, fill = tissue)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = tissue.color) +
    coord_flip()


# Enrichment analysis

brain.l1 <- cannonical_teIslands_composition %>% 
    filter(super_family == 'L1', tissue == 'brain') %>% 
    pull(family)

df <- binfamilydf <- binoRich(cannonical_teIslands_composition %>% filter(tissue == 'brain'),
               te.annotation, c(target = 'family',
                                    background = 'family'),
               FDR.cap = 1e-10)
df <- df %>% 
    filter(n.target >= 10) %>% 
    filter(!grepl("[?]", category))


cat.Sort <- df %>% 
    group_by(category) %>% 
    summarise(meanOdds = sum(ratio)) %>% 
    arrange(desc(meanOdds))


x.max <- max(df$ratio) + 0.2

df$category <-
    factor(df$category, levels = cat.Sort$category)

ggplot(df %>% filter(category %in% brain.l1), aes(ratio, category, size = log10.padj)) +
    geom_point() +
    #scale_color_manual(values = tissue.color) +
    scale_size(range = c(3, 12), name = "log10(FDR)") +
    labs(y = 'TE family', 
         x = 'log(odds ratio)') +
    xlim(c(-x.max, x.max)) +
    geom_vline(xintercept = 0, linetype = 1) +
    guides(size = guide_legend(nrow = 1)) +
    #facet_grid(category~., scale = 'free_y') +
    theme_rob(base_size = 12, base_family = 'arial') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          legend.position = 'None',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, size = 1.2),
          legend.box = "vertical") 


foo <- cannonical_teIslands_composition %>% dplyr::select(super_family, 
                                                          Kimura)

foo$cat <- 'te_island'

foo.bg <- te.annotation %>% dplyr::select(super_family, 
                                                          Kimura)

foo.bg$cat <- 'background'


foo <- rbind(foo, foo.bg)

foo$Kimura <- as.numeric(foo$Kimura)

ggplot(foo, aes(Kimura, cat, fill = cat)) +
    geom_violin()

ggplot(foo, aes(cat,Kimura, fill = cat)) +
    geom_violin(trim = FALSE,
                #scale = "count"
    ) +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 linewidth = 0.1,
                 outlier.size = 0.7,
                 outlier.color = "white",
                 outlier.fill = "black",
                 outlier.stroke = 1,
                 outlier.alpha = 0.6)
