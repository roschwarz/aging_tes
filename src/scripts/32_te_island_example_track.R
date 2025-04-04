library(Gviz)
library(GenomicRanges)
library(blackRcloud)
library(plyranges)
library(BRGenomics)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
source('../lib/styles.R')
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


start_coord = 125248029
end_coord = 125255257
gtrack <- GenomeAxisTrack()





##########################
# Gene and TE-Annotation #
##########################
# Genes are not annotated within this region, hence the track is ignored

te.region.file <- paste0("../../data/shared/mm10_TE_region_5prime_extended.bed")

teRegionRanges <- genomation::readGeneric(te.region.file, 
                                          strand = 6, 
                                          meta.cols = list(names = 4))

teRanges <- blackRcloud::getRanges("mmusculus", 102, "te", 1)

#################
# Blood example #
#################

chrom = 'chr1'
strand = '-'
start_coord = 21322152
end_coord = 21326292
gtrack <- GenomeAxisTrack()


teRegio <- AnnotationTrack(teRegionRanges[teRegionRanges@elementMetadata$names %in% c("TE_Cluster_10350")],
                           name = 'te region',
                           id = c("TE island 10350"),
                           featureAnnotation = "id",
                           fontcolor.item = "white",
                           fill = "black",
                           size = 2)

TEs_of_interest = data.frame(teRanges) %>% 
    filter(seqnames == chrom, 
           strand == strand,
           dplyr::between(start, start_coord, end_coord)) %>% 
    pull(te.id) %>% 
    unique()

ids = data.frame(teRanges) %>% 
    filter(seqnames == chrom, 
           strand == strand,
           dplyr::between(start, start_coord, end_coord)) %>% 
    blackRcloud::splitTEID("te.id") %>% 
    pull(super_family)

#ids = c("L1", "L1", "L1","L1", "B1", "L1")
#ids = TEs_of_interest
fam_col = c(rep("#14273F", 7), 
            rep("#831335", 1), "#14273F", "#831335",
            rep("#8493AE",2), "#831335",  rep("#14273F", 4), "#8493AE", "#14273F", "#831335")

teInstance <- AnnotationTrack(teRanges[teRanges@elementMetadata$te.id %in% TEs_of_interest & !duplicated(teRanges@elementMetadata$te.id)],
                              featureAnnotation = "id",
                              #fontcolor.item = "white",
                              id = ids,
                              stacking = "dense",
                              fill = fam_col,
                              col = 'white',
                              size = 2)

#################
# Cage big wig  #
#################

cage_limits = c(0,80)
ticks =  c(0,20,40,60, 80)

cage_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/blood.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))

cage_old_track <- DataTrack(range = cage_old.bw,
                            type = "h",
                            name = 'old',
                            col = "darkgreen",
                            ylim = cage_limits,
                            yTicksAt = ticks)    

cage_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/blood.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))

cage_young_track <- DataTrack(range = cage_young.bw,
                              type = "h",
                              name = 'young',
                              col = "darkgreen", 
                              yTicksAt = ticks,
                              ylim = cage_limits)    

cage_peaks <- genomation::readGeneric('../../results/cage_RS/blood_gclipped/raw_peaks/segemehl/peakachu/blood.reverse.bed', 
                                      strand = 6, 
                                      meta.cols = list(names = 4))

cagePeaks <- AnnotationTrack(filter(cage_peaks, strand == strand, start > start_coord, end < end_coord, seqnames == chrom),
                             name = 'CAGE Peaks',
                             #id = c("TE island 1300406", "TE island 1300408"),
                             #featureAnnotation = "id",
                             fontcolor.item = "white",
                             col = "darkgreen",
                             fill = "darkgreen",
                             size = 2)


#################
# RNA-Seq track #
#################

rna_limits = c(0,300)
rna_ticks =  c(0,100,200,300)

rna_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/blood.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


rna_old_track <- DataTrack(range = rna_old.bw, 
                           type = "hist", 
                           name = 'old',
                           col = "black",
                           window = -1,
                           fill.histogram = "black",
                           col.boxplotFrame = "black",
                           ylim = rna_limits,
                           yTicksAt = rna_ticks
)    

rna_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/blood.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


rna_young_track <- DataTrack(
    range = rna_young.bw,
    type = "hist",
    name = 'young',
    col = "black",
    fill.histogram = "black",
    col.boxplotFrame = "black",
    window = -1,
    ylim = rna_limits,
    yTicksAt = rna_ticks
) 


#################
# quant-Seq track #
#################

quant_limits = c(0,10)
quant_ticks =  c(0,5,10)

quant_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/blood.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


quant_young_track <- DataTrack(
    range = quant_young.bw,
    type = "h",
    name = 'young',
    col = "red",
    ylim = quant_limits,
    yTicksAt = quant_ticks
) 


quant_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/blood.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


quant_old_track <- DataTrack(range = quant_old.bw, 
                             type = "h", 
                             name = 'old',
                             col = "red",
                             ylim = quant_limits,
                             yTicksAt = quant_ticks
)  

# Need the peaks of the negative strand because of the sequencing protocol
quant_peaks <- genomation::readGeneric('../../results/quant_RS/blood/raw_peaks/segemehl/peakachu/blood.forward.bed', 
                                       strand = 6, 
                                       meta.cols = list(names = 4))

quantPeaks <- AnnotationTrack(filter(quant_peaks, strand == "+", start > start_coord, end < end_coord, seqnames == chrom),
                              name = 'Quant Peaks',
                              fontcolor.item = "white",
                              fill = "red",
                              col = "red",
                              size = 2)

####################
# Highlight tracks #
####################


############
# ideogram #
############

ideoTrack <- IdeogramTrack(genome = "mm10", chromosome = "chr6")

###############################
# Panel Construction  Panel_7 #
###############################
cairo_pdf(paste0('../../results/figures/', '32_te_island_example_blood_2.5_7.8_300.pdf'), 
          height = 2.5,  width = 7.8)

Gviz::plotTracks(list(#ideoTrack,
    gtrack,
    teRegio,
    teInstance,
    cage_young_track,
    cage_old_track,
    cagePeaks,
    rna_young_track,
    rna_old_track,
    quantPeaks,
    quant_old_track,
    quant_young_track
),
add35 = TRUE,
labelPos = 'below',
scale = 0.25,
shape = 'box',
fontface.main = 'arial',
chromosome = chrom,
from = start_coord,
to = end_coord,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")


dev.off()




################
# Brain Example #
################

start_coord = 80997539 - 500
end_coord = 81005311 + 1000
chrom = 'chr13'
strand = '+'
gtrack <- GenomeAxisTrack()


teRegio <- AnnotationTrack(teRegionRanges[teRegionRanges@elementMetadata$names %in% c("TE_Cluster_419117")],
                           name = 'te region',
                           id = c("TE island 419117"),
                           featureAnnotation = "id",
                           fontcolor.item = "white",
                           fill = "black",
                           size = 2)

TEs_of_interest = data.frame(teRanges) %>% 
    filter(seqnames == chrom, 
           strand == strand,
           dplyr::between(start, start_coord, end_coord)) %>% 
    pull(te.id) %>% 
    unique()

ids = data.frame(teRanges) %>% 
    filter(seqnames == chrom, 
           strand == strand,
           dplyr::between(start, start_coord, end_coord)) %>% 
    blackRcloud::splitTEID("te.id") %>% 
    pull(super_family)

#ids = c("L1", "L1", "L1","L1", "B1", "L1")
#ids = TEs_of_interest
fam_col = c("#8493AE", rep("#831335", 5), "#8493AE", rep("#14273F", 2))

teInstance <- AnnotationTrack(teRanges[teRanges@elementMetadata$te.id %in% TEs_of_interest & !duplicated(teRanges@elementMetadata$te.id)],
                              featureAnnotation = "id",
                              #fontcolor.item = "white",
                              id = ids,
                              stacking = "dense",
                              fill = fam_col,
                              col = 'white',
                              size = 2)

#################
# Cage big wig  #
#################

cage_limits = c(0,30)
ticks =  c(0,10,20,30)

cage_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/brain.old.forward.down.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))

cage_old_track <- DataTrack(range = cage_old.bw,
                            type = "h",
                            name = 'old',
                            col = "darkgreen",
                            ylim = cage_limits,
                            yTicksAt = ticks)    

cage_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/brain.young.forward.down.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))

cage_young_track <- DataTrack(range = cage_young.bw,
                              type = "h",
                              name = 'young',
                              col = "darkgreen", 
                              yTicksAt = ticks,
                              ylim = cage_limits)    

# cage_peaks <- genomation::readGeneric('../../results/cage_RS/brain_gclipped/raw_peaks/segemehl/peakachu/brain.forward.bed', 
#                                       strand = 6, 
#                                       meta.cols = list(names = 4))

cage_peaks <- genomation::readGeneric('../../results/cage_RS/brain_downsampled_gclipped/raw_peaks/segemehl/peakachu/brain_downsampled.forward.bed', 
                                      strand = 6, 
                                      meta.cols = list(names = 4))

cagePeaks <- AnnotationTrack(filter(cage_peaks, strand == "+", start > start_coord, end < end_coord, seqnames == chrom),
                             name = 'CAGE Peaks',
                             #id = c("TE island 1300406", "TE island 1300408"),
                             #featureAnnotation = "id",
                             fontcolor.item = "white",
                             col = "darkgreen",
                             fill = "darkgreen",
                             size = 2)


#################
# RNA-Seq track #
#################

rna_limits = c(0,60)
rna_ticks =  c(0,20,40,60)

rna_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/brain.old.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


rna_old_track <- DataTrack(range = rna_old.bw, 
                           type = "hist", 
                           name = 'old',
                           col = "black",
                           window = -1,
                           fill.histogram = "black",
                           col.boxplotFrame = "black",
                           ylim = rna_limits,
                           yTicksAt = rna_ticks
)    

rna_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/brain.young.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


rna_young_track <- DataTrack(
    range = rna_young.bw,
    type = "hist",
    name = 'young',
    col = "black",
    fill.histogram = "black",
    col.boxplotFrame = "black",
    window = -1,
    ylim = rna_limits,
    yTicksAt = rna_ticks
) 


#################
# quant-Seq track #
#################

quant_limits = c(0,80)
quant_ticks =  c(0,20,40,60,80)

quant_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/brain.young.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


quant_young_track <- DataTrack(
    range = quant_young.bw,
    type = "h",
    name = 'young',
    col = "red",
    ylim = quant_limits,
    yTicksAt = quant_ticks
) 


quant_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/brain.old.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chrom, ':', start_coord, '-', end_coord), "GRanges"))


quant_old_track <- DataTrack(range = quant_old.bw, 
                             type = "h", 
                             name = 'old',
                             col = "red",
                             ylim = quant_limits,
                             yTicksAt = quant_ticks
)  

# Need the peaks of the negative strand because of the sequencing protocol
quant_peaks <- genomation::readGeneric('../../results/quant_RS/brain/raw_peaks/segemehl/peakachu/brain.reverse.bed', 
                                       strand = 6, 
                                       meta.cols = list(names = 4))

quantPeaks <- AnnotationTrack(filter(quant_peaks, strand == "-", start > start_coord, end < end_coord, seqnames == chrom),
                              name = 'Quant Peaks',
                              fontcolor.item = "white",
                              fill = "red",
                              col = "red",
                              size = 2)

####################
# Highlight tracks #
####################


############
# ideogram #
############

ideoTrack <- IdeogramTrack(genome = "mm10", chromosome = "chr6")

###############################
# Panel Construction  Panel_7 #
###############################
cairo_pdf(paste0('../../results/figures/', '32_te_island_example_brain_5.3_2.4_300.pdf'), 
          height = 2.4,  width = 5.3)

Gviz::plotTracks(list(#ideoTrack,
    gtrack,
    teRegio,
    teInstance,
    cage_young_track,
    cage_old_track,
    cagePeaks,
    rna_young_track,
    rna_old_track,
    quantPeaks,
    quant_old_track,
    quant_young_track
),
add35 = TRUE,
labelPos = 'below',
scale = 0.2,
shape = 'box',
chromosome = chrom,
from = start_coord,
to = end_coord,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")


dev.off()

################
# Skin Example #
################

teRegio <- AnnotationTrack(teRegionRanges[teRegionRanges@elementMetadata$names %in% c("TE_Cluster_1300408", "TE_Cluster_1300406")],
                           name = 'te region',
                           id = c("TE island 1300406", "TE island 1300408"),
                           featureAnnotation = "id",
                           fontcolor.item = "white",
                           fill = "black",
                           size = 2)

TEs_of_interest = data.frame(teRanges) %>% 
    filter(seqnames == "chr6", 
           strand == "-",
           dplyr::between(start, start_coord, end_coord)) %>% 
    pull(te.id) %>% 
    unique()

ids = data.frame(teRanges) %>% 
    filter(seqnames == "chr6", 
           strand == "-",
           dplyr::between(start, start_coord, end_coord)) %>% 
    blackRcloud::splitTEID("te.id") %>% 
    pull(super_family)

#ids = c("L1", "L1", "L1","L1", "B1", "L1")
#ids = TEs_of_interest
fam_col = c("#8493AE", "#831335", rep("#14273F",7), "#831335", rep("#14273F",8), "#8493AE")


teInstance <- AnnotationTrack(teRanges[teRanges@elementMetadata$te.id %in% TEs_of_interest & !duplicated(teRanges@elementMetadata$te.id)],
                              #featureAnnotation = "id",
                              #fontcolor.item = "white",
                              #id = ids,
                              stacking = "dense",
                              fill = fam_col,
                              col = 'white',
                              size = 2)

#################
# Cage big wig  #
#################

cage_limits = c(0,20)
ticks =  c(0,10,20)

cage_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinII.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr6:', start_coord, '-', end_coord), "GRanges"))

cage_old_track <- DataTrack(range = cage_old.bw,
                            type = "h",
                            name = 'old',
                            col = "darkgreen",
                            ylim = cage_limits,
                            yTicksAt = ticks)    

cage_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinII.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr6:', start_coord, '-', end_coord), "GRanges"))

cage_young_track <- DataTrack(range = cage_young.bw,
                              type = "h",
                              name = 'young',
                              col = "darkgreen", 
                              yTicksAt = ticks,
                              ylim = cage_limits)    

#################
# RNA-Seq track #
#################

rna_limits = c(0,30)
rna_ticks =  c(0,10,20,30)

rna_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/skinII.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr6:', start_coord, '-', end_coord), "GRanges"))


rna_old_track <- DataTrack(range = rna_old.bw, 
                           type = "hist", 
                           name = 'old',
                           col = "black",
                           window = -1,
                           fill.histogram = "black",
                           col.boxplotFrame = "black",
                           ylim = rna_limits,
                           yTicksAt = rna_ticks
)    

rna_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/skinII.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr6:', start_coord, '-', end_coord), "GRanges"))


rna_young_track <- DataTrack(
    range = rna_young.bw,
    type = "hist",
    name = 'young',
    col = "black",
    fill.histogram = "black",
    col.boxplotFrame = "black",
    window = -1,
    ylim = rna_limits,
    yTicksAt = rna_ticks
) 

cage_peaks <- Gviz::AnnotationTrack(GRanges(seqnames = c('chr6', 'chr6'),
                                          ranges = IRanges(start = c(125253975, 125252291),
                                                           end = c(125254287,  125252406),
                                                           names = rep('chr6', 2)),
                                          strand = rep("-",2),
                                          name = rep("CAGE", 2)
                                          ),
                                  name = "CAGE",
                                  size = 2,
                                  col = "darkgreen",
                                  fill = "darkgreen")

#################
# quant-Seq track #
#################

quant_limits = c(0,60)
quant_ticks =  c(0,20,40,60)

quant_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/skinII.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr6:', start_coord, '-', end_coord), "GRanges"))


quant_young_track <- DataTrack(
    range = quant_young.bw,
    type = "h",
    name = 'young',
    col = "red",
    ylim = quant_limits,
    yTicksAt = quant_ticks
) 


quant_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/skinII.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr6:', start_coord, '-', end_coord), "GRanges"))


quant_old_track <- DataTrack(range = quant_old.bw, 
                             type = "h", 
                             name = 'old',
                             col = "red",
                             ylim = quant_limits,
                             yTicksAt = quant_ticks
)  


quant_peaks <- Gviz::AnnotationTrack(GRanges(seqnames = c('chr6', 'chr6', 'chr6'),
                                          ranges = IRanges(start = c(125253171, 125251046, 125249046),
                                                           end = c(125253281,  125251153, 125249187),
                                                           names = rep('chr6', 3)),
                                          strand = rep("-",3),
                                          name = rep("Quant", 3)
                                          ),
                                  name = "Quant",
                                  size = 2,
                                  col = "red",
                                  fill = "red")

####################
# Highlight tracks #
####################


############
# ideogram #
############

ideoTrack <- IdeogramTrack(genome = "mm10", chromosome = "chr6")

###############################
# Panel Construction  Panel_7 #
###############################
cairo_pdf(paste0('../../results/figures/', '32_te_island_example_2.4_7.8_300.pdf'), 
    height = 2.4,  width = 7.8)

Gviz::plotTracks(list(#ideoTrack,
    gtrack,
    teRegio,
    teInstance,
    cage_young_track,
    #cage_old_track,
    cage_peaks,
    rna_young_track,
    #rna_old_track,
    quant_peaks,
    #quant_old_track,
    quant_young_track
),
add35 = TRUE,
labelPos = 'below',
scale = 0.3,
shape = 'box',
chromosome = "chr6",
from = start_coord,
to = end_coord,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")
dev.off()

ggsave(track,
       filename = paste0('../../results/figures/', '32_te_island_example_5.3_2.4_300.pdf'),
       device = cairo_pdf,
       width = 5.3,
       height = 2.4,
       units = "in",
       dpi = 300
)




# ggsave(track,
#        filename = paste0('../../results/figures/', '32_te_island_example_5.3_2.4_300.pdf'),
#        device = cairo_pdf,
#        width = 5.3,
#        height = 2.4,
#        units = "in",
#        dpi = 300
# )