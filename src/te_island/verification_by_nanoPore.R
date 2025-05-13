library(Gviz)
library(GenomicRanges)
library(blackRcloud)
library(plyranges)
library(BRGenomics)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_annotations()
aging_tes::load_te_island_env()

load_te_ranges()


# Search for TE islands that are covered by Nanopore
brain_te_islands <- indie_te_island_bed$brain

brain_counts <- merge(brain_counts, brain_intergenic_te_islands, by.x = 'Name', by.y = 'V4')


# Genome Browser shots

# TE island 305170

## TE_Cluster_305170
chromosome = 'chr12'
start_coord = 21558236
end_coord = 21563700


genomeBrowserNanopore(
    c("TE_Cluster_305170"),
    chromosome = 'chr12',
    start_coord = 21559250,
    end_coord = 21563700,
    te_island_annotation_file = te_island_file_5_prime_extended,
    bam_file = bam_file_dataset_1
)


genomeBrowserNanopore(
    c("TE_Cluster_399728"),
    chromosome = 'chr13',
    start_coord = 49472663,
    end_coord = 49472663,
    te_island_annotation_file = te_island_file_5_prime_extended,
    bam_file = bam_file_dataset_1
)


############# Old needs to be cleaned or serves as a template #########

# TE_Cluster_399728"
chromosome = 'chr13'
start_coord = 49472663
end_coord = 49478404

gtrack <- GenomeAxisTrack()



te_island_ranges <- genomation::readGeneric(te_region_file, 
                                          strand = 6, 
                                          meta.cols = list(names = 4))


teRegio <- AnnotationTrack(teRegionRanges[teRegionRanges@elementMetadata$names %in% c("TE_Cluster_399728")],
                           name = 'te region',
                           id = c("TE island 399728"),
                           featureAnnotation = "id",
                           fontcolor.item = "white",
                           fill = "black",
                           size = 2)

TEs_of_interest = data.frame(teRanges) %>% 
    filter(seqnames == chromosome, 
           strand == strand,
           dplyr::between(start, start_coord, end_coord)) %>% 
    pull(te.id) %>% 
    unique()

ids = data.frame(teRanges) %>% 
    filter(seqnames == chromosome, 
           strand == strand,
           dplyr::between(start, start_coord, end_coord)) %>% 
    blackRcloud::splitTEID("te.id") %>% 
    pull(super_family)

teInstance <- AnnotationTrack(teRanges[teRanges@elementMetadata$te.id %in% TEs_of_interest & !duplicated(teRanges@elementMetadata$te.id)],
                              featureAnnotation = "id",
                              #fontcolor.item = "white",
                              id = ids,
                              stacking = "dense",
                              col = 'white',
                              size = 2)


nanopore_reads <- Gviz::AnnotationTrack(
    range = 'results/te_island/long_read_verification/SRR24578270.bam',
    genome = "mm10",
    name = "reads",
    fill = "lightgray",
    chromosome = chromosome,
    id = c()
)

#################
# Cage big wig  #
#################

cage_limits = c(0,80)
ticks =  c(0,20,40,60, 80)

cage_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/brain.old.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chromosome, ':', start_coord, '-', end_coord), "GRanges"))

cage_old_track <- DataTrack(range = cage_old.bw,
                            type = "h",
                            name = 'old',
                            col = "darkgreen",
                            ylim = cage_limits,
                            yTicksAt = ticks)    

cage_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/brain.young.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chromosome, ':', start_coord, '-', end_coord), "GRanges"))

cage_young_track <- DataTrack(range = cage_young.bw,
                              type = "h",
                              name = 'young',
                              col = "darkgreen", 
                              yTicksAt = ticks,
                              ylim = cage_limits)    

cage_peaks <- genomation::readGeneric('results/cage_seq_gclipped/brain_downsampled_gclipped/raw_peaks/segemehl/peakachu/brain_downsampled.forward.bed', 
                                      strand = 6, 
                                      meta.cols = list(names = 4))

cagePeaks <- AnnotationTrack(filter(cage_peaks, strand == strand, start > start_coord, end < end_coord, seqnames == chromosome),
                             name = 'CAGE Peaks',
                             #id = c("TE island 1300406", "TE island 1300408"),
                             #featureAnnotation = "id",
                             fontcolor.item = "white",
                             col = "darkgreen",
                             fill = "darkgreen",
                             size = 2)


#################
# quant-Seq track #
#################

quant_limits = c(0,10)
quant_ticks =  c(0,5,10)

quant_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/brain.young.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chromosome, ':', start_coord, '-', end_coord), "GRanges"))


quant_young_track <- DataTrack(
    range = quant_young.bw,
    type = "h",
    name = 'young',
    col = "red",
    ylim = quant_limits,
    yTicksAt = quant_ticks
) 


quant_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/quant/mm10/brain.old.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chromosome, ':', start_coord, '-', end_coord), "GRanges"))


quant_old_track <- DataTrack(range = quant_old.bw, 
                             type = "h", 
                             name = 'old',
                             col = "red",
                             ylim = quant_limits,
                             yTicksAt = quant_ticks
)  

# Need the peaks of the negative strand because of the sequencing protocol
quant_peaks <- genomation::readGeneric('../mCQuaRna/results/quant_RS/brain/raw_peaks/segemehl/peakachu/brain.forward.bed', 
                                       strand = 6, 
                                       meta.cols = list(names = 4))

quantPeaks <- AnnotationTrack(filter(quant_peaks, strand == "+", start > start_coord, end < end_coord, seqnames == chromosome),
                              name = 'Quant Peaks',
                              fontcolor.item = "white",
                              fill = "red",
                              col = "red",
                              size = 2)



Gviz::plotTracks(list(#ideoTrack,
    gtrack,
    teRegio,
    cagePeaks,
    quantPeaks,
    cage_young_track,
    cage_old_track,
    nanopore_reads,
    teInstance
),
add35 = TRUE,
labelPos = 'below',
scale = 0.25,
shape = 'box',
fontface.main = 'arial',
chromosome = chromosome,
from = start_coord,
to = end_coord,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")

