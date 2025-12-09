# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

###################################################
# Coordinates of the window that you want to show #
###################################################

start_coord = 37467000
end_coord = 37476340

#start_coord_pcdhb9 = 37396504
#end_coord_pcdhhb9 = 37401086

gtrack <- GenomeAxisTrack()

##########################
# Gene and TE-Annotation #
##########################

te_region_file <-  paste0("data/shared/mm10_TE_region.bed")

teRegionRanges <- genomation::readGeneric(te_region_file, 
                                          strand = 6, 
                                          meta.cols = list(names = 4))

teRanges <- blackRcloud::getRanges("mmusculus", 102, "te", 1)

trk_pcdhb15 <- GeneRegionTrack(txdb, 
                               chromosome = 18, 
                               start = start_coord, 
                               end = end_coord,
                               transcriptAnnotation = "gene",
                               fill = "black",
                               size = 2,
                               name = 'Gene')

# get rid of the other annotations and keep only Pcdhb15 and rename it
trk_pcdhb15 <- trk_pcdhb15[attributes(trk_pcdhb15)$range$gene == "93886"]
attributes(trk_pcdhb15)$range$gene <- "Pcdhb15" 


trk_te_islands <- AnnotationTrack(teRegionRanges[teRegionRanges@elementMetadata$names %in% c("TE_Cluster_728983", "TE_Cluster_728986")],
                           name = 'te island',
                           id = c("TE island 728239", "TE island 728242"),
                           featureAnnotation = "id",
                           fontcolor.item = "white",
                           fill = "black",
                           size = 2)


TEs_of_interest = data.frame(teRanges) %>% 
    filter(seqnames == "chr18", 
           strand == "+",
           dplyr::between(start, start_coord, end_coord)) %>% 
    pull(te.id) %>% 
    unique()

ids = c("L1", "L1", "L1","L1", "B1", "L1")
fam_col = c(rep("darkred",4), "darkgray", "darkred")


trk_te_instances <- AnnotationTrack(teRanges[teRanges@elementMetadata$te.id %in% TEs_of_interest & !duplicated(teRanges@elementMetadata$te.id)],
                              featureAnnotation = "id",
                              fontcolor.item = "white",
                              id = ids,
                              stacking = "dense",
                              fill = fam_col,
                              size = 2)

##################
# Primer & siRNA #
##################

trk_primer_siRNA <- Gviz::AnnotationTrack(GRanges(seqnames = rep('chr18', 6),
                                          ranges = IRanges(start = c(37467484, 37469362, 37470788, 37474847, 37474202, 37469183 ),
                                                           end = c(37467590, 37469491, 37470887, 37474990, 37474222, 37469203),
                                                           names = rep('chr18', 6)),
                                          strand = rep('+', 6),
                                          name = c("TE_up", "TE_inter", "TE_down", "Pcdhb15", "siPcdhb15", 'siTE')),
                                  name = "Primer & Motifs",
                                  size = 2,
                                  fill = c("black", "black", "black", "black", "purple", "purple"))



#############
# Sox motif #
#############

trk_sox <- Gviz::AnnotationTrack(GRanges(seqnames = c('chr18'),
                                          ranges = IRanges(start = c(37468799),
                                                           end = c(37468809),
                                                           names = c('chr18')),
                                          strand = c("+"),
                                          name = c("Sox5")),
                                  name = "Sox Motif",
                                  size = 2,
                                  fill = c("white"))



#################
# Cage big wig  #
#################

cage_limits = c(0,300)
ticks =  c(0,100,200,300)

cage_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/brain.old.forward.down.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr18:', start_coord, '-', end_coord), "GRanges"))

trk_cage_old <- DataTrack(range = cage_old.bw,
                            type = "h",
                            name = 'old',
                            col = "darkgreen",
                            ylim = cage_limits,
                            yTicksAt = ticks)    

cage_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/brain.young.forward.down.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr18:', start_coord, '-', end_coord), "GRanges"))

trk_cage_young <- DataTrack(range = cage_young.bw,
                              type = "h",
                              name = 'young',
                              col = "darkgreen", 
                              yTicksAt = ticks,
                              ylim = cage_limits)    

#################
# RNA-Seq track #
#################

rna_limits = c(0,150)
rna_ticks =  c(0,50,100,150)

rna_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/brain.old.forward.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0('chr18:', start_coord, '-', end_coord), "GRanges"))


trk_rna_old <- DataTrack(range = rna_old.bw, 
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
    plyranges::filter_by_overlaps(as(paste0('chr18:', start_coord, '-', end_coord), "GRanges"))


trk_rna_young <- DataTrack(
    range = rna_young.bw,
    type = "hist",
    name = 'young',
    col = "black",
    fill.histogram = "black",
    col.boxplotFrame = "black",
    widnow = -1,
    ylim = rna_limits,
    yTicksAt = rna_ticks
) 


###################
# Overlay tracks  #
###################

trk_feature_anno <- OverlayTrack(trackList = list(trk_te_islands, trk_sox, trk_pcdhb15), name = 'features') 
############
# ideogram #
############

ideoTrack <- IdeogramTrack(genome = "mm10", chromosome = "chr18")

###############################
# Panel Construction  Panel_5 #
###############################


Gviz::plotTracks(list(#ideoTrack,
    gtrack, 
    trk_primer_siRNA,
    trk_feature_anno,
    trk_te_instances,
    trk_cage_young,
    trk_cage_old,
    trk_rna_young,
    trk_rna_old
),
shape = 'box',
chromosome = "chr18",
from = start_coord,
to = end_coord,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")


#########################
# Table of TE instances #
#########################

data.frame(teRanges) %>% 
    filter(seqnames == "chr18", 
           strand == "+",
           dplyr::between(start, start_coord, end_coord)) %>% 
    filter(!duplicated(te.id)) %>% 
    blackRcloud::splitTEID("te.id") %>% 
    dplyr::select(te.id, order, super_family, family, Kimura)
