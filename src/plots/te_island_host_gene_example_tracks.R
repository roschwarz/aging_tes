if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

##########################
# Gene and TE-Annotation #
##########################

te_region_file <-  paste0("data/shared/mm10_TE_region.bed")

teRegionRanges <- genomation::readGeneric(te_region_file, 
                                          strand = 6, 
                                          meta.cols = list(names = 4))

teRanges <- blackRcloud::getRanges("mmusculus", 102, "te", 1)

###################################################
# Coordinates of the window that you want to show #
###################################################

# Skint5 example
#                                                                                                                                                                                                                
# Skint5: chr4:113477891-113999503                                                                                                                                                                        
# TE_Cluster_1087370: chr4:113905037-113911698

# chr4:113476891-114000503                                                                                                                                                                        

# Window to show::

chro = 'chr4'
str = '-'
start_coord = 113477891 - 5000
end_coord =  113999503 + 5000

#start_coord_pcdhb9 = 37396504
#end_coord_pcdhhb9 = 37401086

gtrack <- GenomeAxisTrack()


trk_skint5 <- GeneRegionTrack(txdb, 
                               chromosome = chro, 
                               start = start_coord, 
                               end = end_coord,
                               transcriptAnnotation = "gene",
                               fill = "black",
                               size = 2,
                               name = 'Gene')

trk_skint5 <- trk_skint5[attributes(trk_skint5)$range$transcript == "ENSMUST00000169631.7"]
attributes(trk_skint5)$range$gene <- "Skint5"

trk_te_islands <- AnnotationTrack(teRegionRanges[teRegionRanges@elementMetadata$names %in% c("TE_Cluster_1087370")],
                                  name = 'te island',
                                  id = c("TE island 129043"),
                                  featureAnnotation = "id",
                                  fontcolor.item = "white",
                                  fill = "black",
                                  size = 2)


TEs_of_interest = data.frame(teRanges) %>% 
    filter(seqnames == chro , 
           strand == str,
           dplyr::between(start, start_coord, end_coord)) %>% 
    pull(te.id) %>% 
    unique()

trk_te_instances <- AnnotationTrack(teRanges[teRanges@elementMetadata$te.id %in% TEs_of_interest & !duplicated(teRanges@elementMetadata$te.id)],
                                    featureAnnotation = "id",
                                    fontcolor.item = "white",
                                    #id = ids,
                                    stacking = "dense",
                                    #fill = fam_col,
                                    size = 2)

#################
# Cage big wig  #
#################

cage_limits = c(0,50)
ticks =  c(0,10,20, 30, 40, 50)

cage_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinI.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord, '-', end_coord), "GRanges"))

trk_cage_old <- DataTrack(range = cage_old.bw,
                          type = "h",
                          name = 'old',
                          col = "darkgreen",
                          ylim = cage_limits,
                          yTicksAt = ticks)    

cage_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinI.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord, '-', end_coord), "GRanges"))

trk_cage_young <- DataTrack(range = cage_young.bw,
                            type = "h",
                            name = 'young',
                            col = "darkgreen", 
                            yTicksAt = ticks,
                            ylim = cage_limits)    



Gviz::plotTracks(list(#ideoTrack,
    gtrack, 
    trk_skint5,
    trk_te_islands
),
shape = 'box',
chromosome = chro,
from = start_coord,
to = end_coord,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")


#######################
# Beginning of Skint5 #
#######################

chro = 'chr4'
str = '-'
start_coord_skin5_part = 113999503 - 5000
end_coord_skin5_part =  113999503 + 100


gtrack <- GenomeAxisTrack()


trk_skint5_beginning <- GeneRegionTrack(txdb, 
                              chromosome = chro, 
                              start = start_coord_skin5_part, 
                              end = end_coord_skin5_part,
                              transcriptAnnotation = "gene",
                              fill = "black",
                              size = 2,
                              name = 'Gene')

trk_skint5_beginning <- trk_skint5_beginning[attributes(trk_skint5_beginning)$range$transcript == "ENSMUST00000169631.7"]
attributes(trk_skint5_beginning)$range$gene <- "Skint5"


#################
# Cage big wig  #
#################

cage_limits = c(0,400)
ticks =  c(0, 200, 400)

cage_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinI.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord_skin5_part, '-', end_coord_skin5_part), "GRanges"))

trk_cage_old <- DataTrack(range = cage_old.bw,
                          type = "h",
                          name = 'old',
                          col = "darkgreen",
                          ylim = cage_limits,
                          yTicksAt = ticks)    

cage_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinI.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord_skin5_part, '-', end_coord_skin5_part), "GRanges"))

trk_cage_young <- DataTrack(range = cage_young.bw,
                            type = "h",
                            name = 'young',
                            col = "darkgreen", 
                            yTicksAt = ticks,
                            ylim = cage_limits)    

#################
# RNA-Seq track #
#################


rna_old.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/skinII.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord_skin5_part, '-', end_coord_skin5_part), "GRanges"))

rna_young.bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/skinII.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro , ':', start_coord_skin5_part, '-', end_coord_skin5_part), "GRanges"))

rna_limits = c(0,300)
rna_ticks =  c(0,50, 100, 150, 200, 250, 300)


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



Gviz::plotTracks(list(#ideoTrack,
    gtrack, 
    trk_skint5_beginning,
    trk_cage_young,
    trk_cage_old,
    trk_rna_young,
    trk_rna_old
),
shape = 'box',
chromosome = chro,
from = start_coord_skin5_part,
to = end_coord_skin5_part,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")


#####################
# TE island 1087370 #
#####################

chro = 'chr4'
str = '-'
start_coord_te_island = 113905037 - 500
end_coord_te_island  =  113911698 + 500


gtrack <- GenomeAxisTrack()


trk_skint5_te_island <- GeneRegionTrack(txdb, 
                              chromosome = chro, 
                              start = start_coord_te_island, 
                              end = end_coord_te_island,
                              transcriptAnnotation = "gene",
                              fill = "black",
                              size = 2,
                              name = 'Gene')

trk_skint5_beginning <- trk_skint5_beginning[attributes(trk_skint5_beginning)$range$transcript == "ENSMUST00000169631.7"]
attributes(trk_skint5_beginning)$range$gene <- "Skint5"


#################
# Cage big wig  #
#################

cage_old_te_island_bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinI.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord_te_island, '-', end_coord_te_island), "GRanges"))

cage_young_te_island_bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/cage/mm10/skinI.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord_te_island, '-', end_coord_te_island), "GRanges"))

cage_limits = c(0,40)
ticks =  c(0, 20, 40)

trk_cage_old_te_island <- DataTrack(range = cage_old_te_island_bw,
                          type = "h",
                          name = 'old',
                          col = "darkgreen",
                          ylim = cage_limits,
                          yTicksAt = ticks)    

trk_cage_young_te_island <- DataTrack(range = cage_young_te_island_bw,
                            type = "h",
                            name = 'young',
                            col = "darkgreen", 
                            yTicksAt = ticks,
                            ylim = cage_limits)    

#################
# RNA-Seq track #
#################


rna_old_te_island_bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/skinII.old.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro, ':', start_coord_te_island, '-', end_coord_te_island), "GRanges"))

rna_young_te_island_bw <- rtracklayer::import.bw('/misc/paras/data/www/robert/mCQuaRna/rna/mm10/skinII.young.reverse.bw') %>% 
    plyranges::filter_by_overlaps(as(paste0(chro , ':', start_coord_te_island, '-', end_coord_te_island), "GRanges"))

rna_limits = c(0,100)
rna_ticks =  c(0,50, 100)


trk_rna_old_te_island <- DataTrack(range = rna_old_te_island_bw, 
                         type = "hist", 
                         name = 'old',
                         col = "black",
                         window = -1,
                         fill.histogram = "black",
                         col.boxplotFrame = "black",
                         ylim = rna_limits,
                         yTicksAt = rna_ticks
)    



trk_rna_young_te_island <- DataTrack(
    range = rna_young_te_island_bw,
    type = "hist",
    name = 'young',
    col = "black",
    fill.histogram = "black",
    col.boxplotFrame = "black",
    widnow = -1,
    ylim = rna_limits,
    yTicksAt = rna_ticks
) 



Gviz::plotTracks(list(#ideoTrack,
    gtrack, 
    trk_skint5,
    trk_te_islands,
    trk_cage_young_te_island,
    trk_cage_old_te_island,
    trk_rna_young_te_island,
    trk_rna_old_te_island
),
shape = 'box',
chromosome = chro,
from = start_coord_te_island,
to = end_coord_te_island,
col.axis = "black",
col.title = "black",
cex.feature = 0.7,
background.title = "white",
fontcolor.legend = "black")

