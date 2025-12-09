# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()

#################
# pcdhb15 locus #
#################

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


start_coord = 37468000
end_coord = 37476340

start_coord_pcdhb9 = 37396504
end_coord_pcdhhb9 = 37401086

gtrack <- GenomeAxisTrack()

##########################
# Gene and TE-Annotation #
##########################


te.region.file <- paste0("data/shared/mm10_TE_region.bed")

genes <- blackRcloud::getGeneRanges()

teRegionRanges <- genomation::readGeneric(te.region.file, 
                                          strand = 6, 
                                          meta.cols = list(names = 4))

teRanges <- blackRcloud::getRanges("mmusculus", 102, "te", 1)

pcdhb15_raw <- GeneRegionTrack(txdb, 
                               chromosome = 18, 
                               start = start_coord, 
                               end = end_coord,
                               transcriptAnnotation = "gene",
                               fill = "black",
                               size = 2,
                               name = 'Gene')

# get rid of the other annotations and keep only Pcdhb15 and rename it
pcdhb15_geneTRack <- pcdhb15_raw[attributes(pcdhb15_raw)$range$gene == "93886"]
attributes(pcdhb15_geneTRack)$range$gene <- "Pcdhb15" 


teRegio <- AnnotationTrack(teRegionRanges[teRegionRanges@elementMetadata$names %in% c("TE_Cluster_728983", "TE_Cluster_728986")],
                           name = 'te region',
                           id = c("TE Region 728239", "TE Region 728242"),
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



teInstance <- AnnotationTrack(teRanges[teRanges@elementMetadata$te.id %in% TEs_of_interest & !duplicated(teRanges@elementMetadata$te.id)],
                              featureAnnotation = "id",
                              fontcolor.item = "white",
                              id = ids,
                              stacking = "dense",
                              fill = fam_col,
                              size = 2)

##################################
# Panel Construction  Supplement #
##################################

######################
# TSS spanning reads #
######################

tss_spanning_reads <- Gviz::AnnotationTrack(
    range = 'results/pcdhb/protocadherin_cluster/pcdhb15_tss_reads.bam',
    genome = "mm10",
    name = "reads",
    fill = "lightgray",
    chromosome = "chr18",
    id = c()
)

##################
# Nanopore reads #
##################

nanopore_reads <- Gviz::AnnotationTrack(
    range = 'results/te_island/long_read_verification/dataset_1_PMC10862843/SRR24578270.bam',
    genome = "mm10",
    name = "nanopore reads",
    fill = "lightgray",
    chromosome = 'chr18',
    id = c()
)



##############
# Sanger Seq #
##############

sanger_primer <- Gviz::AnnotationTrack(GRanges(seqnames = c('chr18', 'chr18'),
                                               ranges = IRanges(start = c(37472806, 37473764),
                                                                end = c(37472829, 37473786),
                                                                names = c('chr18', 'chr18')),
                                               strand = c("+", "-"),
                                               name = c("Sanger_fwd", "Sanger_rev")),
                                       name = "Primer",
                                       fill = "darkgreen",
                                       fontcolor.legend = "green",
                                       size = 3,
                                       col = "black")

sanger_te_pcdhb15 <- Gviz::AnnotationTrack(GRanges(seqnames = c('chr18', 'chr18', 'chr18'),
                                                   ranges = IRanges(start = c(37472856, 37473189, 37473150),
                                                                    end = c(37473163, 37473729, 37473165),
                                                                    names = c('chr18', 'chr18', 'chr18')),
                                                   strand = c("+", "-", "-"),
                                                   name = c("Sanger_fwd", "Sanger_rev_block1", "Sanger_rev_block2")),
                                           name = "Sanger",
                                           fill = "darkgray",
                                           fontcolor.legend = "green",
                                           size = 3,
                                           col = "black")

#######################
# Highlight CAGE-peak #
#######################

ht <- HighlightTrack(trackList = c(nanopore_reads,
                                   tss_spanning_reads,
                                   sanger_te_pcdhb15),
                     start = c(37473513),
                     width = 154,
                     fill = 'transparent',
                     chromosome = 'chr18',
                     inBackground = FALSE)


Gviz::plotTracks(list(gtrack, 
                      pcdhb15_geneTRack,
                      teRegio,
                      teInstance,
                      #nanopore_reads,
                      ht,
                      sanger_primer),
                 shape = 'box',
                 chromosome = "chr18",
                 from = start_coord + 4000, # 4800
                 to = end_coord - 2500,
                 col.axis = "black",
                 col.title = "black",
                 cex.feature = 0.7,
                 background.title = "white",
                 fontcolor.legend = "black")
