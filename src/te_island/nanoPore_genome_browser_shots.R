# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_te_island_env()

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

aging_tes::load_annotations()
aging_tes::load_te_island_env()

load_te_ranges()

# Find interesting candidates by looking at read counts
counts <- fread("./results/te_island/long_read_verification/dataset_1_PMC10862843/count_brain/quant.sf")


####################
# TE_island_305170 #
####################

te_island_305170 <- counts %>% 
    filter(Name == "TE_Cluster_305170")

te_island_305170 <- merge(te_island_305170, brain_intergenic_te_islands, by.x = 'Name', by.y = 'V4')

trk_te_island_305170 = GRanges(
    seqnames = te_island_305170$V1,
    ranges = IRanges(start = te_island_305170$V2 - 1000, end = te_island_305170$V3 + 1000),
    strand = te_island_305170$V6
)

trk_te_island <- te_island_track(c("TE_Cluster_305170"),
                                    te_island_file_5_prime_extended)

trk_genome_axis <- GenomeAxisTrack()

trk_te_instances <- local_te_track(teRanges, trk_te_island_305170)

trk_nanopore <- nanopore_track(trk_te_island_305170,
                               bam_file = bam_file_dataset_1,
                               max_reads = 5000,
                               name = "Nanopore reads") 

trk_axis <- GenomeAxisTrack(scale = 0.25, labelPos = 'above', col = 'black')

trk_cage <- cage_seq_tracks(trk_te_island_305170)

te_island_305170_gbs <- Gviz::plotTracks(list(trk_axis,
                      trk_te_island,
                      trk_te_instances,
                      trk_cage,
                      trk_nanopore),
                 shape = 'box',
                 chromosome = "chr12",
                 from = te_island_305170$V2 - 1000,
                 to = te_island_305170$V3 + 1000,
                 col.axis = "black",
                 col.title = "black",
                 cex.feature = 0.7,
                 background.title = "white",
                 fontcolor.legend = "black")



cairo_pdf('results/figures/nanopore_example_te_island_305170.pdf', 
          height = 3,  width = 7)

Gviz::plotTracks(list(trk_axis,
                      trk_te_island,
                      trk_te_instances,
                      trk_cage,
                      trk_nanopore),
                 shape = 'box',
                 chromosome = "chr12",
                 from = te_island_305170$V2 - 1000,
                 to = te_island_305170$V3 + 1000,
                 col.axis = "black",
                 col.title = "black",
                 cex.feature = 0.7,
                 background.title = "white",
                 fontcolor.legend = "black")

dev.off()

#####################
# TE_island_1236382 #
#####################


te_island_1236382 <- counts %>% 
    filter(Name == "TE_Cluster_1236382")

te_island_1236382 <- merge(te_island_1236382, brain_intergenic_te_islands, by.x = 'Name', by.y = 'V4')

trk_te_island_1236382 = GRanges(
    seqnames = te_island_1236382$V1,
    ranges = IRanges(start = te_island_1236382$V2 - 1500, end = te_island_1236382$V3 + 5000),
    strand = te_island_1236382$V6
)

trk_te_island <- te_island_track(c("TE_Cluster_1236382"),
                                 te_island_file_5_prime_extended)

trk_genome_axis <- GenomeAxisTrack()

trk_te_instances <- local_te_track(teRanges, trk_te_island_1236382)


trk_nanopore <- nanopore_track(trk_te_island_1236382,
                               bam_file = bam_file_dataset_1,
                               max_reads = 10,
                               name = "Nanopore reads",
                               sampling_strat = "by_length") 


trk_cage <- cage_seq_tracks(trk_te_island_1236382)

trk_genes <- GeneRegionTrack(txdb, 
                               chromosome = 6, 
                               start = te_island_1236382$V2 - 1500, 
                               end = te_island_1236382$V3 + 5000,
                               strand = "+",
                               transcriptAnnotation = "gene",
                               fill = "black",
                               size = 2,
                               name = 'Gene')

trk_axis <- GenomeAxisTrack(scale = 0.2, labelPos = 'above', col = 'black')

te_island_1236382_gbs <- Gviz::plotTracks(list(trk_axis,
                      #trk_genome_axis,
                      trk_genes,
                      trk_te_island,
                      trk_te_instances,
                      trk_cage,
                      trk_nanopore),
                 shape = 'box',
                 chromosome = "chr6",
                 from = te_island_1236382$V2 - 5000,
                 to = te_island_1236382$V3 + 5000,
                 col.axis = "black",
                 col.title = "black",
                 cex.feature = 0.7,
                 background.title = "white",
                 fontcolor.legend = "black")

cairo_pdf('results/figures/nanopore_example_te_island_1236382.pdf', 
          height = 3,  width = 7)

te_island_1236382_gbs <- Gviz::plotTracks(list(trk_axis,
                                               #trk_genome_axis,
                                               trk_genes,
                                               trk_te_island,
                                               trk_te_instances,
                                               trk_cage,
                                               trk_nanopore),
                                          shape = 'box',
                                          chromosome = "chr6",
                                          from = te_island_1236382$V2 - 1500,
                                          to = te_island_1236382$V3 + 5000,
                                          col.axis = "black",
                                          col.title = "black",
                                          cex.feature = 0.7,
                                          background.title = "white",
                                          fontcolor.legend = "black")

dev.off()
