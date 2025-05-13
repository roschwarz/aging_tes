library(blackRcloud)
library(GenomicRanges)

buildTEannotation(bedFile = "data/shared/mm10_transposome.bed",
                  organism = "mmusculus",
                  gene.version = 102,
                  te.version = 1)