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


# Find interesting candidates by looking at read counts
counts <- fread("./results/te_island/long_read_verification/dataset_1_PMC10862843/count_brain/quant.sf")

# Select top 20 candidates that are longer than 1500 bp
top_100 <- counts %>% 
    filter(Length > 1500, Name %in% brain_intergenic_te_islands$V4) %>%
    arrange(desc(TPM)) %>% 
    head(100)


top_100 <- merge(top_100, brain_intergenic_te_islands, by.x = 'Name', by.y = 'V4')

# add coordinated as column like chr1:1000-2000

top_100 <- top_100 %>%
    mutate(coordinates = paste0(V1, ":", V2, "-", V3)) #%>%
    #select(Name, Length, TPM, coordinates, everything())


# Iterate over top 20 candidates and create genome browser shots and store them in 
# results/te_island/long_read_verification/genome_browser_shots

for (island in top_100$Name) {
    print(island)
    te_island_info <- top_100 %>% filter(Name == island)
    
    te_island = GRanges(
        seqnames = te_island_info$V1,
        ranges = IRanges(start = te_island_info$V2 - 1000, end = te_island_info$V3 + 1000),
        strand = te_island_info$V6
    )
    
    
    
    outfile <- paste0("results/te_island/long_read_verification/genome_browser_shots/", island, ".png")
    
    tracks <- genomeBrowserNanoporeTracks(
        c(island),
        te = te_island,
        te_island_annotation_file = te_island_file_5_prime_extended,
        bam_file = bam_file_dataset_1
    )
    
    if (!is.null(tracks)) {
        png(outfile, width = 10, height = 6, units = "in", res = 300)
           plotGenomeBrowserNanopore(tracks, te_island)
        dev.off()
    } else {
        # Optionally message or log that no plot was created for this island
        message("No plot created for island index ", island)
    }
   
}

