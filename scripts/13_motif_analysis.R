# --------------------------------- Notes --------------------------------------
# 
# This script is used to do a motif analysis for TE-CAGE peaks of the different
# TE regions in the different quadrants.
#
# Output data (../data/cage; ../tables):
#

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    source('./01_load_environment.R')
} else if ( !ENVIRONMENT_LOADED ) {
    
    source('./01_load_environment.R')
    
}

# ----------------------- for all expressed TE regions -------------------------

# The autonomous TEs are intersected with TSS to define the promoter region, 
# size of 500 bp beginning at the TSS an going towards the 5'-direction.
# The motif search and extraction is done by motifAnalysis.sh

#system(paste('PATH=/gsc/biosw/src/HOMER/bin/:$PATH; bash ./motifAnalysis.sh', chrom.size))


# ------------------- for all TE regions with a host gene ----------------------


auto.TE.cluster.files <- list(brain = "data/shared/brain_downsampled_independent_TE_regions.bed",
                              blood = "data/shared/blood_independent_TE_regions.bed",
                              skin = "data/shared/skinII_independent_TE_regions.bed"
)

Gene_host_TE_regions <- sapply(names(auto.TE.cluster.files), simplify = F, function(x){
    
    auto.df <- readGeneric(auto.TE.cluster.files[[x]], strand = 6, meta.cols = list(names = 4))
    
    auto.gene.df <- intersectGranger(geneRanges, auto.df, 'all') %>% 
        dplyr::select(subject.seqnames, 
                      subject.start, 
                      subject.end, 
                      subject.names, 
                      subject.width, 
                      subject.strand) %>% 
        dplyr::rename(chromosome = subject.seqnames,
                      start = subject.start,
                      end = subject.end,
                      te_region_id = subject.names,
                      width = subject.width,
                      strand = subject.strand) %>% 
        filter(!duplicated(te_region_id))
    
    return(auto.gene.df)
    
})

sapply(tissues, function(x){

    writeBedFile(Gene_host_TE_regions[[x]], 
                 file.name = paste0(x, '.intergenic.TE.regions.bed'),
                 directory = '../data/motif/intergenic_TE_regions/')
    
})

# ----------------------- Run Homer for intergenic TE regions ------------------

ctss_files <- list(brain = 'results/cage_RS/brain_downsampled_gclipped/tss/brain_downsampled.tss.bed',
                   skin = 'results/cage_RS/skinII_gclipped/tss/skin.tss.bed',
                   blood = 'results/cage_RS/blood_gclipped/tss/blood.tss.bed')

TE_region_files <- list(brain = '../data/motif/intergenic_TE_regions/brain.intergenic.TE.regions.bed',
                   skin = '../data/motif/intergenic_TE_regions/skin.intergenic.TE.regions.bed',
                   blood = '../data/motif/intergenic_TE_regions/blood.intergenic.TE.regions.bed')

sapply(tissues, function(x){

    system(paste('PATH=/gsc/biosw/src/HOMER/bin/:$PATH; bash ./runHomer.sh', 
                 TE_region_files[[x]],
                 ctss_files[[x]],
                 x,
                 '../data/motif/intergenic_TE_regions/'),
           wait = F)
    
})

# -------------------------------- Quadrant I ----------------------------------

Gene_host_TE_correlation <- read.csv('../tables/10_auto_TE_Host_correlation.csv')

res_I <- '../data/motif/intergenic_TE_regions/quadrant_I/'

if(!dir.exists(res_I)){
    
    dir.create(res_I)
}

Gene_host_TE_regions_I <- sapply(names(auto.TE.cluster.files), simplify = F, function(x){
    
    auto.df <- read.csv(auto.TE.cluster.files[[x]], 
                        sep = '\t',
                        header = F)
    
    te_regions <-  Gene_host_TE_correlation %>% 
        filter(tissue == x, rna_quadrant == 1, !duplicated(te_region_id)) %>% 
        dplyr::pull(te_region_id)
    
    auto.df <- auto.df %>% filter(V4 %in% te_regions)  
    
    
    
    return(auto.df)
})

TE_region_files_I <- sapply(tissues, simplify = F, function(x){
    
    writeBedFile(Gene_host_TE_regions_I[[x]], 
                 file.name = paste0(x, '.intergenic.TE.regions.I.bed'),
                 directory = res_I)
    
    return(paste0(res_I, x, '.intergenic.TE.regions.I.bed'))
    
})

sapply(tissues, function(x){

    system(paste('PATH=/gsc/biosw/src/HOMER/bin/:$PATH; bash ./runHomer.sh', 
                 TE_region_files_I[[x]],
                 ctss_files[[x]],
                 x,
                 res_I),
           wait = F)
    
})

# -------------------------------- Quadrant II --------------------------------

Gene_host_TE_correlation <- read.csv('../tables/10_auto_TE_Host_correlation.csv')

res_II <- '../data/motif/intergenic_TE_regions/quadrant_II/'

if(!dir.exists(res_II)){
    dir.create(res_II)
}

Gene_host_TE_regions_II <- sapply(names(auto.TE.cluster.files), simplify = F, function(x){
    
    auto.df <- read.csv(auto.TE.cluster.files[[x]], 
                        sep = '\t',
                        header = F)
    
    te_regions <-  Gene_host_TE_correlation %>% 
        filter(tissue == x, rna_quadrant == 2, !duplicated(te_region_id)) %>% 
        dplyr::pull(te_region_id)
    
    auto.df <- auto.df %>% filter(V4 %in% te_regions)  
    
    return(auto.df)
})

TE_region_files_II <- sapply(tissues, simplify = F, function(x){
    
    writeBedFile(Gene_host_TE_regions_II[[x]], 
                 file.name = paste0(x, '.intergenic.TE.regions.II.bed'),
                 directory = res_II)
    
    return(paste0(res_II, x, '.intergenic.TE.regions.II.bed'))
    
})

sapply(tissues, function(x){

    system(paste('PATH=/gsc/biosw/src/HOMER/bin/:$PATH; bash ./runHomer.sh', 
                 TE_region_files_II[[x]],
                 ctss_files[[x]],
                 x,
                 res_II),
           wait = F)
    
})

# -------------------------------- Quadrant III --------------------------------

Gene_host_TE_correlation <- read.csv('../tables/10_auto_TE_Host_correlation.csv')

res_III <- '../data/motif/intergenic_TE_regions/quadrant_III/'

if(!dir.exists(res_III)){
    dir.create(res_III)
}

Gene_host_TE_regions_III <- sapply(names(auto.TE.cluster.files), simplify = F, function(x){
    
    auto.df <- read.csv(auto.TE.cluster.files[[x]], 
                        sep = '\t',
                        header = F)
    
    te_regions <-  Gene_host_TE_correlation %>% 
        filter(tissue == x, rna_quadrant == 3, !duplicated(te_region_id)) %>% 
        dplyr::pull(te_region_id)
    
    auto.df <- auto.df %>% filter(V4 %in% te_regions)  
    
    return(auto.df)
})

TE_region_files_III <- sapply(tissues, simplify = F, function(x){
    
    writeBedFile(Gene_host_TE_regions_III[[x]], 
                 file.name = paste0(x, '.intergenic.TE.regions.III.bed'),
                 directory = res_III)
    
    return(paste0(res_III, x, '.intergenic.TE.regions.III.bed'))
    
})

sapply(tissues, function(x){

    system(paste('PATH=/gsc/biosw/src/HOMER/bin/:$PATH; bash ./runHomer.sh', 
                 TE_region_files_III[[x]],
                 ctss_files[[x]],
                 x,
                 res_III),
           wait = F)
    
})

# -------------------------------- Quadrant IV --------------------------------

Gene_host_TE_correlation <- read.csv('../tables/10_auto_TE_Host_correlation.csv')

res_IV <- '../data/motif/intergenic_TE_regions/quadrant_IV/'

if(!dir.exists(res_IV)){
    dir.create(res_IV)
}

Gene_host_TE_regions_IV <- sapply(names(auto.TE.cluster.files), simplify = F, function(x){
    
    auto.df <- read.csv(auto.TE.cluster.files[[x]], 
                        sep = '\t',
                        header = F)
    
    te_regions <-  Gene_host_TE_correlation %>% 
        filter(tissue == x, rna_quadrant == 4, !duplicated(te_region_id)) %>% 
        dplyr::pull(te_region_id)
    
    auto.df <- auto.df %>% filter(V4 %in% te_regions)  
    
    return(auto.df)
})

TE_region_files_IV <- sapply(tissues, simplify = F, function(x){
    
    writeBedFile(Gene_host_TE_regions_IV[[x]], 
                 file.name = paste0(x, '.intergenic.TE.regions.IV.bed'),
                 directory = res_IV)
    
    return(paste0(res_IV, x, '.intergenic.TE.regions.IV.bed'))
    
})

sapply(tissues, function(x){

    system(paste('PATH=/gsc/biosw/src/HOMER/bin/:$PATH; bash ./runHomer.sh', 
                 TE_region_files_IV[[x]],
                 ctss_files[[x]],
                 x,
                 res_IV),
           wait = F)
    
})
