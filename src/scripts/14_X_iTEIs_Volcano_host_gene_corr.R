

# Load Environment
if ( !exists("ENVIRONMENT_LOADED") ) {
    
    source('./01_load_environment.R')
    
} else if ( !ENVIRONMENT_LOADED ) {
    
    source('./01_load_environment.R')
    
}


# ======================= Load autonomous TE regions ===========================

auto.TE.cluster.files <- list(brain = "data/shared/brain_downsampled_independent_TE_regions.bed",
                              blood = "data/shared//blood_independent_TE_regions.bed",
                              skin = "data/shared/skinII_independent_TE_regions.bed"
)

auto.TE.regions <- sapply(tissues, simplify = F, function(x){
    df <- read.csv(auto.TE.cluster.files[[x]], 
                   sep = '\t', 
                   header = F) %>% 
        dplyr::pull(V4)
})


# ===================== Load DESeq2 results =================================

te.region.deseq.results <- read.csv('results/tables/02_deseq_results_te_region.csv')
gene.region.deseq.results <- read.csv('results/tables/02_deseq_results_gene.csv')

