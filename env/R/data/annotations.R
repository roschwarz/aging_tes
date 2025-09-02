###################
# gene annotation #
###################

#' @export
load_gene_ranges <- function(){
    
    if (!exists("geneRanges")) {
        geneRanges <<- blackRcloud::getGeneRanges(type = 'gene')
        gene.exonRanges <<- blackRcloud::getGeneRanges(type = "exon")
    }
    
}

###################
# te annotation #
###################

#' @export
load_te_ranges <- function(){
    if (!exists("teRanges")) {
        tryCatch(
            {
                logmsg('Load te ranges...')
                teRanges <<- getRanges("mmusculus", 102, "te", 1)
            },
            error = function(cond) {
                logmsg('The te annotation does not appear to exists. Please run Rscript src/annotation/build_TE_annotation')
            }
        )
    }
    
}

#' @export
load_te_annotation <- function(){
    if (!file.exists(paste0(common_data, "te.annotation.Rdata"))) {
    
        te_annotation <<- transGrange(teRanges) %>%
            dplyr::rename(te_id = te.id) %>%
            filter(!duplicated(te_id)) %>%
            dplyr::select("te_id", "position", "ensembl_gene_id", "external_gene_name") %>%
            splitTEID("te_id")
    
    
        save(te_annotation, file = paste0(common_data, "te.annotation.Rdata"))
    
    }else{
    
        te_annotation <<- loadRdata(paste0(common_data, "te.annotation.Rdata"))
    
    }
}
# =============================== TE island ====================================
#
# Create and load TE island annotations. TEs within 500 bps are merged to TE
# islands, which are extended by 500 bp in either 5'-prime, 3'-prime, or both (full)
# directions. When the annotation files are not available than they are created 
# with the respective shell script.

#' @export
load_te_island_annotation <- function(){
    
    chrom_size <- paste0(common_data, "GRCm38.p6.chrom.sizes")
    # te_island_file <- paste0(common_data, "mm10_TE_island_5prime_extended.bed") # get rid of this file, however, you have to adapt the whole code base

    # Updated version of TE region files for the integration of Quant data
    te_island_file = 'data/shared/mm10_TE_island.bed'
    te_island_5prime_file = 'data/shared/mm10_TE_island_5prime_extended.bed'
    te_island_3prime_file = 'data/shared/mm10_TE_island_3prime_extended.bed'

    te_annotation_file <- paste0(common_data, "mm10_transposome_sorted.bed")


    if (!file.exists(te_island_file)) {
    
        system(paste("bash ./src/annotation/build_TE_island_annotation.sh",
                     te_annotation_file, 
                     chrom_size,
                     common_data), 
               wait = TRUE)
    
        te_islandRanges <<- readGeneric(te_island_file, 
                                      strand = 6, 
                                      meta.cols = list(names = 4))
    
    
        te_island_5primeRanges <<-
            readGeneric(te_island_5prime_file,
                        strand = 6,
                        meta.cols = list(names = 4)) 
    
        te_island_3primeRanges <<-
            readGeneric(te_island_3prime_file,
                    strand = 6,
                    meta.cols = list(names = 4)) 
    
    }else{
    
        te_islandRanges <<- readGeneric(te_island_file,
                                      strand = 6,
                                      meta.cols = list(names = 4)) 
        
        
        te_island_5primeRanges <<- readGeneric(te_island_5prime_file,
                        strand = 6,
                        meta.cols = list(names = 4)) 
        
        te_island_3primeRanges <<- readGeneric(te_island_3prime_file,
                        strand = 6,
                        meta.cols = list(names = 4)) 
    }
    
}




#' @export
load_te_island_instance_annotation <- function(te_island_file = 'data/shared/mm10_TE_island_5prime_extended.bed'){
    
    logmsg(paste0("Load annotation of instances overlapping with te island. File: ", te_island_file))
    
    te_island_anno <- readGeneric(te_island_file,
                                              strand = 6,
                                              meta.cols = list(names = 4))
    
    
    te_island_instances <- blackRcloud::intersectGranger(te_island_anno, teRanges, tab = 'all') 
    
    # TEs can intersect with multiple genes, so that they occur multiple times in                                           
    # the teRange object. Therefore a unique is applied at the end of the pipeline                                          
    te_island_instances <<- te_island_instances %>%                                                                          
        dplyr::select(query.names, subject.te.id, subject.position) %>%                                                     
        filter(!duplicated(subject.te.id)) %>%                                                                              
        unique()                                                                                                            
    
    #needs te_region_id as header for blackRcloud functions
    names(te_island_instances) <- c("te_island_id", "te_id", "te_position")   

}

#' @export
te_island_gene_rel <- function(){
    
    te_island_gene_relation <- intersectGranger(te_island_5primeRanges, 
                                                geneRanges, 
                                                "all") %>% 
        dplyr::select("query.names", # TE region id 
                      "subject.ensembl_gene_id", 
                      "subject.external_gene_name")
    
    names(te_island_gene_relation) <- c("te_island_id",
                                        "ensembl_gene_id",
                                        "external_gene_name")
    
    return(te_island_gene_relation)
}