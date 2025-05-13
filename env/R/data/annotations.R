###################
# gene annotation #
###################

#' @export
load_gene_ranges <- function(){
    
    if (!exists("geneRanges")) {
        geneRanges <- blackRcloud::getGeneRanges()
        gene.exonRanges <- blackRcloud::getGeneRanges(type = "exon")
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
    
        te.annotation <- transGrange(teRanges) %>%
            dplyr::rename(te_id = te.id) %>%
            filter(!duplicated(te_id)) %>%
            dplyr::select("te_id", "position", "ensembl_gene_id", "external_gene_name") %>%
            splitTEID("te_id")
    
    
        save(te.annotation, file = paste0(common_data, "te.annotation.Rdata"))
    
    }else{
    
        te.annotation <- loadRdata(paste0(common_data, "te.annotation.Rdata"))
    
    }
}
