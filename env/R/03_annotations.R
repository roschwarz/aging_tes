###################
# gene annotation #
###################

if (!exists("geneRanges")) {
    geneRanges <- blackRcloud::getGeneRanges()
    gene.exonRanges <- blackRcloud::getGeneRanges(type = "exon")
}

###################
# te annotation #
###################

if (!exists("teRanges")) {


  #buildTEannotation(bedFile = "../../common_data/mm10/mm10_transposome.bed",
  #                  organism = "mmusculus",
  #                  gene.version = 102,
  #                  te.version = 1)

    teRanges <- getRanges("mmusculus", 102, "te", 1)
}


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