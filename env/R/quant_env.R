#' Quant-Seq Environment
#' 
#' Loads the Quant-Seq specific environment
#' - counts of peaks
#' - paths for:
#'  - results
#'  - deseq
#' @export
load_quant_seq_env <- function(){
    
    logmsg("→ Loading quant-seq environment...")
    
    # Count tables
    quant_counts <<- list(
        brain = 'results/quant_seq/brain/counts/brain_peak_counts.csv',
        blood = 'results/quant_seq/blood/counts/blood_peak_counts.csv',
        skin = 'results/quant_seq/skinII/counts/skinII_peak_counts.csv')    
    
    # Directories & Files
    quant_results_dir <<- 'results/quant_seq/'
    quant_deseq_dir <<- paste0(quant_results_dir, 'deseq2/')
    
    if (!dir.exists(quant_deseq_dir)) {
        dir.create(quant_deseq_dir, recursive = TRUE)
    }
    
    quant_dds <<- paste0(quant_deseq_dir, "dds_quant_segemehl.Rdata")
    quant_deseq_res <<- paste0(quant_deseq_dir, "deseq_quant_segemehl.Rdata")
    quant_deseq_res_csv <<- "deseq_quant_all.csv"
    
    
}


#' Loading the annotation of peaks as a gRange object
#'
#' @export
load_quant_peak_annotation <- function(){
    
    logmsg("→ Loading annotation of CAGE peaks as gRange object...")
    
    quant_beds <- list('brain' ='results/quant_seq/brain/raw_peaks/brain.quant_peaks.bed',
                       'skin' = 'results/quant_seq/skinII/raw_peaks/skin.quant_peaks.bed',
                       'blood' = 'results/quant_seq/blood/raw_peaks/blood.quant_peaks.bed')
    
    
    quantRanges <<- sapply(names(quant_beds), simplify = F, function(tissue){
        
        
        file <- quant_beds[[tissue]]
        if (file.exists(file)){
            return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
        }else{
            logmsg(paste('Raw quant peak annotation does not exist for', tissue, '. Please check quant_store_annotation.R'))
        }
        
    })
    
    base_path <- paste0(getwd(), "/env/R/load_data")
    
    annotation_files <- c("annotations.R")
    
    for (f in annotation_files){
        source(file.path(base_path, f))
    }
    
}