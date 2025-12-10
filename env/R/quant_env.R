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