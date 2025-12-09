
#' @export
load_deseq_merged <- function(){
    if (!exists("deseq.te.merged")) {
        
        if (file.exists(paste0(table_dir, deseq_results_te_csv))) {
            message("Load Deseq results for te instances")
            deseq.te.merged <- read.csv(paste0(table_dir, deseq_results_te_csv))
        }else{
            message(paste(table_dir, deseq_results_te_csv, 'does not exists.',
                           'Check your DESeq analysis'))
        }
        
    }
    
    if (!exists("deseq.mixed.merged")) {
        
        if (file.exists(paste0(table_dir, deseq_results_mixed_csv))) {
            message("Load Deseq results for te instances and genes")
            deseq.mixed.merged <- read.csv(paste0(table_dir, deseq_results_mixed_csv))
        }else{
            message(paste(table_dir, deseq_results_te_csv, 'does not exists.',
                          'Check your DESeq analysis'))
        }
        
    }
    
}

#' @export
load_deseq_results <- function(){
    rna_deseq_res_te <<- blackRcloud::loadRdata(paste0(rna_seq_deseq_dir, deseq_results_te))
}