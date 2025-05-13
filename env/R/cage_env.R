#' CAGE-Seq 
#'
#'Loads the CAGE-Specific environment
#' - counts of peaks
#' - paths for:
#'  - results
#'  - deseq
#'  
#' @export
load_cage_seq_env <- function(){
    
    logmsg("→ Loading  cage-seq environment...")
    
    # Count tables
    cage_counts <<- list(
        brain = 'results/cage_seq_gclipped/brain_downsampled_gclipped/counts/brain_downsampled_peak_counts.csv',
        blood = 'results/cage_seq_gclipped/blood_gclipped/counts/blood_peak_counts.csv',
        skin = 'results/cage_seq_gclipped/skinII_gclipped/counts/skinII_peak_counts.csv')    
    
    # Directories & Files
    cage_results_dir <<- 'results/cage_seq_gclipped/'
    cage_deseq_dir <<- paste0(cage_results_dir, 'deseq2/')
    
    if (!dir.exists(cage_deseq_dir)) {
        dir.create(cage_deseq_dir, recursive = TRUE)
    }
    
    cage_dds <<- paste0(cage_deseq_dir, "dds_cage_segemehl.Rdata")
    cage_deseq_res <<- paste0(cage_deseq_dir, "deseq_cage_segemehl.Rdata")
    cage_deseq_res_csv <<- "deseq_cage_all.csv"
    
    
}

#' Loading the annotation of peaks as a gRange object
#'
#' @export
load_cage_peak_annotation <- function(){
    
    logmsg("→ Loading annotation of CAGE peaks as gRange object...")
    
    raw_cage_counts <- list('brain' ='results/cage_seq_gclipped/brain_downsampled_gclipped/raw_peaks/brain_cage_peaks.bed',
                            'skin' = 'results/cage_seq_gclipped/skinII_gclipped/raw_peaks/skin_cage_peaks.bed',
                            'blood' = 'results/cage_seq_gclipped/blood_gclipped/raw_peaks/blood_cage_peaks.bed')
    
    
    cageRanges <<- sapply(names(raw_cage_counts), simplify = F, function(tissue){
        
        
        file <- raw_cage_counts[[tissue]]
        if (file.exists(file)){
            return(readGeneric(file, strand = 6, meta.cols = list(names = 4)))
        }else{
            logmsg(paste('Raw cage peak annotation does not exist for', tissue, '. Please check cage_store_annotation.R'))
        }
        
    })
    
}