#' @export
load_te_island_env <- function(){
    
    logmsg("→ Loading te island environment...")
    
    source("src/lib/te_island_nanopore_genomebrowser.R")
    
    # ==============================================================================
    # Directories
    # ==============================================================================
    
    tables_and_co <<- "results/te_island/tables_and_co"
    
    if ( !dir.exists(tables_and_co) )
        dir.create(tables_and_co, recursive = TRUE)
    
    # ==============================================================================
    # Annotations
    # ==============================================================================
    te_island_file_5_prime_extended <<- "data/processed/annotation/mm10_TE_island_5prime_extended.bed"
    brain_intergenic_te_islands <<- read.csv('results/te_island/long_read_verification/dataset_1_PMC10862843/intergenic.teilands.bed',
                                      sep = '\t',
                                      header = FALSE)
    
    # ==============================================================================
    # indie TE island annotation tables
    # TE islands that overlap with at least one cage peak and where at least on 
    # transposable element is considered as expressed (p adjusted is not NA)
    # ==============================================================================
    indie_te_island_bed <<- list(brain = "./results/te_island/brain_downsampled_indie_te_island.bed",
                                 skin = "./results/te_island/skin_indie_te_island.bed",
                                 blood = "./results/te_island/blood_indie_te_island.bed")
    
    
    # ==============================================================================
    # TE island transcripts s
    # TE islands that overlap with at least one cage- and quant-peak and where at 
    # least one transposable element is considered as expressed (p adjusted is not NA)
    # The script to annotate them can be found in src/te_island_transcripts.R
    # ==============================================================================
    
    te_island_transcripts <<- "results/TEItx/canonical_te_island_transcripts.csv"
    
    # ==============================================================================
    # Verification of TE islands by nanopore
    # intergenic TEs in brain
    # ==============================================================================
    brain_counts <<- read.csv('results/te_island/long_read_verification/dataset_1_PMC10862843/count_brain/quant.sf',
                       sep = '\t')
    
    bam_file_dataset_1 <<- 'results/te_island/long_read_verification/dataset_1_PMC10862843/SRR24578270.bam'
    
    cage_fwd_bw <<- 'results/te_island/brain_forward.bw'
    cage_rev_bw <<- 'results/te_island/brain_reverse.bw'
    cage_old_rev_bw <<- '/misc/paras/data/www/robert/mCQuaRna/cage/mm10/brain.old.reverse.bw'
    
    
    
}
