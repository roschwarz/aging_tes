# ------------------------------------------------------------------------------
# RNA-Seq 
# ------------------------------------------------------------------------------
#' @export
load_rna_seq_env <- function(){
    
    logmsg("→ Loading rna seq env...")
    
    # Count tables
    counts_rna <<- list(brain = "./results/rna_seq/brain/alignment_SalmonTE/EXPR.csv",
                        skin = "./results/rna_seq/skinII/alignment_SalmonTE/EXPR.csv",
                        blood = "./results/rna_seq/blood/alignment_SalmonTE/EXPR.csv")
    
    # Directories & Files
    rna_seq_results_dir <<- 'results/rna_seq/'
    rna_seq_deseq_dir <<- paste0(rna_seq_results_dir, 'deseq2/')
    
    deseq_dds_te <<- "dds_TE_instances_salmonTE.Rdata"
    deseq_results_te <<- "deseq_TE_instances_salmonTE.Rdata"
    deseq_results_te_csv <<- "02_deseq_results_te_instances.csv"
    
    deseq_dds_gene <<- "dds_genes_salmonTE.Rdata"
    deseq_results_gene <<- "deseq_genes_salmonTE.Rdata"
    deseq_results_gene_csv <<- "02_deseq_results_genes.csv"
    
    deseq_dds_mixed <<- "dds_mixed_salmonTE.Rdata"
    deseq_results_mixed <<- "deseq_mixed_salmonTE.Rdata"
    deseq_results_mixed_csv <<- "02_deseq_results_mixed.csv"
    
    base_path <- paste0(getwd(), "/env/R/data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    
}
