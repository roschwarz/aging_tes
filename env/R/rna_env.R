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

# ------------------------------------------------------------------------------
# RNA-Seq female
# ------------------------------------------------------------------------------
#' @export
load_rna_seq_female_env <- function(){
    
    logmsg("→ Loading rna seq env...")
    
    # Count tables
    counts_rna <<- list(brain = "./results/rna_seq/female/EXPR_brain.csv",
                        skin = "./results/rna_seq/female/EXPR_skin.csv")
    
    # Directories & Files
    rna_seq_results_dir <<- 'results/rna_seq/female/'
    rna_seq_deseq_dir <<- paste0(rna_seq_results_dir, 'deseq2/')
    
    deseq_dds_te <<- "dds_TE_instances_salmonTE.Rdata"
    deseq_results_te <<- "deseq_TE_instances_salmonTE.Rdata"
    deseq_results_te_csv <<- "02_deseq_results_te_female_instances.csv"
    
    deseq_dds_gene <<- "dds_genes_salmonTE.Rdata"
    deseq_results_gene <<- "deseq_genes_salmonTE.Rdata"
    deseq_results_gene_csv <<- "02_deseq_results_female_genes.csv"
    
    deseq_dds_mixed <<- "dds_mixed_salmonTE.Rdata"
    deseq_results_mixed <<- "deseq_mixed_salmonTE.Rdata"
    deseq_results_mixed_csv <<- "02_deseq_results_female_mixed.csv"
    
    base_path <- paste0(getwd(), "/env/R/data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    
}


# ------------------------------------------------------------------------------
# RNA-Seq public data set
# ------------------------------------------------------------------------------
#' @export
load_rna_seq_public_data_env <- function(){
    
    logmsg("→ Loading rna seq env for public data PMC11529319 (Liver, Gastrocnemius muscle [gm], white adipose tissue [wat])")
    
    # Public data from PMC11529319
    
    # Tissues:
    #   - liver
    #   - gm - Gastrocnemius muscle
    #   - wat - white adipose tissue
    
    # Count tables
    
    
    counts_rna <<- list(female = "./results/rna_seq/public/PMCID_PMC11529319/female/EXPR.csv",
                        male = "./results/rna_seq/public/PMCID_PMC11529319/male/EXPR.csv")
    
    meta_data <<- list(female = "./data/public/PMCID_PMC11529319/female/SraRunTable.csv",
                       male = "./data/public/PMCID_PMC11529319/male/SraRunTable.csv")
    
    # Directories & Files
    rna_seq_results_dir <<- 'results/rna_seq/public/PMCID_PMC11529319/'
    rna_seq_deseq_dir <<- paste0(rna_seq_results_dir, 'deseq2/')
    
    
    deseq_files <<- list(female = list(
        deseq_dds_te_liver = "dds_TE_instances_salmonTE_liver_female.Rdata",
        deseq_results_te_liver = "deseq_TE_instances_salmonTE_liver_female.Rdata",
        
        deseq_dds_gene_liver = "dds_genes_salmonTE_liver_female.Rdata",
        deseq_results_gene_liver = "deseq_genes_salmonTE_liver_female.Rdata",
        
        deseq_dds_te_gm = "dds_TE_instances_salmonTE_gm_female.Rdata",
        deseq_results_te_gm = "deseq_TE_instances_salmonTE_gm_female.Rdata",
        
        deseq_dds_gene_gm = "dds_genes_salmonTE_gm_female.Rdata",
        deseq_results_gene_gm = "deseq_genes_salmonTE_gm_female.Rdata",
        
        deseq_dds_te_wat = "dds_TE_instances_salmonTE_wat_female.Rdata",
        deseq_results_te_wat = "deseq_TE_instances_salmonTE_wat_female.Rdata",
        
        deseq_dds_gene_wat = "dds_genes_salmonTE_wat_female.Rdata",
        deseq_results_gene_wat = "deseq_genes_salmonTE_wat_female.Rdata"
        
    ),
    
    male = list(
        deseq_dds_te_liver = "dds_TE_instances_salmonTE_liver_male.Rdata",
        deseq_results_te_liver = "deseq_TE_instances_salmonTE_liver_male.Rdata",
        
        deseq_dds_gene_liver = "dds_genes_salmonTE_liver_male.Rdata",
        deseq_results_gene_liver = "deseq_genes_salmonTE_liver_male.Rdata",
        
        deseq_dds_te_gm = "dds_TE_instances_salmonTE_gm_male.Rdata",
        deseq_results_te_gm = "deseq_TE_instances_salmonTE_gm_male.Rdata",
        
        deseq_dds_gene_gm = "dds_genes_salmonTE_gm_male.Rdata",
        deseq_results_gene_gm = "deseq_genes_salmonTE_gm_male.Rdata",
        
        deseq_dds_te_wat = "dds_TE_instances_salmonTE_wat_male.Rdata",
        deseq_results_te_wat = "deseq_TE_instances_salmonTE_wat_male.Rdata",
        
        deseq_dds_gene_wat = "dds_genes_salmonTE_wat_male.Rdata",
        deseq_results_gene_wat = "deseq_genes_salmonTE_wat_male.Rdata"
        
        
        )
    )
    
    
    base_path <- paste0(getwd(), "/env/R/data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    