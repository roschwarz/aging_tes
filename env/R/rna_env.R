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
    rna_seq_deseq_dir_female <<- paste0(rna_seq_results_dir, 'female/deseq2/')  
    
    deseq_dds_te <<- "dds_TE_instances_salmonTE.Rdata"
    deseq_results_te <<- "deseq_TE_instances_salmonTE.Rdata"
    deseq_results_te_csv <<- "02_deseq_results_te_instances.csv"
    deseq_results_te_csv_male <<- "02_deseq_results_te_instances.csv"
    deseq_results_te_csv_female <<- "02_deseq_results_te_female_instances.csv"
    
    deseq_dds_gene <<- "dds_genes_salmonTE.Rdata"
    deseq_results_gene <<- "deseq_genes_salmonTE.Rdata"
    deseq_results_gene_csv <<- "02_deseq_results_genes.csv"
    
    deseq_dds_mixed <<- "dds_mixed_salmonTE.Rdata"
    deseq_results_mixed <<- "deseq_mixed_salmonTE.Rdata"
    deseq_results_mixed_csv <<- "02_deseq_results_mixed.csv"
    
    base_path <- paste0(getwd(), "/env/R/load_data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    
}

# ------------------------------------------------------------------------------
# RNA-Seq female
# ------------------------------------------------------------------------------
#' @export
load_rna_seq_female_old_env <- function(){
    
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
    
    base_path <- paste0(getwd(), "/env/R/load_data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    
}


# ----------------------------------------------------------------------------
# All RNA-Seq experiments
#
# TODO: MAKE IT MORE GENERIC FOR THE WHOLE DATE WHEN YOU HAVE TIME
# ------------------------------------------------------------------------------
#@export
#' Load all RNA-Seq environments
load_all_rna_seq_envs <- function(sex, base_results_path = "./results/rna_seq"){
    
    logmsg("→ Loading all rna seq envs...")
    
    if (!sex %in% c("female", "male")) {
        stop("sex must be either 'female' or 'male'")
    }
    
    # New environment for data set
    env <- new.env(parent = .GlobalEnv)
    
    logmsg(paste0("→ Creating rna seq namespace for ", sex, "..."))
    
    # Directories & Files 
    if (sex == 'male') {
        env$count_tables <- list(brain = "./results/rna_seq/brain/alignment_SalmonTE/EXPR.csv",
                                 skin = "./results/rna_seq/skinII/alignment_SalmonTE/EXPR.csv",
                                 blood = "./results/rna_seq/blood/alignment_SalmonTE/EXPR.csv")
        
    } else if (sex == 'female') {
        env$count_tables <- "./results/rna_seq/female/detector/EXPRs.csv"
    }
    
    # Helper function for standard file names
    create_file_mapping <- function(sex, analysis_types = c("te", "gene", "mixed")) {
        file_mapping <- list()

        for (type in analysis_types) {
            base_name <- switch(type,
                                "te" = "TE_instances",
                                "gene" = "genes",
                                "mixed" = "mixed")

            file_mapping[[paste0("deseq_dds_", type)]] <-
                paste0("dds_", base_name, "_salmonTE.Rdata")

            file_mapping[[paste0("deseq_results_", type)]] <-
                paste0("deseq_", base_name, "_salmonTE.Rdata")

            file_mapping[[paste0("deseq_results_", type, "_csv")]] <-
                paste0("02_deseq_results_", sex, "_",
                       ifelse(type == "te", "instances",
                              ifelse(type == "gene", "genes", "mixed")), ".csv")
        }

        return(file_mapping)
    }
    
    env$results_dir <- paste0(base_results_path, "_", sex, "/")
    env$deseq_dir <- paste0(env$results_dir, "deseq2/")
    
    env$files <- create_file_mapping(sex)
    
    base_path <- paste0(getwd(), "/env/R/load_data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R", "annotations.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    
    return(env)
    
}



#' Generic RNA-Seq Environment Loader  
#' @param sex Character string indicating sex ("female" or "male")
#' @param base_results_path Character string for base results directory
#' @export
load_foo_rna_seq_env <- function(sex, base_results_path = "./results/rna_seq/") {
    
    # Validierung der Eingabe
    if (!sex %in% c("female", "male")) {
        stop("sex must be either 'female' or 'male'")
    }
    
    logmsg(paste0("→ Loading rna seq env for ", sex, "..."))
    
    # Helper function für standardisierte Dateinamen
    create_file_mapping <- function(sex, analysis_types = c("te", "gene", "mixed")) {
        file_mapping <- list()
        
        for (type in analysis_types) {
            base_name <- switch(type,
                                "te" = "TE_instances",
                                "gene" = "genes", 
                                "mixed" = "mixed")
            
            file_mapping[[paste0("deseq_dds_", type)]] <- 
                paste0("dds_", base_name, "_salmonTE.Rdata")
            
            file_mapping[[paste0("deseq_results_", type)]] <- 
                paste0("deseq_", base_name, "_salmonTE.Rdata")
            
            file_mapping[[paste0("deseq_results_", type, "_csv")]] <- 
                paste0("02_deseq_results_", sex, "_", 
                       ifelse(type == "te", "instances", 
                              ifelse(type == "gene", "genes", "mixed")), ".csv")
        }
        
        return(file_mapping)
    }
    
    # Pfade erstellen
    sex_results_dir <- paste0(base_results_path, sex, "/")
    sex_deseq_dir <- paste0(sex_results_dir, "deseq2/")
    
    # Globale Variablen zuweisen
    assign("counts_rna", paste0(sex_results_dir, "detector/EXPRs.csv"), envir = .GlobalEnv)
    assign("rna_seq_results_dir", sex_results_dir, envir = .GlobalEnv) 
    assign("rna_seq_deseq_dir", sex_deseq_dir, envir = .GlobalEnv)
    
    # Alle Datei-Mappings erstellen und zuweisen
    file_mappings <- create_file_mapping(sex)
    for (var_name in names(file_mappings)) {
        assign(var_name, file_mappings[[var_name]], envir = .GlobalEnv)
    }
    
    # Zusätzliche R-Dateien sourcen
    base_path <- paste0(getwd(), "/env/R/load_data")
    rna_files <- c("rna_seq_load_deseq_tes.R", "annotations.R")
    
    for (f in rna_files) {
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
    counts_rna <<- "./results/rna_seq/female/detector/EXPRs.csv"
    
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
    
    base_path <- paste0(getwd(), "/env/R/load_data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R", "annotations.R")
    
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
    
    
    base_path <- paste0(getwd(), "/env/R/load_data")
    
    rna_files <- c("rna_seq_load_deseq_tes.R")
    
    for (f in rna_files){
        source(file.path(base_path, f))
    }
    
}