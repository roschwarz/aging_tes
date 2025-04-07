
###########
# general #
###########

load_project_fonts <- function() {
    message("Loading fonts...")
    extrafont::loadfonts()
}

create_project_dirs <- function() {
    dirs <- c("data/processed",
              "results/figures",
              "results/tables",
              "data/shared")
    for (d in dirs) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    }
}

common_data = 'data/shared/'
table_dir = 'results/tables/'
figure_dir = 'results/figures/'

###########
# RNA-Seq #
###########

counts_rna <- list(brain = "./results/rna_seq/brain/alignment_SalmonTE/EXPR.csv",
                   skin = "./results/rna_seq/skinII/alignment_SalmonTE/EXPR.csv",
                   blood = "./results/rna_seq/blood/alignment_SalmonTE/EXPR.csv")

rna_seq_results_dir = 'results/rna_seq/'
rna_seq_deseq_dir = paste0(rna_seq_results_dir, 'deseq2/')

deseq_dds = "dds_TE_instances_salmonTE.Rdata"
deseq_results = "deseq_TE_instances_salmonTE.Rdata"



