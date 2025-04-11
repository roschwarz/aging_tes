
# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_cage_seq_env()


# ------------------------------------------------------------------------------
# All CAGE-Peaks
# ------------------------------------------------------------------------------


data <- sapply(names(cage_counts),
               simplify = F,
               function(x) getFeatureCountTab(cage_counts[[x]], filter = F))

# store peak annotation
for (tissue in names(data)) {
    
    meta <-  data[[tissue]][['meta']]
    directory <- strsplit(cage_counts[[tissue]], "/")[[1]][3]
    
    # set the length column to zero to mimic the 5th column in a bed file 
    meta$Length <- 0 
    
    write.table(meta, file = paste0(cage_results_dir,
                                    directory,
                                    '/raw_peaks/',
                                    tissue,
                                    '_cage_peaks.bed'),
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t')
}
