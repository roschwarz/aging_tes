# Load custom environment if not already loaded
if (!"aging_tes" %in% loadedNamespaces()) {
    devtools::load_all("env/")
}

aging_tes::load_quant_env()


data <- sapply(names(quant_counts),
               simplify = F,
               function(x) getFeatureCountTab(quant_counts[[x]], filter = F))

for (tissue in names(data)) {
    
    meta <-  data[[tissue]][['meta']]
    directory <- strsplit(quant_counts[[tissue]], "/")[[1]][3]
    
    # set the length column to zero to mimic the 5th column in a bed file 
    meta$Length <- 0 
    
    write.table(meta, file = paste0(quant_results_dir,
                                    directory,
                                    '/raw_peaks/',
                                    tissue,
                                    '_quant_peaks.bed'),
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t')
}
