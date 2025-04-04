
portionCalc <- function(df, object = 'super_family', counting = TRUE){
    
    # This function takes a data frame and counts for a specific objects and
    # calculates the percentage of categories. Additionally, the possibility 
    # for lengths exists .
    
    # To-Do:
    # - Check if a length column is in the data frame
    # - transfer to blackRCloud
    
    if(counting){
            
        portion.df <- dplyr::count(df, 
                                 !!as.symbol(object), 
                                 name = 'count')
        
        portion.df <- portion.df %>% 
            mutate(percentage = count/sum(count))
        
    }else{
        
        portion.df <- df %>% 
            group_by(!!as.symbol(object)) %>% 
            summarise(sum = sum(length))
    
        portion.df <- portion.df %>% 
            mutate(percentage = sum/sum(sum))
        
    }
    
    return(portion.df)
    
}


bino.pValue <- function(x,n,p){
    # two sided is the default
    binom.test(x,n,p)$p.value
}

prepEnrichment <- function(target, background, object = 'super_family'){
    
    names(target) <- paste0('target.', names(target)) 
    names(background) <- paste0('background.', names(background)) 
    
    
    df <- merge(target, 
                background, 
                by.x = paste0('target.', object), 
                by.y = paste0('background.', object), 
                all.x = T)
    
    df$ratio <- df$target.percentage/df$background.percentage
    
    # if("target.length" %in% names(df)) {
    #     df$length.ratio <- df$target.per.length/df$background.per.length
    # }
    #remove lines with NA
    
    names(df)[grepl(object, names(df))] <- object
    names(df)[grepl('tissue', names(df))] <- 'tissue'
    
    return(na.omit(df))
}


doEnrichment <- function(target, genome, object = 'super_family' ){
    
    
    df <- prepEnrichment(target, genome, object = object)
    
    variable = 'count' # if the data frame is based on counts, otherwise it can be base on a sum of length for example
    
    if(any(sapply(names(df), function(x) {grepl('sum', x)}))){
        variable = 'sum'
    }
    
    df$p.value <- mapply(bino.pValue, 
                         df[[paste0('target.', variable)]], 
                         sum(df[[paste0('target.', variable)]]), 
                         df$background.percentage)
    
    
    df$padj <- p.adjust(df$p.value, method = 'fdr')
    
    df$log10.padj <- -log10(df$padj)
    
    # If the p-value isn't detectable because it is to low than set it to 300
    df <- mutate(df, log10.padj = case_when(log10.padj == Inf ~ 0.0000001,
                                            log10.padj < 0.0000001 ~ 0.0000001,
                                            TRUE ~ log10.padj))
    
    # the log10.padj is set to negative for all depleted features
    # df <- mutate(df, log10.padj = case_when(ratio < 1 ~ -log10.padj,
    #                                         TRUE ~ log10.padj))
    
    
    # df <- mutate(df, ratio = case_when(n.target < 5 ~ NA_real_,
    #                                    TRUE ~ ratio),
    #              df, log10.padj = case_when(n.target < 5 ~ NA_real_,
    #                                         TRUE ~ log10.padj))
    
    return(df)
}

teCounter <- function(df, group = NULL, category = 'super_family'){
    
    if(!is.null(group)){
        group <- sym(group)
    }
    
    count <- df %>% 
        dplyr::count(!!group, !!sym(category), name = 'n.target')
    
    count <- count %>% 
        group_by(!!group) %>% 
        mutate(per.target = n.target/sum(n.target))
    
    count %>% filter(!!sym(category) != 'NA') 
}

TEEnrichment <- function(target, background, group = 'tissue'){
    
    groups <- levels(factor(target[['tissue']]))
    enrichment <- prepEnrichment(target, background)
    
    enrichment <- do.call('rbind', sapply(groups, simplify = F, function(x){
        
        df <- enrichment %>% filter(tissue == x)
        
        df$p.value <- mapply(bino.pValue, 
                             df$n.target, 
                             sum(df$n.target), 
                             df$per.genome)
        
        df$padj <- p.adjust(df$p.value)
        
        df$log10.padj <- -log10(df$padj)
        
        df <- mutate(df, log10.padj = case_when(log10.padj == Inf ~ 300,
                                                 TRUE ~ log10.padj))
        
        df <- mutate(df, log10.padj = case_when(ratio < 1 ~ -log10.padj,
                                                TRUE ~ log10.padj))
        
        df <- mutate(df, ratio = case_when(n.target < 5 ~ NA_real_,
                                           TRUE ~ ratio),
                     df, log10.padj = case_when(n.target < 5 ~ NA_real_,
                                                TRUE ~ log10.padj))
        
        df$label <- paste(round(df$ratio, 1), '(', df$n.target, ')')
        return(df)
    }))
    
    return(enrichment)
    
}



