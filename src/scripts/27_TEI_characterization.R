library(data.table)
library(tidyverse)
library(blackRcloud)

setwd("/misc/paras/data/rschwarz/projects/mCQuaRna/src/scripts")

mm10_te_island_instances <- fread('../../data/shared/mm10_TE_region_instances.csv')

mm10_te_island_instances <- splitTEID(mm10_te_island_instances, 'te_id')

n_tei_member <- data.frame(table(mm10_te_island_instances$te_region_id))

data.frame(table(n_tei_member$Freq)) %>% view()

n_teis <- nrow(n_tei_member)

ggplot(n_tei_member %>% filter(Freq <= 30), aes(Freq)) +
    geom_bar() +
    geom_text(x = 20, y = 5e+05, label = paste0("n = ", n_teis)) +
    theme_bw()


