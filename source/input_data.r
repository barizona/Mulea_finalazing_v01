library(tidyverse)
library(magrittr)
#xxxxxxxxxxxxxxxxxxxxxx

# logFC data from GEO GSE7305 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7305
# 10 Ovary-Disease Endometrium vs 10 Normal Endometrium

# read the data
GSE7305_logFC_tab <- read_tsv("input/GSE7305.geo2R.result.logFC.table.tsv")

# clean the data
GSE7305_logFC_tab %<>% 
  # remove ID col
  select(-ID) %>% 
  # remove where Gene.symbol is empty
  drop_na(Gene.symbol) %>% 
  # remove where Gene.symbol contains ///
  filter(!grepl('///', Gene.symbol)) %>%
  # remove where Gene.symbol starts with LOC
  filter(!grepl('^LOC', Gene.symbol)) %>%
  # remove where Gene.symbol starts with LINC
  filter(!grepl('^LINC', Gene.symbol))

# calculate the mean logFC by Gene.symbol
GSE7305_logFC_tab %<>% 
  group_by(Gene.symbol) %>% 
  summarise(logFC = mean(logFC)) %>% 
  # arrange by logFC
  arrange(desc(logFC))
# where high logFC means overexpression in Ovary-Disease Endometrium

write_tsv(GSE7305_logFC_tab, "input/GSE7305.Gene.names.logFC.tsv")

#xxxxxxxxxxxxxxxxx
# TODO: read TFLink gmt (small) and test