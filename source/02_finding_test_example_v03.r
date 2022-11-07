library(MulEA)
library(tidyverse)
#xxxxxxxxxxxxxxxxxxxxxx
# Data --------------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxx

# 2 értelmes GO: GSE68107.top.table_control_vs_cipro_30h.tsv
# 1 értelmes GO: GSE68106.top.table.tsv

# E.col treated ciprofloxacyn vs control in IST
Geo2R_reult_tab <- read_tsv("input/GSE55662.top.table_wt_cipro.tsv")
# positive logFC means underexpression when cipro

# Reformatting the Geo2R result table
Geo2R_reult_tab <- Geo2R_reult_tab %>% 
  # renaming Gene.symbol to Gene.symbol.multiple
  rename(Gene.symbol.multiple = Gene.symbol) %>% 
  # extracting the first gene from the Gene.symbol coumn
  mutate(Gene.symbol = str_remove(string = Gene.symbol.multiple,
                                  pattern = "\\/.*")) %>% 
  # removing rows where Gene.symbol is NA
  filter(!is.na(Gene.symbol)) %>% 
  # ordering by logFC
  arrange(desc(logFC))

# available GMTs
GO_GMT <- read_gmt("/home/barizona/Dropbox/MulEA/Databases/GO/update/GO_Escherichia_coli_Marton.gmt")
# filter
GO_GMT_filtered <- filter_ontology(gmt = GO_GMT, min_nr_of_elements = 10, max_nr_of_elements = 500)
KEGG_GMT <- read_gmt("/home/barizona/Dropbox/MulEA/Databases/KEGG Pathways/08_2020/KEGG_Escherichia_coli_gene_symbol_Leila.gmt")
# Operon_GMT <- read_gmt("/home/barizona/Dropbox/MulEA/Databases/Operon/Operon_Marton_Escherichia_coli_kegg.gmt")
# not gene symbol

#xxxxxxxxxxxxxxxxxxxxx
# ORA GO ---------------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxx

# Sign. overexpressed genes
E.coli_test_element_names <- Geo2R_reult_tab %>% 
  filter(adj.P.Val < 0.05
         & logFC < 0) %>% 
  select(Gene.symbol) %>% 
  pull() %>% 
  unique() %>% 
  sort()

# background genes
E.coli_background_element_names <- Geo2R_reult_tab %>% 
  select(Gene.symbol) %>% 
  pull() %>% 
  unique() %>% 
  sort()

ora_model <- ora(
  gmt = GO_GMT_filtered,
  element_names = E.coli_test_element_names, 
  background_element_names = E.coli_background_element_names,
  p_value_adjustment_method = "PT",
  number_of_permutations = 10000
)

ora_results <- run_test(ora_model)

ora_results %>% 
  filter(adjustedPValueEmpirical < 0.05) %>% 
  nrow()
# 228 GO filtered
# 28 KEGG

#xxxxxxxxxxxxxxxxx
# * Plotting ----------------------------------------------------------------
#xxxxxxxxxxxxxxxxx

#xxxxxxxx
# ** Initialization ----
#xxxxxxxx
ora_reshaped_results <- reshape_results(
  model = ora_model, 
  model_results = ora_results, 
  # ontology_id_colname = "ontologyName",
  # ontology_element_colname = "genIdInOntology",
  p_value_type_colname='adjustedPValueEmpirical'
)

# adding ontologyName to the ora_reshaped_results
ora_reshaped_results <- ora_reshaped_results %>% 
  left_join(., select(ora_results, ontologyId, ontologyName),
            by = "ontologyId")

#xxxxxxxx
# ** Graph ----
#xxxxxxxx
plot_graph(
  reshaped_results = ora_reshaped_results,
  ontology_id_colname = "ontologyName",
  p_value_max_threshold = 0.005,
  p_value_type_colname = "adjustedPValueEmpirical"
)
# Error in `[.data.table`(ontologies, i, ontologyId) : 
# j (the 2nd argument inside [...]) is a single symbol but column name 'ontologyId' is not found. Perhaps you intended DT[, ..ontologyId]. This difference to data.frame is deliberate and explained in FAQ 1.1.
# TODO: solve the above error

# TODO: order of arguments

#xxxxxxxx
# ** Barplot ----
#xxxxxxxx
plot_barplot(
  reshaped_results = ora_reshaped_results,
  p_value_max_threshold = 0.05,
  p_value_type_colname = "adjustedPValueEmpirical"
)

#xxxxxxxx
# ** Heatmap ----
#xxxxxxxx
plot_heatmap(
  reshaped_results=ora_reshaped_results,
  p_value_max_threshold = 0.05,
  p_value_type_colname = 'adjustedPValueEmpirical'
)

#xxxxxxxxxxxxxxxxxxxxx
# GSEA ---------------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxx

# if there are duplicated Gene.symbols keep the first one only
Geo2R_reult_tab_filt <- Geo2R_reult_tab %>% 
  group_by(Gene.symbol) %>%
  arrange(adj.P.Val) %>%
  filter(row_number()==1) %>% 
  ungroup()
  
gsea_model <- gsea(
  gmt = GO_GMT_filtered,
  element_names = Geo2R_reult_tab_filt$Gene.symbol,
  element_scores = Geo2R_reult_tab_filt$logFC,
  element_score_type = "neg"
)

gsea_results <- run_test(gsea_model)

#xxxxxxxxxxxxxxxxx
# * Plotting ----------------------------------------------------------------
#xxxxxxxxxxxxxxxxx
# TODO: itt tartok (a KEGG uolyan, mint az ora-nál, mosz GO filtered van bent)
#xxxxxxxx
# ** Initialization ----
#xxxxxxxx
gsea_reshaped_results <- reshape_results(
  model = gsea_model, 
  model_results = gsea_results, 
  ontology_id_colname = 'ontologyId'
)

#xxxxxxxx
# ** Graph ----
#xxxxxxxx
plot_graph(
  reshaped_results = gsea_reshaped_results,
  p_value_max_threshold = 1.00
)
# TODO: change font in the legend to default

#xxxxxxxx
# ** Barplot ----
#xxxxxxxx
plot_barplot(
  reshaped_results = gsea_reshaped_results,
  p_value_max_threshold = 1.00
)
# TODO: is the scale color OK? seems like the legend and the scale doesn't match. I need to further test it with the values of gsea_reshaped_results_read_in

# classes of the original gsea_reshaped_results
gsea_reshaped_results %>% class()
# [1] "data.table" "data.frame"
gsea_reshaped_results$ontologyId %>% class()
# [1] "character"
gsea_reshaped_results$genIdInOntology %>% class()
# [1] "character"
gsea_reshaped_results$adjustedPValue %>% class()
# [1] "numeric"

write_tsv(gsea_reshaped_results, "tmp/gsea_reshaped_results.tsv")

# classes of the gsea_reshaped_results_read_in
gsea_reshaped_results_read_in <- read.delim("tmp/gsea_reshaped_results_v2.tsv", sep = "\t") %>% 
  data.table::data.table()
gsea_reshaped_results_read_in %>% class()
# [1] "data.table" "data.frame"
gsea_reshaped_results_read_in$ontologyId %>% class()
# [1] "character"
gsea_reshaped_results_read_in$genIdInOntology %>% class()
# [1] "character"
gsea_reshaped_results_read_in$adjustedPValue %>% class()
# [1] "numeric"

plot_barplot(
  reshaped_results = gsea_reshaped_results_read_in,
  p_value_max_threshold = 1.00
)
# Error in `levels<-`(`*tmp*`, value = as.character(levels)) : 
# factor level [5] is duplicated

#xxxxxxxx
# ** Heatmap ----
#xxxxxxxx
plot_heatmap(
  reshaped_results = gsea_reshaped_results,
  p_value_max_threshold = 1.00
)
