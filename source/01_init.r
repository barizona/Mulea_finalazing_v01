# library(devtools)
# install_github("https://github.com/koralgooll/MulEA.git")
library(MulEA) # TODO: rename MulEA -> mulea
library(tidyverse)

#xxxxxxxxxxxxxxxxxxxxxx
# Data --------------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxx

# import example ontology
# import other ontology from a GMT file using read_gmt()
data(geneSet)
# TODO: rename geneSet -> example_ontology
# TODO: colnames in example_ontology should be changed:
# ontologyId -> ontology_id
# ontologyName -> ontology_name
# listOfValues -> list_of_elements

# TODO: import test set: df containing all genes in the expariment ordered by FC and BH corrected p-val -> can be applied for ORA and GSEA as well

# import selected genes for ORA
data(selectDf)

# import background gene set for ORA
data(poolDf)

#xxxxxxxxxxxxxxxxxxxxx
# ORA ---------------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxx

ora_model <- ora(
  gmt = geneSet,
  element_names = selectDf$select, 
  background_element_names = poolDf$background_element_names,
  p_value_adjustment_method = "PT",
  number_of_permutations = 1000
)

ora_model
# TODO: change 
# Slot "test":
#   function (setBasemodel) 
#   {
#     setBasedTestRes <- NULL 
# etc
# to
# "Over-representation analysis (ORA)"

ora_results <- run_test(ora_model)

# TODO: colnames should be changed:
# ontologyId -> ontology_id
# ontologyName -> ontology_name
# nrCommonGenesOntologySet -> nr_common_with_tested_elements
# nrCommonGenesOntologyBackground -> nr_common_with_backgound_elements
# pValue -> p_value
# adjustedPValue -> adjusted_p_value
# adjustedPValueEmpirical -> eFDR

#xxxxxxxxxxxxxxxxx
# * Plotting ----------------------------------------------------------------
#xxxxxxxxxxxxxxxxx

#xxxxxxxx
# ** Initialization ----
#xxxxxxxx
ora_reshaped_results <- reshape_results(
  model = ora_model, 
  model_results = ora_results, 
  p_value_type_colname='adjustedPValueEmpirical'
)

#xxxxxxxx
# ** Graph ----
#xxxxxxxx
plot_graph(
  reshaped_results = ora_reshaped_results,
  p_value_max_threshold = 1.00,
  p_value_type_colname = "adjustedPValueEmpirical"
)
# TODO: change font in the legend to default

#xxxxxxxx
# ** Barplot ----
#xxxxxxxx
plot_barplot(
  reshaped_results = ora_reshaped_results,
  p_value_max_threshold = 1.00,
  p_value_type_colname = "adjustedPValueEmpirical"
)

#xxxxxxxx
# ** Heatmap ----
#xxxxxxxx
plot_heatmap(
  reshaped_results=ora_reshaped_results,
  p_value_max_threshold=1.00,
  p_value_type_colname = 'adjustedPValueEmpirical'
)

#xxxxxxxxxxxxxxxxxxxxx
# GSEA ---------------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxx

# Define an S4 object of class RankedBasedTest, run ranked based test (Subramanian method) and reshape results:

# TODO: rename ranked_model to gsea_model

gsea_model <- gsea(
    gmt = geneSet,
    element_names = selectDf$select,
    element_scores = selectDf$score
  )

# TODO: rename ranked_results to gsea_results
gsea_results <- run_test(gsea_model)

gsea_results
# TODO: why there are NAs in the result tab?
# ontologyId                                         ontologyName nrCommonGenesOntologySet nrCommonGenesOntologyBackground    pValue adjustedPValue
# 1 ID:0000001                          "mitochondrion inheritance"                        1                               2 0.3384030       0.891635
# 2 ID:0000002                   "mitochondrial genome maintenance"                        3                               9 0.6237113       0.891635
# 3 ID:0000009             "alpha-1,6-mannosyltransferase activity"                       NA                               2        NA             NA
# 4 ID:0000010          "trans-hexaprenyltranstransferase activity"                        2                               4 0.6865149       0.891635
# 5 ID:0000012                         "single strand break repair"                        1                               4 0.8916350       0.891635
# 6 ID:0000014 "single-stranded DNA endodeoxyribonuclease activity"                        5                               6 0.2360061       0.891635
# 7 ID:0000015                  "phosphopyruvate hydratase complex"                        3                               5 0.8573883       0.891635

#xxxxxxxxxxxxxxxxx
# * Plotting ----------------------------------------------------------------
#xxxxxxxxxxxxxxxxx

#xxxxxxxx
# ** Initialization ----
#xxxxxxxx
# TODO: ranked_reshaped_results -> gsea_reshaped_results
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
