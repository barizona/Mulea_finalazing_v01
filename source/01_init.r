# library(devtools)
# install_github("https://github.com/koralgooll/MulEA.git")
library(MulEA)
library(tidyverse)
library(magrittr)
#xxxxxxxxxxxxxxxxxxxxx

# import example gene set
# import other gene sets from a GMT file using read_gmt()
data(geneSet)
# TODO: rename: geneSet -> example_ontology
# TODO: colnames should be changed:
# ontologyId -> ontology_id
# ontologyName -> ontology_name
# listOfValues -> list_of_elements

# import "WHAT ARE WE IMPORTING HERE?"
data(selectDf)
# TODO: rename: selectDf -> example_tested_elements (only their names)

# import "WHAT ARE WE IMPORTING HERE?"
data(poolDf)
# TODO rename: poolDf -> example_background_elements (only their names)


ora_model <- ora(
  gmt = geneSet,
  element_names = selectDf$select, 
  background_element_names = poolDf$background_element_names,
  p_value_adjustment_method = "PT",
  number_of_permutations = 1000
)

ora_results <- run_test(ora_model)

ora_reshaped_results <- reshape_results(
  model = ora_model, 
  model_results = ora_results, 
  p_value_type_colname='adjustedPValueEmpirical'
)

# TODO: colnames should be changed:
# ontologyId -> ontology_id
# ontologyName -> ontology_name
# nrCommonGenesOntologySet -> nr_common_with_tested_elements
# nrCommonGenesOntologyBackground -> nr_common_with_backgound_elements
# pValue -> p_value
# adjustedPValue -> adjusted_p_value
# adjustedPValueEmpirical -> adjusted_p_value_empirical


