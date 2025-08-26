
save.image("pair_encode_predictions.rda")

# Load libraries
suppressPackageStartupMessages(
	library(tidyverse)
)

# Load in the predictions
encode_predictions <- read_tsv(snakemake@input$predictions)

# Pair up the elements
paired_predictions <- encode_predictions %>% 
  filter(distanceToTSS.Feature < snakemake@params$distance_threshold, isSelfPromoter == F) %>% 
  # Calculate signed distance to determine which side of gene each element is on
  mutate(
    element_center = (start + end) / 2, 
    signed_distance = element_center - TargetGeneTSS
  ) %>%
  # Create all possible pairs of elements targeting the same gene using self-join
  # This automatically duplicates all columns with _element1 and _element2 suffixes
  inner_join(., ., by = "TargetGene", suffix = c("_element1", "_element2")) %>%
  # Remove pairs where an element is paired with itself
  filter(name_element1 != name_element2) %>%
  # Only keep pairs where both elements are on the same side of the gene
  filter(sign(signed_distance_element1) == sign(signed_distance_element2)) %>%
  # Remove duplicate pairs by keeping only pairs where element1 is further from TSS
  # This ensures each unique pair appears only once (A-B, but not B-A)
  filter(distanceToTSS.Feature_element1 > distanceToTSS.Feature_element2) %>%
  # Add the pair distance
  mutate(distance_between_elements = distanceToTSS.Feature_element1 - distanceToTSS.Feature_element2)
  
write_tsv(paired_predictions, snakemake@output$paired_predictions)