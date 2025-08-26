
# This file is to compute predicted synergistic interactions
suppressPackageStartupMessages({
	library(tidyverse)
	library(data.table)
})

p300_threshold <- snakemake@params$p300_threshold
EE_distance_threshold <- snakemake@params$EE_distance_threshold

paired_predictions_with_hic <- fread(snakemake@input$paired_predictions_with_hic)

synergistic_pairs <- paired_predictions_with_hic %>%
  # p300 gating (keep your asymmetry: e1 > thresh, e2 < thresh)
  filter(
    normalizedEP300_enh.Feature_element1 > p300_threshold,
    normalizedEP300_enh.Feature_element2 < p300_threshold
  ) %>%
  # keep short-range OR long-range with a real contact value
  filter(
    distance_between_elements <= EE_distance_threshold | 
      (distance_between_elements > EE_distance_threshold & is.finite(e1_e2_obs_scale))
  ) %>%
  # label synergy: short-range ⇒ TRUE; long-range ⇒ e1_e2_obs_scale > 40
  mutate(
    synergistic = distance_between_elements < EE_distance_threshold |
      (distance_between_elements >= EE_distance_threshold & e1_e2_obs_scale > 40)
  )