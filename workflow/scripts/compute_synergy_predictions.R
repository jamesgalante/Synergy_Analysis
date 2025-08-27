
save.image("compute_synergy_predictions.rda")

# This file is to compute predicted synergistic interactions
suppressPackageStartupMessages({
	library(tidyverse)
	library(data.table)
  source(paste0(snakemake@scriptdir, "/plotting_functions.R"))
})

# Load params and paired ENCODE predictions
p300_threshold <- snakemake@params$p300_threshold
EE_distance_threshold <- snakemake@params$EE_distance_threshold
EE_contact_threshold <- snakemake@params$EE_contact_threshold
ENCODE_Score_threshold <- snakemake@params$ENCODE_Score_threshold
paired_predictions_with_hic <- fread(snakemake@input$paired_predictions_with_hic)

# Filter for synergistic pairs
synergistic_pairs <- paired_predictions_with_hic %>%
  # --- P300 filtering: E1 high, E2 low ---
  filter(
    normalizedEP300_enh.Feature_element1 > p300_threshold,
    normalizedEP300_enh.Feature_element2 < p300_threshold
  ) %>%
  # --- Distance/contact filtering ---
  filter(
    distance_between_elements <= EE_distance_threshold |    # always keep short-range
      (distance_between_elements > EE_distance_threshold &  # if long-range:
         !is.na(e1_e2_obs_scale) &                         #   must have contact value
         e1_e2_obs_scale >= EE_contact_threshold)          #   and it must pass threshold
  )

# Plot stats
contact_statistics <- plot_contact_statistics(
  all_df = paired_predictions_with_hic,
  syn_df = synergistic_pairs,
  EE_contact_threshold = EE_contact_threshold,
  cell_type = snakemake@wildcards$cell_type
)

p300_statistics <- plot_p300_statistics(
  all_df = paired_predictions_with_hic,
  syn_df = synergistic_pairs,
  p300_threshold = p300_threshold,
  cell_type = snakemake@wildcards$cell_type,
  box_width = 0.4
)

distance_statistics <- plot_distance_statistics(
  all_df = paired_predictions_with_hic,
  syn_df = synergistic_pairs,
  cell_type = snakemake@wildcards$cell_type,
  max_rows_per_panel = 5000
)

# Save output
ggsave(plot = contact_statistics, filename = snakemake@output$contact_statistics, device = "pdf", width = 10, height = 4)
ggsave(plot = p300_statistics, filename = snakemake@output$p300_statistics, device = "pdf", width = 6, height = 4)
ggsave(plot = distance_statistics, filename = snakemake@output$distance_statistics, device = "pdf", width = 10, height = 4)
write_tsv(synergistic_pairs, snakemake@output$synergy_predictions)