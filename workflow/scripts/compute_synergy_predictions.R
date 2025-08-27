
# This file is to compute predicted synergistic interactions
suppressPackageStartupMessages({
	library(tidyverse)
	library(data.table)
})

# Load data
p300_threshold <- snakemake@params$p300_threshold
EE_distance_threshold <- snakemake@params$EE_distance_threshold
EE_contact_threshold <- snakemake@params$EE_contact_threshold
paired_predictions_with_hic <- fread(snakemake@input$paired_predictions_with_hic)

# Plot contact statistics
contact_statistics <- paired_predictions_with_hic %>%
  filter(!is.na(e1_e2_obs_scale)) %>%
  { 
    n_non_na <- nrow(.)
    ggplot(., aes(x = e1_e2_obs_scale)) +
      geom_histogram(color = "black", fill = "lightblue", bins = 100) +
      geom_vline(xintercept = as.integer(EE_contact_threshold), color = "red", linetype = "dashed") +
      scale_x_log10(
        breaks = c(scales::log_breaks()(range(.$e1_e2_obs_scale)), EE_contact_threshold),
        labels = scales::label_number()
      ) +
      theme_classic() +
      labs(
        title = "Histogram of non-NA E1-E2 contact values",
        subtitle = paste0(
          snakemake@wildcards$cell_type, 
          "; n = ", format(n_non_na, big.mark = ","), 
          "; May include duplicates if both Enhancers target more than 1 of the same gene"
        ),
        y = "Count",
        x = "Contact"
      )
  }

# Plot P300 statistics - HAVE TO UNDERSTAND WHAT DUPLICATES MEAN HERE
p300_statistics <- paired_predictions_with_hic %>%
  {
    n_total <- nrow(.)   # just use one n
    
    # ---- prepare y values (for breaks only) ----
    yvals <- c(.$normalizedEP300_enh.Feature_element1,
               .$normalizedEP300_enh.Feature_element2)
    yvals_pos <- yvals[is.finite(yvals) & yvals > 0]
    
    base_breaks <- if (length(yvals_pos) > 0) {
      scales::log_breaks()(range(yvals_pos))
    } else {
      NULL
    }
    
    thr <- suppressWarnings(as.numeric(p300_threshold))
    add_thr <- if (is.finite(thr) && thr > 0) thr else NULL
    
    # ---- filter JUST FOR PLOTTING (remove zeros for log transform) ----
    df_pos <- dplyr::filter(.,
                            normalizedEP300_enh.Feature_element1 > 0,
                            normalizedEP300_enh.Feature_element2 > 0
    )
    
    ggplot(df_pos, aes()) +
      geom_boxplot(aes(x = "E1", y = normalizedEP300_enh.Feature_element1),
                   fill = "lightblue", outlier.shape = NA, width = 0.6) +
      geom_boxplot(aes(x = "E2", y = normalizedEP300_enh.Feature_element2),
                   fill = "lightblue", outlier.shape = NA, width = 0.6) +
      geom_hline(yintercept = p300_threshold, color = "red", linetype = "dashed") +
      scale_y_log10(
        breaks = sort(unique(c(base_breaks, add_thr))),
        labels = scales::label_number()
      ) +
      theme_classic() +
      labs(
        title = "P300 distribution by element",
        subtitle = paste0(
          snakemake@wildcards$cell_type,
          "; n = ", format(n_total, big.mark = ","),
          "; May include duplicates"
        ),
        x = "Element",
        y = "Normalized P300"
      )
  }

# Filter for synergistic pairs
synergistic_pairs <- paired_predictions_with_hic %>%
  # p300 filtering 
  filter(
    normalizedEP300_enh.Feature_element1 > p300_threshold,
    normalizedEP300_enh.Feature_element2 < p300_threshold
  ) %>%
  # keep short-range OR long-range with a real contact value
  filter(
    distance_between_elements <= EE_distance_threshold | 
      (distance_between_elements > EE_distance_threshold & !is.na(e1_e2_obs_scale))
  ) %>%
  # filter for synergy: short-range ⇒ TRUE; long-range ⇒ e1_e2_obs_scale > 40
  filter(
    distance_between_elements < EE_distance_threshold | 
      (distance_between_elements >= EE_distance_threshold & e1_e2_obs_scale >= EE_contact_threshold)
  )

# Save output
ggsave(plot = contact_statistics, filename = snakemake@output$contact_statistics, device = "pdf", width = 8, height = 4)
ggsave(plot = p300_statistics, filename = snakemake@output$p300_statistics, device = "pdf", width = 4, height = 4)
write_tsv(synergistic_pairs, snakemake@output$synergy_predictions)