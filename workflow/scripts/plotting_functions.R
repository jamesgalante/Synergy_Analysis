# ===============================
# Plotting functions
# ===============================

# CONTACT histogram
plot_contact_statistics <- function(all_df,
                                    syn_df = NULL,
                                    EE_contact_threshold,
                                    cell_type = NULL) {
  df_all <- all_df %>% transmute(group = "All", e1_e2_obs_scale = e1_e2_obs_scale)
  df <- if (!is.null(syn_df)) {
    df_syn <- syn_df %>% transmute(group = "Synergistic", e1_e2_obs_scale = e1_e2_obs_scale)
    bind_rows(df_all, df_syn)
  } else df_all
  
  df <- df %>% filter(!is.na(e1_e2_obs_scale), is.finite(e1_e2_obs_scale), e1_e2_obs_scale > 0)
  
  # n counts for strip labels
  n_counts <- df %>%
    group_by(group) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(group_label = paste0(group, " (n = ", format(n, big.mark = ","), ")"))
  df <- df %>% left_join(n_counts %>% select(group, group_label), by = "group")
  
  # Breaks + threshold
  all_vals   <- df$e1_e2_obs_scale
  base_breaks <- if (length(all_vals) > 0) scales::log_breaks()(range(all_vals)) else NULL
  thr         <- suppressWarnings(as.numeric(EE_contact_threshold))
  add_thr     <- if (is.finite(thr) && thr > 0) thr else NULL
  
  # Explicit ticks: add 0.01 and 1
  extra_ticks <- c(0.01, 1, 100)
  
  # Custom labeler: decimals < 1 keep 3 decimals, >= 1 no decimals
  tick_labeler <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return(NA_character_)
      if (v < 1) {
        format(v, trim = TRUE, scientific = FALSE, nsmall = 3)
      } else {
        format(v, trim = TRUE, scientific = FALSE, big.mark = "", nsmall = 0)
      }
    })
  }
  
  ggplot(df, aes(x = e1_e2_obs_scale)) +
    geom_histogram(color = "black", fill = "lightblue", bins = 75) +
    geom_vline(xintercept = thr, color = "red", linetype = "dashed") +
    scale_x_log10(
      breaks = sort(unique(c(base_breaks, add_thr, extra_ticks))),
      labels = tick_labeler
    ) +
    theme_classic() +
    facet_wrap(~ group_label, nrow = 1, scales = "free_y") +
    labs(
      title = paste0(toupper(cell_type), ": Histogram of E1-E2 contact values"),
      subtitle = "May include duplicates if both Enhancers target the same gene",
      x = "Contact (observed; scale-normalized)",
      y = "Count"
    ) +
    theme(
      plot.title = element_text(size = rel(1.3)),
      axis.title = element_text(size = rel(1.2)),
      axis.text  = element_text(size = rel(1.0)),
      strip.text = element_text(size = rel(1.0))
    )
}

# P300 boxplots
plot_p300_statistics <- function(all_df,
                                 syn_df = NULL,
                                 p300_threshold,
                                 cell_type = NULL,
                                 box_width = 0.4) {
  to_long <- function(df, group_name) {
    df %>%
      transmute(group = group_name,
                E1 = normalizedEP300_enh.Feature_element1,
                E2 = normalizedEP300_enh.Feature_element2) %>%
      tidyr::pivot_longer(cols = c(E1, E2), names_to = "Element", values_to = "P300")
  }
  
  df_all <- to_long(all_df, "All")
  df_long <- if (!is.null(syn_df)) bind_rows(df_all, to_long(syn_df, "Synergistic")) else df_all
  df_pos <- df_long %>% filter(is.finite(P300), P300 > 0)
  
  # n counts for strip labels (use rows of original data, not doubled)
  n_all <- nrow(all_df)
  n_syn <- if (!is.null(syn_df)) nrow(syn_df) else NA_integer_
  strip_map <- tibble::tibble(
    group = c("All", if (!is.null(syn_df)) "Synergistic"),
    group_label = c(
      paste0("All (n = ", format(n_all, big.mark = ","), ")"),
      if (!is.null(syn_df)) paste0("Synergistic (n = ", format(n_syn, big.mark = ","), ")")
    )
  )
  df_pos <- df_pos %>% left_join(strip_map, by = "group")
  
  base_breaks <- if (nrow(df_pos) > 0) scales::log_breaks()(range(df_pos$P300)) else NULL
  thr <- suppressWarnings(as.numeric(p300_threshold))
  add_thr <- if (is.finite(thr) && thr > 0) thr else NULL
  
  ggplot(df_pos, aes(x = Element, y = P300)) +
    geom_boxplot(fill = "lightblue", outlier.shape = NA, width = box_width) +
    geom_hline(yintercept = thr, color = "red", linetype = "dashed") +
    scale_y_log10(
      breaks = sort(unique(c(base_breaks, add_thr))),
      labels = scales::label_number()
    ) +
    theme_classic() +
    facet_wrap(~ group_label, nrow = 1) +
    labs(
      title = paste0(toupper(cell_type), ": P300 distribution by element"),
      subtitle = "May include duplicates",
      x = "Element",
      y = "Normalized P300"
    ) +
    theme(
      plot.title = element_text(size = rel(1.3)),  # larger, plain
      axis.title = element_text(size = rel(1.2)),  # larger axis titles
      axis.text  = element_text(size = rel(1.0)),  # default tick size
      strip.text = element_text(size = rel(1.0))   # facet labels plain
    )
}