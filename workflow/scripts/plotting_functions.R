# ===============================
# Plotting functions
# ===============================

# All vs Synergistic: CONTACT histogram
plot_contact_statistics <- function(all_df,
                                    syn_df = NULL,
                                    EE_contact_threshold,
                                    cell_type = NULL) {
  # Prepare combined data with a group column for faceting
  df_all <- all_df %>%
    transmute(group = "All", e1_e2_obs_scale = e1_e2_obs_scale)
  
  if (!is.null(syn_df)) {
    df_syn <- syn_df %>%
      transmute(group = "Synergistic", e1_e2_obs_scale = e1_e2_obs_scale)
    df <- bind_rows(df_all, df_syn)
  } else {
    df <- df_all
  }
  
  df <- df %>%
    filter(!is.na(e1_e2_obs_scale), is.finite(e1_e2_obs_scale), e1_e2_obs_scale > 0)
  
  # Counts per facet for subtitle
  n_counts <- df %>%
    group_by(group) %>%
    summarise(n = n(), .groups = "drop")
  n_text <- paste0(n_counts$group, " n = ", format(n_counts$n, big.mark = ","), collapse = " | ")
  
  # Breaks + threshold
  all_vals <- df$e1_e2_obs_scale
  base_breaks <- if (length(all_vals) > 0) scales::log_breaks()(range(all_vals)) else NULL
  thr <- suppressWarnings(as.numeric(EE_contact_threshold))
  add_thr <- if (is.finite(thr) && thr > 0) thr else NULL
  
  ggplot(df, aes(x = e1_e2_obs_scale)) +
    geom_histogram(color = "black", fill = "lightblue", bins = 100) +
    geom_vline(xintercept = thr, color = "red", linetype = "dashed") +
    scale_x_log10(
      breaks = sort(unique(c(base_breaks, add_thr))),
      labels = scales::label_number()
    ) +
    theme_classic() +
    facet_wrap(~ group, nrow = 1, scales = "free_y") +
    labs(
      title = "Histogram of E1–E2 contact values",
      subtitle = paste0(
        if (!is.null(cell_type)) paste0(cell_type, "; ") else "",
        n_text,
        "; May include duplicates if both Enhancers target the same gene"
      ),
      x = "Contact (Hi-C, obs/exp)",
      y = "Count"
    )
}


# All vs Synergistic: P300 boxplots
plot_p300_statistics <- function(all_df,
                                 syn_df = NULL,
                                 p300_threshold,
                                 cell_type = NULL,
                                 box_width = 0.4) {
  # Build long-form data for E1/E2 boxplots
  to_long <- function(df, group_name) {
    df %>%
      transmute(
        group = group_name,
        E1 = normalizedEP300_enh.Feature_element1,
        E2 = normalizedEP300_enh.Feature_element2
      ) %>%
      pivot_longer(cols = c(E1, E2), names_to = "Element", values_to = "P300")
  }
  
  df_all <- to_long(all_df, "All")
  if (!is.null(syn_df)) {
    df_syn <- to_long(syn_df, "Synergistic")
    df_long <- bind_rows(df_all, df_syn)
  } else {
    df_long <- df_all
  }
  
  # n per group (rows, not doubled by E1/E2)
  n_counts <- df_long %>%
    group_by(group) %>%
    summarise(n = n() / 2, .groups = "drop")
  n_text <- paste0(n_counts$group, " n = ", format(as.integer(n_counts$n), big.mark = ","), collapse = " | ")
  
  # For plotting: remove zeros
  df_pos <- df_long %>% filter(is.finite(P300), P300 > 0)
  
  # Breaks + threshold
  yvals_pos <- df_pos$P300
  base_breaks <- if (length(yvals_pos) > 0) scales::log_breaks()(range(yvals_pos)) else NULL
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
    facet_wrap(~ group, nrow = 1) +
    labs(
      title = "P300 distribution by element",
      subtitle = paste0(
        if (!is.null(cell_type)) paste0(cell_type, "; ") else "",
        n_text,
        "; May include duplicates"
      ),
      x = "Element",
      y = "Normalized P300"
    )
}