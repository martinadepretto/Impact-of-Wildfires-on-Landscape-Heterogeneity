# ==============================================================================
# Post-fire spatial heterogeneity analysis
# ==============================================================================
#
# Purpose
# -------
# Quantify changes in within-site spatial dispersion of four Ecosystem
# Functional Attributes (EFAs) following the 2005 wildfire and test whether
# those changes are associated with wildfire extent.
#
# Main response variables
#   - ΔMAD: change in within-site spatial dispersion
#   - %ΔMAD: relative change in spatial dispersion
#
# Primary relationship
#   - total burned area vs ΔMAD
#
# Secondary analyses
#   - total burned area vs %ΔMAD
#   - ΔMean vs ΔMAD
#   - exploratory Structural Sensitivity Ratio (SSR)
#
# Required objects from the preceding local EFA-processing step
#   heterogeneity_summary_df
#   selected_focal_sites_with_area
#   mean_mad_trajectory_df       (for the trajectory figure)
#   df_veg_analysis              (for the land-cover figure)
#   sintesi_correlazione_veg     (optional table export)
#   sintesi_dominanza_estesa     (optional table export)
#
# The focal sites were deliberately selected to represent contrasting NDVI
# trajectories within strata. Results therefore describe the selected sites
# and should not be interpreted as unbiased estimates for the whole region.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(ggrepel)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

output_dir <- "Local_Analysis_Output"
table_dir <- file.path(output_dir, "Tables")
figure_dir <- file.path(output_dir, "Figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

efa_order <- c("NDVI", "LST", "TCTB", "TCTW")
fire_year <- 2005
pre_period <- "2000-2004"
post_period <- "2006-2024"
epsilon_mean <- 1e-6

required_objects <- c(
  "heterogeneity_summary_df",
  "selected_focal_sites_with_area"
)
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "The following objects must be created by the preceding EFA-processing step: ",
    paste(missing_objects, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 2. Assemble the focal-site analysis dataset
# ------------------------------------------------------------------------------

area_lookup <- selected_focal_sites_with_area %>%
  st_drop_geometry() %>%
  transmute(
    site_id = trimws(as.character(id)),
    size_class,
    total_fire_area_ha
  ) %>%
  distinct(site_id, .keep_all = TRUE)

analysis_df <- heterogeneity_summary_df %>%
  mutate(
    site_id = trimws(as.character(site_id)),
    metric = toupper(metric)
  ) %>%
  inner_join(area_lookup, by = "site_id") %>%
  mutate(
    delta_mad = post_mad_mean - pre_mad_baseline_mean,
    pct_delta_mad = if_else(
      is.finite(pre_mad_baseline_mean) & pre_mad_baseline_mean != 0,
      100 * delta_mad / abs(pre_mad_baseline_mean),
      NA_real_
    ),
    delta_mean = post_mean_mean - pre_mean_baseline_mean,
    pct_delta_mean = if_else(
      is.finite(pre_mean_baseline_mean) & pre_mean_baseline_mean != 0,
      100 * delta_mean / abs(pre_mean_baseline_mean),
      NA_real_
    ),
    heterogeneity_response = case_when(
      is.na(delta_mad) ~ NA_character_,
      delta_mad < 0 ~ "Decreased spatial dispersion",
      delta_mad > 0 ~ "Increased spatial dispersion",
      TRUE ~ "No change"
    ),
    heterogeneity_trajectory = case_when(
      is.na(pct_delta_mad) ~ NA_character_,
      pct_delta_mad < -15 ~ "Strong homogenization",
      pct_delta_mad < -5 ~ "Moderate homogenization",
      pct_delta_mad <= 5 ~ "Stable pattern",
      pct_delta_mad <= 15 ~ "Moderate heterogenization",
      TRUE ~ "Strong heterogenization"
    )
  )

# ------------------------------------------------------------------------------
# 3. Summary statistics by EFA
# ------------------------------------------------------------------------------

heterogeneity_summary <- analysis_df %>%
  group_by(metric) %>%
  summarise(
    n_sites = n_distinct(site_id),
    mean_pre_MAD = mean(pre_mad_baseline_mean, na.rm = TRUE),
    sd_pre_MAD = sd(pre_mad_baseline_mean, na.rm = TRUE),
    mean_post_MAD = mean(post_mad_mean, na.rm = TRUE),
    sd_post_MAD = sd(post_mad_mean, na.rm = TRUE),
    mean_delta_MAD = mean(delta_mad, na.rm = TRUE),
    sd_delta_MAD = sd(delta_mad, na.rm = TRUE),
    median_delta_MAD = median(delta_mad, na.rm = TRUE),
    mean_pct_delta_MAD = mean(pct_delta_mad, na.rm = TRUE),
    sd_pct_delta_MAD = sd(pct_delta_mad, na.rm = TRUE),
    median_pct_delta_MAD = median(pct_delta_mad, na.rm = TRUE),
    n_decreased = sum(delta_mad < 0, na.rm = TRUE),
    n_increased = sum(delta_mad > 0, na.rm = TRUE),
    pct_decreased = 100 * mean(delta_mad < 0, na.rm = TRUE),
    pct_increased = 100 * mean(delta_mad > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(metric, levels = efa_order))

trajectory_summary <- analysis_df %>%
  count(metric, heterogeneity_trajectory, name = "n") %>%
  group_by(metric) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup() %>%
  arrange(factor(metric, levels = efa_order), heterogeneity_trajectory)

# ------------------------------------------------------------------------------
# 4. Statistical helper
# ------------------------------------------------------------------------------

spearman_by_metric <- function(data, x, y) {
  data %>%
    group_by(metric) %>%
    group_modify(~ {
      complete <- .x %>%
        filter(is.finite(.data[[x]]), is.finite(.data[[y]]))

      if (nrow(complete) < 3) {
        return(tibble(n = nrow(complete), rho = NA_real_, p_value = NA_real_))
      }

      test <- cor.test(
        complete[[x]],
        complete[[y]],
        method = "spearman",
        exact = FALSE
      )

      tibble(
        n = nrow(complete),
        rho = unname(test$estimate),
        p_value = test$p.value
      )
    }) %>%
    ungroup() %>%
    arrange(factor(metric, levels = efa_order))
}

# ------------------------------------------------------------------------------
# 5. Primary analysis: wildfire extent vs change in spatial dispersion
# ------------------------------------------------------------------------------

fire_area_delta_mad <- spearman_by_metric(
  analysis_df,
  "total_fire_area_ha",
  "delta_mad"
)

fire_area_pct_delta_mad <- spearman_by_metric(
  analysis_df,
  "total_fire_area_ha",
  "pct_delta_mad"
)

fire_heterogeneity_relationship <- fire_area_delta_mad %>%
  rename(
    rho_delta_MAD = rho,
    p_delta_MAD = p_value
  ) %>%
  left_join(
    fire_area_pct_delta_mad %>%
      select(metric, rho, p_value) %>%
      rename(
        rho_pct_delta_MAD = rho,
        p_pct_delta_MAD = p_value
      ),
    by = "metric"
  )

# ------------------------------------------------------------------------------
# 6. Secondary analysis: functional-state change vs spatial dispersion
# ------------------------------------------------------------------------------

mean_mad_relationship <- spearman_by_metric(
  analysis_df,
  "delta_mean",
  "delta_mad"
) %>%
  rename(
    rho_delta_mean_delta_mad = rho,
    p_delta_mean_delta_mad = p_value
  )

# ------------------------------------------------------------------------------
# 7. Exploratory Structural Sensitivity Ratio
# ------------------------------------------------------------------------------
# SSR is a descriptive ratio, not a causal sensitivity measure. Values can become
# unstable when the denominator approaches zero.

analysis_df <- analysis_df %>%
  mutate(
    SSR = if_else(
      is.finite(delta_mean) & is.finite(delta_mad) & abs(delta_mean) > epsilon_mean,
      abs(delta_mad) / abs(delta_mean),
      NA_real_
    ),
    relative_SSR = if_else(
      is.finite(pct_delta_mean) & is.finite(pct_delta_mad) &
        abs(pct_delta_mean) > epsilon_mean,
      abs(pct_delta_mad) / abs(pct_delta_mean),
      NA_real_
    )
  )

SSR_summary <- analysis_df %>%
  group_by(metric) %>%
  summarise(
    n_valid_SSR = sum(is.finite(SSR)),
    mean_SSR = mean(SSR, na.rm = TRUE),
    median_SSR = median(SSR, na.rm = TRUE),
    mean_relative_SSR = mean(relative_SSR, na.rm = TRUE),
    median_relative_SSR = median(relative_SSR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(metric, levels = efa_order))

# ------------------------------------------------------------------------------
# 8. Plotting helpers
# ------------------------------------------------------------------------------

plot_scatter_by_efa <- function(data, x, y, x_label, y_label, title) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2.8, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~metric, scales = "free_y") +
    labs(title = title, x = x_label, y = y_label) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

p_fire_delta_mad <- plot_scatter_by_efa(
  analysis_df,
  "total_fire_area_ha",
  "delta_mad",
  "Total burned area (ha)",
  expression(Delta * " MAD"),
  "Wildfire extent and post-fire change in spatial dispersion"
)

p_fire_pct_delta_mad <- plot_scatter_by_efa(
  analysis_df,
  "total_fire_area_ha",
  "pct_delta_mad",
  "Total burned area (ha)",
  "% change in MAD",
  "Wildfire extent and relative change in spatial dispersion"
)

p_delta_mean_mad <- analysis_df %>%
  ggplot(aes(delta_mean, delta_mad)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~metric, scales = "free") +
  labs(
    title = "Coupling between functional-state change and spatial dispersion",
    x = expression(Delta * " Mean"),
    y = expression(Delta * " MAD")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 9. Mean–MAD trajectory figure
# ------------------------------------------------------------------------------

if (exists("mean_mad_trajectory_df")) {
  trajectory_df <- mean_mad_trajectory_df %>%
    mutate(metric = toupper(metric))

  plot_trajectory <- function(metric_name) {
    ggplot(
      filter(trajectory_df, metric == metric_name),
      aes(delta_mean, delta_mad, color = trajectory_concordance, shape = cluster)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_point(size = 3, alpha = 0.85) +
      ggrepel::geom_text_repel(
        aes(label = site_id),
        size = 2.8,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      labs(
        title = metric_name,
        x = expression(Delta * " Mean"),
        y = expression(Delta * " MAD"),
        color = "Trajectory",
        shape = "Cluster"
      ) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
  }

  trajectory_figure <-
    (plot_trajectory("NDVI") | plot_trajectory("LST")) /
    (plot_trajectory("TCTB") | plot_trajectory("TCTW")) +
    plot_annotation(title = "Concordance between changes in EFA mean and spatial dispersion")

  print(trajectory_figure)
  ggsave(
    file.path(figure_dir, "Mean_MAD_Trajectory_Concordance.png"),
    trajectory_figure,
    width = 12,
    height = 9,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 10. Land-cover context figure
# ------------------------------------------------------------------------------
# Land-cover variables are contemporary COS2023 descriptors and are therefore
# interpreted as contextual associations rather than pre-fire explanatory factors.

if (exists("df_veg_analysis")) {
  landcover_plot <- function(metric_name, period = c("pre", "post")) {
    period <- match.arg(period)
    mad_var <- if (period == "pre") "pre_mad_baseline_mean" else "post_mad_mean"
    period_label <- if (period == "pre") "Pre-fire (2000-2004)" else "Post-fire (2006-2024)"

    dat <- df_veg_analysis %>%
      filter(toupper(metric) == metric_name) %>%
      transmute(
        site_id,
        class_richness = n_patches,
        MAD = .data[[mad_var]]
      ) %>%
      filter(is.finite(class_richness), is.finite(MAD))

    ggplot(dat, aes(class_richness, MAD)) +
      geom_point(size = 2.8, alpha = 0.8) +
      geom_smooth(method = "lm", se = FALSE) +
      geom_text_repel(aes(label = site_id), size = 2.7, show.legend = FALSE) +
      labs(
        title = metric_name,
        subtitle = period_label,
        x = "Land-cover class richness",
        y = "Spatial dispersion (MAD)"
      ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
  }

  landcover_figure <-
    (landcover_plot("NDVI", "pre") | landcover_plot("NDVI", "post")) /
    (landcover_plot("LST", "pre") | landcover_plot("LST", "post")) +
    plot_annotation(title = "Land-cover class richness and spatial heterogeneity")

  print(landcover_figure)
  ggsave(
    file.path(figure_dir, "Landcover_Class_Richness_MAD.png"),
    landcover_figure,
    width = 11,
    height = 9,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 11. Export analysis tables
# ------------------------------------------------------------------------------

final_site_level <- analysis_df %>%
  select(
    site_id,
    metric,
    size_class,
    total_fire_area_ha,
    pre_mean_baseline_mean,
    post_mean_mean,
    delta_mean,
    pct_delta_mean,
    pre_mad_baseline_mean,
    post_mad_mean,
    delta_mad,
    pct_delta_mad,
    heterogeneity_response,
    heterogeneity_trajectory,
    SSR,
    relative_SSR
  ) %>%
  arrange(factor(metric, levels = efa_order), site_id)

write_csv(final_site_level, file.path(table_dir, "Site_Level_Postfire_Heterogeneity.csv"))
write_csv(heterogeneity_summary, file.path(table_dir, "Heterogeneity_Summary_By_EFA.csv"))
write_csv(trajectory_summary, file.path(table_dir, "Heterogeneity_Trajectory_Summary.csv"))
write_csv(fire_heterogeneity_relationship, file.path(table_dir, "Fire_Area_Heterogeneity_Relationship.csv"))
write_csv(mean_mad_relationship, file.path(table_dir, "Delta_Mean_Delta_MAD_Correlations.csv"))
write_csv(SSR_summary, file.path(table_dir, "Structural_Sensitivity_Ratio_Summary.csv"))

# Optional exports generated by preceding land-cover analyses.
optional_exports <- list(
  sintesi_correlazione_veg = "Landcover_Richness_Shannon_Correlations.csv",
  sintesi_dominanza_estesa = "Landcover_Dominance_Correlations.csv",
  sintesi_range_frammentazione = "Landcover_Richness_vs_EFA_Ranges.csv",
  sintesi_range_dominanza = "Landcover_Dominance_vs_EFA_Ranges.csv"
)

for (object_name in names(optional_exports)) {
  if (exists(object_name)) {
    write_csv(
      get(object_name) %>% arrange(metric),
      file.path(table_dir, optional_exports[[object_name]])
    )
  }
}

# ------------------------------------------------------------------------------
# 12. Export figures
# ------------------------------------------------------------------------------

ggsave(file.path(figure_dir, "Fire_Area_vs_Delta_MAD.png"), p_fire_delta_mad, width = 10, height = 7, dpi = 300)
ggsave(file.path(figure_dir, "Fire_Area_vs_Percent_Delta_MAD.png"), p_fire_pct_delta_mad, width = 10, height = 7, dpi = 300)
ggsave(file.path(figure_dir, "Delta_Mean_vs_Delta_MAD.png"), p_delta_mean_mad, width = 10, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 13. Console summary
# ------------------------------------------------------------------------------

message("Post-fire heterogeneity analysis completed.")
message("Primary fire-area relationships:")
print(fire_heterogeneity_relationship)
message("Mean–MAD relationships:")
print(mean_mad_relationship)
message("SSR summary:")
print(SSR_summary)
