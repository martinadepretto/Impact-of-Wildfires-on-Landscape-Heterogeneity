# Fire-regime cluster-level relationships -----------------------------------
# Calculates Spearman correlations between significant EFA trends and
# fire-regime variables within each fire-regime cluster.
# Cluster 4 is excluded because of its very small sample size.
# Benjamini-Hochberg correction is applied separately within each cluster.

source("R/00_config.R")

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

master <- read.csv(
  file.path(OUTPUT_DIR, "Final_Master_Dataset.csv"),
  check.names = FALSE
)

cluster_levels <- master |>
  filter(!is.na(cluster), cluster != 4) |>
  distinct(cluster) |>
  arrange(cluster) |>
  pull(cluster)

calculate_cluster_correlation <- function(data, cluster_id, efa_var, p_var, fire_var) {
  dat <- data |>
    filter(
      cluster == cluster_id,
      !is.na(.data[[efa_var]]),
      !is.na(.data[[p_var]]),
      .data[[p_var]] < 0.05
    )

  valid <- complete.cases(dat[[efa_var]], dat[[fire_var]])
  n <- sum(valid)

  if (n < 4 ||
      sd(dat[[efa_var]][valid]) == 0 ||
      sd(dat[[fire_var]][valid]) == 0) {
    return(tibble(
      Cluster = cluster_id,
      EFA = efa_var,
      Fire = fire_var,
      N = n,
      Rho = NA_real_,
      P_val = NA_real_
    ))
  }

  test <- cor.test(
    dat[[efa_var]][valid],
    dat[[fire_var]][valid],
    method = "spearman",
    exact = FALSE
  )

  tibble(
    Cluster = cluster_id,
    EFA = efa_var,
    Fire = fire_var,
    N = n,
    Rho = as.numeric(test$estimate),
    P_val = test$p.value
  )
}

results <- map_dfr(cluster_levels, function(cl) {
  map_dfr(seq_along(EFA_TREND_VARS), function(i) {
    map_dfr(FIRE_VARS, function(fire_var) {
      calculate_cluster_correlation(
        master,
        cluster_id = cl,
        efa_var = EFA_TREND_VARS[i],
        p_var = EFA_PVARS[i],
        fire_var = fire_var
      )
    })
  })
}) |>
  group_by(Cluster) |>
  mutate(p_adjusted = p.adjust(P_val, method = "BH")) |>
  ungroup() |>
  mutate(
    EFA = recode(EFA, !!!setNames(EFA_LABELS, paste0(names(EFA_LABELS), "_sd_sen_slope"))),
    Fire = recode(Fire, !!!FIRE_LABELS),
    EFA = factor(EFA, levels = unname(EFA_LABELS)),
    Fire = factor(Fire, levels = unname(FIRE_LABELS)),
    Cluster = factor(paste0("Cluster ", Cluster), levels = paste0("Cluster ", cluster_levels)),
    Label = if_else(
      is.na(Rho), "",
      if_else(p_adjusted < 0.05, paste0(sprintf("%.2f", Rho), "*"), sprintf("%.2f", Rho))
    )
  )

write.csv(
  results,
  file.path(OUTPUT_DIR, "Table_EFA_Trends_Fire_Variables_by_Cluster.csv"),
  row.names = FALSE
)

p <- ggplot(results, aes(x = Fire, y = EFA, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = Label), size = 3.5) +
  facet_wrap(~ Cluster, ncol = 3) +
  scale_fill_gradient2(
    low = "#d8b365", mid = "white", high = "#5ab4ac",
    midpoint = 0, limits = c(-0.4, 0.4), na.value = "grey90",
    name = "Spearman ρ"
  ) +
  labs(
    x = "Fire-regime variable",
    y = "EFA indicator",
    title = "Relationships between EFA trends and fire-regime characteristics by cluster",
    subtitle = "Analysis restricted to significant EFA trends; * indicates BH-adjusted p < 0.05"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p)

ggsave(
  file.path(OUTPUT_DIR, "EFA_Trends_Fire_Variables_by_Cluster.png"),
  p,
  width = 14,
  height = 9,
  dpi = 300
)
