# National-scale relationships between EFA trends and fire regime -----------
# Calculates Spearman correlations between temporal trends in EFA spatial
# heterogeneity and fire-regime variables across mainland Portugal.
# Benjamini-Hochberg correction is applied across the full set of tests.

source("R/00_config.R")

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

master <- read.csv(
  file.path(OUTPUT_DIR, "Final_Master_Dataset.csv"),
  check.names = FALSE
)

calculate_correlation <- function(data, efa_var, fire_var) {
  valid <- complete.cases(data[[efa_var]], data[[fire_var]])
  n <- sum(valid)

  if (n < 4 ||
      sd(data[[efa_var]][valid]) == 0 ||
      sd(data[[fire_var]][valid]) == 0) {
    return(tibble(
      EFA = efa_var,
      Fire = fire_var,
      N = n,
      Rho = NA_real_,
      P_val = NA_real_
    ))
  }

  test <- cor.test(
    data[[efa_var]][valid],
    data[[fire_var]][valid],
    method = "spearman",
    exact = FALSE
  )

  tibble(
    EFA = efa_var,
    Fire = fire_var,
    N = n,
    Rho = as.numeric(test$estimate),
    P_val = test$p.value
  )
}

results <- crossing(
  EFA = EFA_TREND_VARS,
  Fire = FIRE_VARS
) |>
  pmap_dfr(~ calculate_correlation(master, ..1, ..2)) |>
  mutate(
    p_adjusted = p.adjust(P_val, method = "BH"),
    EFA = recode(EFA, !!!setNames(EFA_LABELS, paste0(names(EFA_LABELS), "_sd_sen_slope"))),
    Fire = recode(Fire, !!!FIRE_LABELS),
    EFA = factor(EFA, levels = unname(EFA_LABELS)),
    Fire = factor(Fire, levels = unname(FIRE_LABELS)),
    Label = if_else(
      is.na(Rho), "",
      if_else(p_adjusted < 0.05, paste0(sprintf("%.2f", Rho), "*"), sprintf("%.2f", Rho))
    )
  )

write.csv(
  results,
  file.path(OUTPUT_DIR, "Table_EFA_Trends_Fire_Variables_National.csv"),
  row.names = FALSE
)

p <- ggplot(results, aes(x = Fire, y = EFA, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = Label), size = 3.5) +
  scale_fill_gradient2(
    low = "#d8b365", mid = "white", high = "#5ab4ac",
    midpoint = 0, limits = c(-0.4, 0.4), na.value = "grey90",
    name = "Spearman ρ"
  ) +
  labs(
    x = "Fire-regime variable",
    y = "EFA indicator",
    title = "Relationships between EFA trends and fire-regime characteristics",
    subtitle = "* indicates Benjamini-Hochberg adjusted p < 0.05"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p)

ggsave(
  file.path(OUTPUT_DIR, "EFA_Trends_Fire_Variables_National.png"),
  p,
  width = 11,
  height = 6,
  dpi = 300
)
