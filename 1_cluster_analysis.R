# ============================================================
# Fire Regime Cluster Analysis - Mainland Portugal
# ============================================================
#
# Purpose:
#   1. Load and prepare the input dataset
#   2. Run k-means clustering on fire and land-cover variables
#   3. Evaluate candidate numbers of clusters using Elbow and
#      Silhouette methods
#   4. Characterize the selected clusters
#   5. Export tables, rankings, figures, and a QGIS-ready dataset
#
# Recommended project structure:
#   project/
#   ├── data/
#   │   └── 7_variables_EU_completed_table.csv
#   ├── outputs/
#   └── scripts/
#       └── cluster_analysis.R
#
# Required packages:
#   tidyverse, factoextra, cluster, openxlsx, GGally
#
# ============================================================


# ------------------------------------------------------------
# 0. PACKAGE SETUP
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "factoextra",
  "cluster",
  "openxlsx",
  "GGally"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(factoextra)
  library(cluster)
  library(openxlsx)
  library(GGally)
})


# ------------------------------------------------------------
# 1. PROJECT PATHS AND SETTINGS
# ------------------------------------------------------------

# Change these paths only if your repository uses a different structure.
input_dir <- "data"
output_dir <- "outputs"

input_file <- "7_variables_EU_completed_table.csv"

# Selected number of clusters
k_final <- 6

# Random seed for reproducibility
seed <- 123

# Number of random starts for k-means
nstart <- 25

# Maximum number of iterations for k-means
iter_max <- 100

# Create output directory if it does not exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------

input_path <- file.path(input_dir, input_file)

if (!file.exists(input_path)) {
  stop(
    "Input file not found: ", input_path,
    "\nPlace the CSV file in the 'data' directory or update 'input_dir'."
  )
}

master_df <- read_csv(
  input_path,
  show_col_types = FALSE
)

message("Input data loaded successfully.")
message("Rows: ", nrow(master_df))
message("Columns: ", ncol(master_df))


# ------------------------------------------------------------
# 3. DEFINE VARIABLES
# ------------------------------------------------------------

clustering_variables <- c(
  "FF_mean",
  "TBA_sum",
  "TSF_median",
  "Forest_per",
  "Grass_perc",
  "Shrub_perc"
)

severity_columns <- names(master_df)[
  grepl("^sev", names(master_df), ignore.case = TRUE)
]

if (length(severity_columns) == 0) {
  stop(
    "No severity columns were found.",
    "\nExpected columns starting with 'sev'."
  )
}

required_columns <- c(clustering_variables, severity_columns)

missing_columns <- setdiff(required_columns, names(master_df))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing:\n",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 4. HANDLE MISSING VALUES
# ------------------------------------------------------------

# IMPORTANT:
# NA values are replaced with zero only in the variables used
# for clustering and severity calculations.
#
# This reproduces the original analysis assumption that missing
# values represent no recorded fire / zero contribution.

analysis_columns <- unique(c(clustering_variables, severity_columns))

master_df <- master_df %>%
  mutate(
    across(
      all_of(analysis_columns),
      ~ replace_na(.x, 0)
    )
  )


# ------------------------------------------------------------
# 5. PREPARE DATA FOR CLUSTERING
# ------------------------------------------------------------

data_for_clustering <- master_df %>%
  select(all_of(clustering_variables), all_of(severity_columns))

# Check that all clustering variables are numeric
non_numeric_columns <- names(data_for_clustering)[
  !vapply(data_for_clustering, is.numeric, logical(1))
]

if (length(non_numeric_columns) > 0) {
  stop(
    "The following clustering variables are not numeric:\n",
    paste(non_numeric_columns, collapse = ", ")
  )
}

# Check for missing or infinite values after cleaning
if (any(!is.finite(as.matrix(data_for_clustering)))) {
  stop("Non-finite values remain in the clustering dataset.")
}


# ------------------------------------------------------------
# 6. STANDARDIZE VARIABLES
# ------------------------------------------------------------

# Standardization is required because variables are measured
# on different scales (e.g. FF vs. severity values).

data_scaled <- scale(data_for_clustering)

# Check for zero-variance variables, which cannot be standardized.
zero_variance <- apply(data_for_clustering, 2, sd) == 0

if (any(zero_variance)) {
  stop(
    "Zero-variance variables detected:\n",
    paste(names(zero_variance)[zero_variance], collapse = ", "),
    "\nRemove or reconsider these variables before clustering."
  )
}


# ------------------------------------------------------------
# 7. CLUSTER NUMBER DIAGNOSTICS
# ------------------------------------------------------------

set.seed(seed)

elbow_plot <- fviz_nbclust(
  data_scaled,
  FUNcluster = kmeans,
  method = "wss",
  k.max = 10
) +
  labs(
    title = "Elbow Method",
    subtitle = "Within-cluster sum of squares"
  )

set.seed(seed)

silhouette_plot <- fviz_nbclust(
  data_scaled,
  FUNcluster = kmeans,
  method = "silhouette",
  k.max = 10
) +
  labs(
    title = "Silhouette Method",
    subtitle = "Average silhouette width"
  )

ggsave(
  filename = file.path(output_dir, "cluster_diagnostic_elbow.png"),
  plot = elbow_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "cluster_diagnostic_silhouette.png"),
  plot = silhouette_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 8. FINAL K-MEANS CLUSTERING
# ------------------------------------------------------------

set.seed(seed)

final_clusters <- kmeans(
  data_scaled,
  centers = k_final,
  nstart = nstart,
  iter.max = iter_max
)

message(
  "K-means clustering completed with k = ",
  k_final,
  "."
)

# Add numeric cluster ID to the original dataset
master_df <- master_df %>%
  mutate(
    cluster = final_clusters$cluster
  )


# ------------------------------------------------------------
# 9. CLUSTER VISUALIZATION IN PCA SPACE
# ------------------------------------------------------------

cluster_plot <- fviz_cluster(
  final_clusters,
  data = data_scaled,
  geom = "point",
  ellipse.type = "convex"
) +
  labs(
    title = paste0("Cluster Visualization for k = ", k_final)
  ) +
  theme_minimal()

ggsave(
  filename = file.path(output_dir, "cluster_pca_visualization.png"),
  plot = cluster_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 10. CALCULATE MAXIMUM SEVERITY
# ------------------------------------------------------------

# Maximum recorded severity for each spatial unit across
# all available severity columns.

master_df <- master_df %>%
  mutate(
    max_sev = apply(
      select(., all_of(severity_columns)),
      1,
      max
    )
  )


# ------------------------------------------------------------
# 11. CLUSTER SUMMARY TABLE
# ------------------------------------------------------------

table_cluster <- master_df %>%
  group_by(cluster) %>%
  summarise(
    Cells_number = n(),
    FF = round(mean(FF_mean, na.rm = TRUE), 2),
    TBA_ha = round(mean(TBA_sum, na.rm = TRUE), 1),
    TSF_yr = round(median(TSF_median, na.rm = TRUE), 1),
    Sev_Mean_Max = round(mean(max_sev, na.rm = TRUE), 1),
    Forest_perc = round(mean(Forest_per, na.rm = TRUE), 1),
    Shrub_perc = round(mean(Shrub_perc, na.rm = TRUE), 1),
    Grass_perc = round(mean(Grass_perc, na.rm = TRUE), 1),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 12. DEFINE CLUSTER LABELS
# ------------------------------------------------------------

labels_regimes <- c(
  "1" = "C1 - High-Recurrence, High-Severity Regime",
  "2" = "C2 - Dynamic Forest-Dominated Regime",
  "3" = "C3 - High-Severity, Low-Recurrence Regime",
  "4" = "C4 - Shrubland-Dominated Hotspots",
  "5" = "C5 - Stable, Low-Fire Regime",
  "6" = "C6 - Moderate-Intensity Transition Regime"
)

table_with_labels <- table_cluster %>%
  mutate(
    Label_Regime = labels_regimes[as.character(cluster)]
  ) %>%
  select(
    cluster,
    Label_Regime,
    everything()
  )


# ------------------------------------------------------------
# 13. EXPORT CLUSTER STATISTICS
# ------------------------------------------------------------

write_csv(
  table_with_labels,
  file.path(output_dir, "Statistics_Clusters.csv")
)


# ------------------------------------------------------------
# 14. RANK CLUSTERS BY VARIABLE
# ------------------------------------------------------------

variables_to_rank <- c(
  "FF",
  "TBA_ha",
  "TSF_yr",
  "Sev_Mean_Max",
  "Forest_perc",
  "Shrub_perc",
  "Grass_perc"
)

ranking_tables <- list()

for (variable in variables_to_rank) {

  ranking_tables[[variable]] <- table_with_labels %>%
    select(
      cluster,
      Label_Regime,
      all_of(variable)
    ) %>%
    arrange(desc(.data[[variable]]))
}


# ------------------------------------------------------------
# 15. EXPORT RANKINGS TO EXCEL
# ------------------------------------------------------------

ranking_workbook <- createWorkbook()

for (variable in variables_to_rank) {

  addWorksheet(
    ranking_workbook,
    sheetName = variable
  )

  writeData(
    ranking_workbook,
    sheet = variable,
    x = ranking_tables[[variable]]
  )
}

saveWorkbook(
  ranking_workbook,
  file.path(output_dir, "Cluster_Rankings_By_Variable.xlsx"),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 16. PREPARE DATA FOR COMPARATIVE FIGURE
# ------------------------------------------------------------

variable_labels <- c(
  "FF" = "FF",
  "TSF_yr" = "TSF (yr)",
  "TBA_ha" = "TBA (ha)",
  "Sev_Mean_Max" = "Maximum Severity",
  "Forest_perc" = "Forest (%)",
  "Shrub_perc" = "Shrub (%)",
  "Grass_perc" = "Grass (%)"
)

table_long <- table_with_labels %>%
  pivot_longer(
    cols = all_of(variables_to_rank),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(
    Variable = recode(
      Variable,
      !!!variable_labels
    )
  )


# ------------------------------------------------------------
# 17. COMPARATIVE CLUSTER FIGURE
# ------------------------------------------------------------

# The colours are kept consistent across all figures so that
# each cluster has the same visual identity throughout the thesis.

palette_regimes <- c(
  "C1 - High-Recurrence, High-Severity Regime" = "#4a148c",
  "C2 - Dynamic Forest-Dominated Regime" = "#8d6e63",
  "C3 - High-Severity, Low-Recurrence Regime" = "#ff5a00",
  "C4 - Shrubland-Dominated Hotspots" = "#b71c1c",
  "C5 - Stable, Low-Fire Regime" = "#81c784",
  "C6 - Moderate-Intensity Transition Regime" = "#ffb300"
)

comparative_plot <- ggplot(
  table_long,
  aes(
    x = factor(cluster),
    y = Value,
    fill = Label_Regime
  )
) +
  geom_col(
    alpha = 0.85,
    colour = "black",
    linewidth = 0.2
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free_y",
    ncol = 3,
    axes = "all"
  ) +
  scale_fill_manual(values = palette_regimes) +
  labs(
    title = "Comparative Characterization of the Six Identified Fire Regimes",
    subtitle = "Comparison of fire-related and land-cover variables across k = 6 clusters",
    x = "Cluster ID",
    y = "Variable value",
    fill = "Fire Regime"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "comparative_fire_regimes.png"),
  plot = comparative_plot,
  width = 12,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 18. PARALLEL COORDINATES PLOT
# ------------------------------------------------------------

pcp_data <- table_with_labels %>%
  select(
    Label_Regime,
    all_of(variables_to_rank)
  )

parallel_plot <- ggparcoord(
  pcp_data,
  columns = 2:8,
  groupColumn = 1,
  scale = "uniminmax",
  showPoints = TRUE,
  alphaLines = 0.6
) +
  scale_color_manual(values = palette_regimes) +
  labs(
    title = "Fire Regime Profiles",
    x = "Fire and Ecological Variables",
    y = "Normalized Value (0 to 1)",
    color = "Fire Regime"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 20,
      vjust = 0.5,
      face = "bold"
    ),
    legend.position = "bottom",
    panel.grid.major.y = element_line(
      colour = "grey90"
    )
  )

ggsave(
  filename = file.path(output_dir, "fire_regime_parallel_coordinates.png"),
  plot = parallel_plot,
  width = 12,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 19. EXPORT QGIS-READY DATASET
# ------------------------------------------------------------

master_df <- master_df %>%
  mutate(
    cluster = as.integer(cluster),
    Label_Regime = labels_regimes[as.character(cluster)]
  )

write_csv(
  master_df,
  file.path(output_dir, "Results_CA_final.csv")
)


# ------------------------------------------------------------
# 20. EXPORT FINAL R OBJECTS (OPTIONAL)
# ------------------------------------------------------------

saveRDS(
  final_clusters,
  file.path(output_dir, "kmeans_k6_model.rds")
)

saveRDS(
  table_with_labels,
  file.path(output_dir, "cluster_statistics.rds")
)


# ------------------------------------------------------------
# 21. CONSOLE SUMMARY
# ------------------------------------------------------------

message("\n==============================================")
message("FIRE REGIME CLUSTER ANALYSIS COMPLETED")
message("==============================================")
message("Selected number of clusters: k = ", k_final)
message("Number of observations: ", nrow(master_df))
message("Number of severity variables: ", length(severity_columns))
message("\nCluster sizes:")

print(
  master_df %>%
    count(cluster, name = "Cells_number")
)

message("\nOutput files saved in: ", normalizePath(output_dir))
message("==============================================")
