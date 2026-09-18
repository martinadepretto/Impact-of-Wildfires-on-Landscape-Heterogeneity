# Project configuration -------------------------------------------------------
# Define project-relative paths and shared analysis settings.
# Raw data should remain outside version control unless licensing permits sharing.

PROJECT_ROOT <- if (requireNamespace("here", quietly = TRUE)) {
  here::here()
} else {
  getwd()
}

GRID_PATH <- file.path(PROJECT_ROOT, "data", "spatial", "Reticolo.gpkg")
FIRE_DATA_PATH <- file.path(PROJECT_ROOT, "data", "fire", "Results_CA.gpkg")
EFA_ROOT <- file.path(PROJECT_ROOT, "data", "EFA")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "outputs")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

EFA_RASTER_DIRS <- list(
  NDVI = file.path(EFA_ROOT, "NDVI", "NDVI_clip"),
  LST = file.path(EFA_ROOT, "LST_Day", "LST_Day_clip"),
  TCTwet = file.path(EFA_ROOT, "TCT_Wetness", "TCT_Wetness_clip"),
  TCTbright = file.path(EFA_ROOT, "TCT_Brightness", "TCT_Brightness_clip")
)

EFA_TREND_VARS <- c(
  "NDVI_sd_sen_slope",
  "LST_sd_sen_slope",
  "TCTwet_sd_sen_slope",
  "TCTbright_sd_sen_slope"
)

EFA_PVARS <- c(
  "NDVI_sd_sen_p_value",
  "LST_sd_sen_p_value",
  "TCTwet_sd_sen_p_value",
  "TCTbright_sd_sen_p_value"
)

FIRE_VARS <- c(
  "TBA_sum",
  "FF_mean",
  "TSF_median",
  "max_severity",
  "mean_severity"
)

EFA_LABELS <- c(
  NDVI = "NDVI",
  LST = "LST",
  TCTwet = "TCT Wetness",
  TCTbright = "TCT Brightness"
)

FIRE_LABELS <- c(
  TBA_sum = "Total burned area",
  FF_mean = "Fire frequency",
  TSF_median = "Time since fire",
  max_severity = "Maximum severity",
  mean_severity = "Mean severity"
)
