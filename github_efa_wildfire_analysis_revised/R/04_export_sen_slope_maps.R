# Export Sen's slope maps for QGIS -------------------------------------------
# Creates a GeoPackage containing EFA trend estimates, significance flags,
# and fire-regime cluster membership for national-scale cartography.

source("R/00_config.R")

library(terra)
library(dplyr)

master <- read.csv(
  file.path(OUTPUT_DIR, "Final_Master_Dataset.csv"),
  check.names = FALSE
)

grid <- vect(GRID_PATH)

map_data <- master |>
  select(id, cluster, all_of(EFA_TREND_VARS), all_of(EFA_PVARS)) |>
  mutate(
    NDVI_significant = NDVI_sd_sen_p_value < 0.05,
    LST_significant = LST_sd_sen_p_value < 0.05,
    TCTwet_significant = TCTwet_sd_sen_p_value < 0.05,
    TCTbright_significant = TCTbright_sd_sen_p_value < 0.05
  )

map_data$id <- as.character(map_data$id)
grid$id <- as.character(grid$id)

map_spatial <- merge(grid, map_data, by = "id")

writeVector(
  map_spatial,
  file.path(OUTPUT_DIR, "EFA_SenSlope_Maps.gpkg"),
  overwrite = TRUE
)

message("QGIS-ready trend layer written to: ", OUTPUT_DIR)
