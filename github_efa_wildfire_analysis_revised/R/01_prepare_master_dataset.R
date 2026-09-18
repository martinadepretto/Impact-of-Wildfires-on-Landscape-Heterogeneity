# Prepare national-scale EFA heterogeneity dataset ---------------------------
# For each EFA and year, calculate spatial standard deviation across pixels
# within each 5 x 5 km grid cell, then estimate temporal trends and merge
# the results with fire-regime variables.

source("R/00_config.R")

library(terra)
library(dplyr)
library(tidyr)
purrr::walk(c("trend"), require, character.only = TRUE)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

process_efa <- function(raster_dir, variable) {
  files <- sort(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
  if (length(files) == 0) stop("No TIFF files found in: ", raster_dir)

  rasters <- rast(files)
  if (nlyr(rasters) != 25) {
    warning(variable, ": expected 25 annual rasters (2000-2024), found ", nlyr(rasters), ".")
  }

  names(rasters) <- paste0(variable, "_", 2000:(2000 + nlyr(rasters) - 1))
  grid <- vect(GRID_PATH)
  grid_ids <- grid$id

  sd_wide <- zonal(rasters, grid, fun = "sd", na.rm = TRUE) |>
    as.data.frame()
  sd_wide$id <- grid_ids

  sd_wide |>
    pivot_longer(
      cols = starts_with(paste0(variable, "_")),
      names_to = "year",
      names_pattern = paste0(variable, "_(\\d{4})"),
      values_to = paste0(variable, "_sd")
    ) |>
    mutate(year = as.integer(year)) |>
    select(id, year, everything())
}

calculate_trends <- function(data, variable) {
  value_col <- paste0(variable, "_sd")

  data |>
    group_by(id) |>
    summarise(
      !!paste0(variable, "_sd_sen_slope") := tryCatch(
        as.numeric(trend::sens.slope(.data[[value_col]])$estimates),
        error = function(e) NA_real_
      ),
      !!paste0(variable, "_sd_sen_p_value") := tryCatch(
        as.numeric(trend::sens.slope(.data[[value_col]])$p.value),
        error = function(e) NA_real_
      ),
      !!paste0(variable, "_sd_mk_tau") := tryCatch(
        as.numeric(trend::mk.test(.data[[value_col]])$estimates["tau"]),
        error = function(e) NA_real_
      ),
      !!paste0(variable, "_sd_mk_p_value") := tryCatch(
        as.numeric(trend::mk.test(.data[[value_col]])$p.value),
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# Annual spatial heterogeneity
# -----------------------------------------------------------------------------

efa_data <- purrr::map2(
  EFA_RASTER_DIRS,
  names(EFA_RASTER_DIRS),
  process_efa
)

names(efa_data) <- names(EFA_RASTER_DIRS)

annual_efa <- purrr::reduce(
  efa_data,
  full_join,
  by = c("id", "year")
)

write.csv(
  annual_efa,
  file.path(OUTPUT_DIR, "EFA_Annual_Spatial_SD.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Temporal trends in spatial heterogeneity
# -----------------------------------------------------------------------------

target_vars <- names(EFA_RASTER_DIRS)

efa_trends <- purrr::map(
  target_vars,
  ~ calculate_trends(annual_efa, .x)
) |>
  purrr::reduce(full_join, by = "id")

# -----------------------------------------------------------------------------
# Fire-regime variables
# -----------------------------------------------------------------------------

fire_data <- terra::vect(FIRE_DATA_PATH) |>
  as.data.frame()

fire_data$id <- as.character(fire_data$id)
efa_trends$id <- as.character(efa_trends$id)

fire_keep <- fire_data |>
  select(
    id,
    cluster,
    TBA_sum,
    FF_mean,
    TSF_median,
    starts_with("sev_")
  )

master <- efa_trends |>
  left_join(fire_keep, by = "id")

severity_cols <- grep("^sev_", names(master), value = TRUE)

if (length(severity_cols) > 0) {
  master <- master |>
    rowwise() |>
    mutate(
      max_severity = {
        x <- c_across(all_of(severity_cols))
        if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
      },
      mean_severity = {
        x <- c_across(all_of(severity_cols))
        x <- x[!is.na(x) & x > 0]
        if (length(x) == 0) NA_real_ else mean(x)
      }
    ) |>
    ungroup() |>
    select(-all_of(severity_cols))
} else {
  master <- master |>
    mutate(
      max_severity = NA_real_,
      mean_severity = NA_real_
    )
}

# Keep only variables used in the final national-scale analyses.
master <- master |>
  select(
    id,
    cluster,
    all_of(EFA_TREND_VARS),
    all_of(EFA_PVARS),
    ends_with("_mk_tau"),
    ends_with("_mk_p_value"),
    all_of(FIRE_VARS)
  )

write.csv(
  master,
  file.path(OUTPUT_DIR, "Final_Master_Dataset.csv"),
  row.names = FALSE
)

# Spatial version for QGIS.
grid <- terra::vect(GRID_PATH)
master_spatial <- merge(grid, master, by = "id")
terra::writeVector(
  master_spatial,
  file.path(OUTPUT_DIR, "Final_Master_Dataset.gpkg"),
  overwrite = TRUE
)

message("National-scale master dataset written to: ", OUTPUT_DIR)
