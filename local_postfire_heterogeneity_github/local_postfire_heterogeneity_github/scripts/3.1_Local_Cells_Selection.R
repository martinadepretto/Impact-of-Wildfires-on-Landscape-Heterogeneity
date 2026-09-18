# ==============================================================================
# Local focal-site selection and spatial data preparation
# ==============================================================================
#
# Purpose
# -------
# Identify focal areas within the 2005 wildfire footprint for the local-scale
# analysis. The workflow:
#   1. reconstructs contiguous 2005 fire scars;
#   2. excludes grid cells that burned again after 2005;
#   3. assigns each fire scar a dominant fire-regime cluster and size class;
#   4. extracts annual NDVI values and selects contrasting trajectories;
#   5. summarizes contemporary COS2023 land-cover composition; and
#   6. exports spatial layers for QGIS and Google Earth Engine.
#
# The focal sites are intentionally selected to represent contrasting NDVI
# trajectories within fire-regime and burned-area strata. They are therefore
# not a random sample of the 2005 burned landscape.
# ==============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

fire_layer_path <- "C:/fires/Portugal/total_fires_grid_intersecation.gpkg"
ndvi_dir <- "C:/fires/Portugal/EFA/Local_Analysis/NDVI/NDVI_median_year/"
cos_path <- "C:/fires/Portugal/land_cover/land_cover_DGT/COS2023v1-S2.gpkg"

study_years <- 2000:2024
pre_years <- 2000:2004
post_years <- 2020:2024
fire_year <- 2005
min_fire_fragment_area_ha <- 10
n_size_classes <- 5

output_dir <- "Local_Analysis_Output"
ggee_dir <- file.path(output_dir, "GEE_Individual_Assets")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ggee_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Load and clean the fire-grid intersection
# ------------------------------------------------------------------------------

message("Loading fire-grid intersection...")

fire_grid <- read_sf(fire_layer_path) %>%
  mutate(
    year = as.integer(year),
    id = as.character(id),
    cluster = as.character(cluster),
    area_ha = as.numeric(gsub(",", ".", as.character(area_ha)))
  ) %>%
  filter(area_ha >= min_fire_fragment_area_ha)

# ------------------------------------------------------------------------------
# 3. Reconstruct 2005 fire scars and remove re-burned cells
# ------------------------------------------------------------------------------

message("Reconstructing contiguous 2005 fire scars...")

fires_2005 <- fire_grid %>% filter(year == fire_year)

fire_scars <- fires_2005 %>%
  st_union(by_feature = FALSE) %>%
  st_cast("POLYGON") %>%
  st_as_sf() %>%
  mutate(global_fire_id = paste0("Fire2005_", row_number()))

fire_fragments <- st_intersection(fires_2005, fire_scars)

reburned_ids <- fire_grid %>%
  filter(year > fire_year) %>%
  distinct(id) %>%
  pull(id)

undisturbed_2005 <- fire_fragments %>%
  filter(!id %in% reburned_ids)

# Total undisturbed burned area and dominant fire-regime cluster are calculated
# once for each reconstructed fire scar.
fire_metrics <- undisturbed_2005 %>%
  st_drop_geometry() %>%
  group_by(global_fire_id) %>%
  summarise(
    total_fire_area_ha = sum(area_ha, na.rm = TRUE),
    .groups = "drop"
  )

fire_clusters <- undisturbed_2005 %>%
  st_drop_geometry() %>%
  group_by(global_fire_id, cluster) %>%
  summarise(cluster_area_ha = sum(area_ha, na.rm = TRUE), .groups = "drop") %>%
  group_by(global_fire_id) %>%
  slice_max(cluster_area_ha, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(global_fire_id, majority_cluster = cluster)

# Size classes are assigned at the fire-scar level rather than to individual
# fragments so that all fragments belonging to the same fire receive the same class.
fire_strata <- fire_metrics %>%
  left_join(fire_clusters, by = "global_fire_id") %>%
  mutate(
    size_quantile = ntile(total_fire_area_ha, n_size_classes),
    size_class = factor(
      size_quantile,
      levels = seq_len(n_size_classes),
      labels = c("Very_Low", "Low", "Medium", "High", "Very_High")
    ),
    cluster = majority_cluster,
    stratum = paste0("Cluster_", cluster, "__Size_", size_class)
  )

stratified_fragments <- undisturbed_2005 %>%
  left_join(fire_strata, by = "global_fire_id", suffix = c("", "_fire")) %>%
  mutate(
    cluster = majority_cluster,
    stratum = paste0("Cluster_", cluster, "__Size_", size_class)
  )

message("Stratification matrix:")
print(table(fire_strata$cluster, fire_strata$size_class))

# ------------------------------------------------------------------------------
# 4. Extract annual NDVI and quantify long-term trajectories
# ------------------------------------------------------------------------------

message("Loading annual NDVI rasters...")

ndvi_files <- list.files(
  ndvi_dir,
  pattern = "^LandsatH_NDVI_Median_.*\\.tif$",
  full.names = TRUE
)

ndvi_years <- as.integer(str_extract(basename(ndvi_files), "\\d{4}"))
ndvi_files <- ndvi_files[order(ndvi_years)]
ndvi_years <- ndvi_years[order(ndvi_years)]

if (!identical(sort(ndvi_years), study_years)) {
  warning("Detected NDVI years differ from the expected 2000:2024 range.")
}

ndvi_stack <- rast(ndvi_files)

if (!terra::same.crs(vect(stratified_fragments), ndvi_stack)) {
  stratified_fragments <- st_transform(stratified_fragments, st_crs(ndvi_stack))
}

ndvi_values <- terra::extract(
  ndvi_stack,
  vect(stratified_fragments),
  fun = median,
  na.rm = TRUE
)

names(ndvi_values)[-1] <- paste0("Y_", ndvi_years)

ndvi_by_fragment <- stratified_fragments %>%
  mutate(extraction_id = row_number()) %>%
  left_join(ndvi_values, by = c("extraction_id" = "ID"))

ndvi_response <- ndvi_by_fragment %>%
  st_drop_geometry() %>%
  rowwise() %>%
  mutate(
    ndvi_pre_fire = median(c_across(all_of(paste0("Y_", pre_years))), na.rm = TRUE),
    ndvi_post_fire = median(c_across(all_of(paste0("Y_", post_years))), na.rm = TRUE),
    delta_ndvi = ndvi_post_fire - ndvi_pre_fire
  ) %>%
  ungroup()

# Each fire scar can contain multiple grid fragments. Site selection is performed
# at the grid-cell level within each fire-stratum, matching the original workflow.
selected_focal_sites <- ndvi_response %>%
  group_by(stratum) %>%
  filter(
    delta_ndvi == min(delta_ndvi, na.rm = TRUE) |
      delta_ndvi == max(delta_ndvi, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  distinct(id, .keep_all = TRUE) %>%
  select(
    id, global_fire_id, cluster, size_class, stratum,
    total_fire_area_ha, ndvi_pre_fire, ndvi_post_fire, delta_ndvi
  ) %>%
  arrange(stratum, delta_ndvi)

message("Selected focal sites: ", nrow(selected_focal_sites))
print(selected_focal_sites)

# Preserve the selected spatial geometries.
selected_sites_sf <- stratified_fragments %>%
  filter(id %in% selected_focal_sites$id) %>%
  select(id, global_fire_id, cluster, size_class, stratum, total_fire_area_ha) %>%
  group_by(id) %>%
  summarise(
    across(c(global_fire_id, cluster, size_class, stratum, total_fire_area_ha), first),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 5. Plot the selected NDVI trajectories
# ------------------------------------------------------------------------------

selection_plot_data <- selected_focal_sites %>%
  mutate(
    site_label = paste0("ID: ", id, " (", size_class, ")"),
    response = if_else(
      delta_ndvi < 0,
      "Lower post-fire NDVI",
      "Higher post-fire NDVI"
    )
  ) %>%
  arrange(stratum, delta_ndvi) %>%
  mutate(site_label = factor(site_label, levels = unique(site_label)))

selection_plot <- ggplot(selection_plot_data, aes(site_label, delta_ndvi, fill = response)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Long-term NDVI trajectories of selected focal sites",
    subtitle = "Extreme trajectories selected within fire-regime and fire-size strata",
    x = "Focal site",
    y = expression(Delta * " NDVI"),
    fill = "Trajectory"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(selection_plot)
ggsave(
  file.path(output_dir, "Focal_Site_NDVI_Selection.png"),
  selection_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 6. Summarise contemporary COS2023 land-cover composition
# ------------------------------------------------------------------------------
# COS2023 is used as contemporary contextual information. It does not represent
# the pre-fire land-cover composition in 2005.

message("Intersecting selected sites with COS2023...")

cos <- read_sf(cos_path, layer = "COS2023v1", quiet = TRUE)

if (st_crs(selected_sites_sf) != st_crs(cos)) {
  cos <- st_transform(cos, st_crs(selected_sites_sf))
}

site_cos <- st_intersection(selected_sites_sf, cos)

focal_landcover <- site_cos %>%
  mutate(fragment_area_ha = as.numeric(st_area(.)) / 10000) %>%
  st_drop_geometry() %>%
  group_by(id, cluster, size_class, COS23_n4_C) %>%
  summarise(area_ha = sum(fragment_area_ha, na.rm = TRUE), .groups = "drop") %>%
  group_by(id) %>%
  mutate(
    area_prop = area_ha / sum(area_ha),
    area_pct = 100 * area_prop
  ) %>%
  ungroup() %>%
  arrange(id, desc(area_pct))

landcover_summary <- focal_landcover %>%
  group_by(id) %>%
  summarise(
    n_landcover_classes = n_distinct(COS23_n4_C),
    shannon_index = -sum(area_prop * log(area_prop), na.rm = TRUE),
    dominant_class = COS23_n4_C[which.max(area_prop)],
    dominant_class_pct = max(area_pct, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 7. Export spatial layers and tabular outputs
# ------------------------------------------------------------------------------

message("Exporting focal-site layers and tables...")

write_csv <- readr::write_csv

write_csv(
  selected_focal_sites,
  file.path(output_dir, "Selected_Focal_Sites.csv")
)

write_csv(
  focal_landcover,
  file.path(output_dir, "Focal_Sites_COS2023_Composition.csv")
)

write_csv(
  landcover_summary,
  file.path(output_dir, "Focal_Sites_COS2023_Summary.csv")
)

# QGIS layer with one geometry per focal site.
st_write(
  selected_sites_sf,
  file.path(output_dir, "Selected_Focal_Sites.gpkg"),
  layer = "focal_sites",
  driver = "GPKG",
  delete_dsn = TRUE,
  quiet = TRUE
)

# KML with a compact land-cover description for visual inspection.
kml_description <- focal_landcover %>%
  group_by(id) %>%
  summarise(
    description = paste(
      paste0(COS23_n4_C, ": ", round(area_pct, 1), "%"),
      collapse = "<br>"
    ),
    .groups = "drop"
  )

selected_sites_kml <- selected_sites_sf %>%
  left_join(kml_description, by = "id") %>%
  st_transform(4326) %>%
  transmute(
    Name = paste0("Fire site ", id),
    Description = description
  )

st_write(
  selected_sites_kml,
  file.path(output_dir, "Selected_Focal_Sites.kml"),
  driver = "KML",
  delete_dsn = TRUE,
  quiet = TRUE
)

# GEE-ready shapefile containing one feature per focal site.
ggee_sites <- selected_sites_sf %>%
  st_cast("MULTIPOLYGON") %>%
  st_transform(4326)

ggee_shp <- file.path(output_dir, "Selected_Focal_Sites_GEE.shp")
st_write(ggee_sites, ggee_shp, driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE)

# Optional individual assets for site-specific GEE tasks.
for (site_id in unique(ggee_sites$id)) {
  st_write(
    filter(ggee_sites, id == site_id),
    file.path(ggee_dir, paste0("fire_site_", site_id, ".shp")),
    driver = "ESRI Shapefile",
    delete_layer = TRUE,
    quiet = TRUE
  )
}

message("Local focal-site preparation completed successfully.")
