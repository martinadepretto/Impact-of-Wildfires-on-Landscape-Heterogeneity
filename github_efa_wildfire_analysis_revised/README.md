# Wildfire and Functional Landscape Heterogeneity

R workflow for the national-scale analysis of wildfire regimes and temporal trends in spatial heterogeneity of Ecosystem Functional Attributes (EFAs) across mainland Portugal.

## Scope of the repository

This repository contains only the analyses retained in the final thesis. Exploratory analyses that were removed from the thesis are intentionally not included, including:

- continuous total-burned-area (TBA) analyses;
- spatial autocorrelation and Conley spatial regression;
- threshold-based trend analyses and sensitivity indices.

The workflow therefore focuses on the three national-scale analytical components retained in the thesis:

1. quantification of annual EFA spatial heterogeneity and temporal trends;
2. relationships between EFA trends and fire-regime variables at the national scale;
3. variation in these relationships across fire-regime clusters.

## Workflow

Run the scripts in the following order:

1. `R/00_config.R` — project paths and shared variable definitions.
2. `R/01_prepare_master_dataset.R` — calculates annual spatial standard deviation for each EFA, estimates Sen's slope and Mann-Kendall statistics, and combines the trend metrics with fire-regime variables.
3. `R/02_national_correlations.R` — calculates national-scale Spearman correlations between EFA temporal trends and fire-regime variables, with Benjamini-Hochberg correction across tests.
4. `R/03_cluster_correlations.R` — calculates cluster-level Spearman correlations, restricted to significant EFA trends and excluding Cluster 4 because of its very small sample size. Benjamini-Hochberg correction is applied within cluster.
5. `R/04_export_sen_slope_maps.R` — exports trend estimates and significance flags to a GeoPackage for cartographic work in QGIS.

## Important data-processing note

The EFA TIFF files used here should be the annual EFA rasters generated in Google Earth Engine. Each annual raster represents the temporal median of available observations calculated independently for each pixel for that calendar year. The R workflow then quantifies spatial variability among pixels within each 5 x 5 km grid cell for each year.

## Expected project structure

```text
project/
├── R/
│   ├── 00_config.R
│   ├── 01_prepare_master_dataset.R
│   ├── 02_national_correlations.R
│   ├── 03_cluster_correlations.R
│   └── 04_export_sen_slope_maps.R
├── data/
│   ├── EFA/
│   │   ├── NDVI/NDVI_clip/*.tif
│   │   ├── LST_Day/LST_Day_clip/*.tif
│   │   ├── TCT_Wetness/TCT_Wetness_clip/*.tif
│   │   └── TCT_Brightness/TCT_Brightness_clip/*.tif
│   ├── spatial/
│   │   └── Reticolo.gpkg
│   └── fire/
│       └── Results_CA.gpkg
├── outputs/
└── README.md
```

## Reproducibility notes

- User-specific absolute paths are not stored in the scripts.
- Shared paths and variable definitions are centralized in `R/00_config.R`.
- Raw datasets should not be committed to the public repository unless their licensing and data-sharing conditions allow it.
- The national-scale correlations use all grid cells with complete data for the variables being tested.
- Cluster-level correlations are restricted to cells with a statistically significant temporal EFA trend (`p < 0.05`) and exclude Cluster 4 because of its very small sample size.
- Multiple-testing correction uses the Benjamini-Hochberg procedure.
