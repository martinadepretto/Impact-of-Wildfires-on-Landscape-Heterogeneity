# Fire Regime Cluster Analysis – Mainland Portugal

## Overview

This repository contains the R workflow used to characterize fire-regime contexts across mainland Portugal through k-means clustering.

The analysis combines fire-related and land-cover variables to identify spatially distinct fire-regime clusters. The selected clustering solution is then characterized using cluster-level summaries and exported datasets suitable for further analysis and visualization in R and QGIS.

## Analysis workflow

The script performs the following steps:

1. Loads the input dataset.
2. Defines the variables used for clustering:
   - Fire frequency (`FF_mean`)
   - Total burned area (`TBA_sum`)
   - Time since fire (`TSF_median`)
   - Forest cover (`Forest_per`)
   - Grassland cover (`Grass_perc`)
   - Shrubland cover (`Shrub_perc`)
   - Available severity variables (`sev*`)
3. Replaces missing values with zero for the variables used in the clustering and severity calculations, following the assumptions of the original analysis.
4. Standardizes all clustering variables.
5. Evaluates candidate numbers of clusters using:
   - Elbow method
   - Silhouette method
6. Applies k-means clustering using the selected solution (`k = 6`).
7. Calculates maximum recorded fire severity for each spatial unit.
8. Characterizes the six clusters using fire-regime and land-cover variables.
9. Assigns descriptive labels to the identified fire-regime clusters.
10. Produces comparative cluster visualizations and a parallel-coordinates plot.
11. Exports cluster statistics, variable rankings, the final QGIS-ready dataset, and the fitted k-means objects.

The script is intended to reproduce the fire-regime classification used as the spatial context for the subsequent ecosystem functional attribute analyses.

## Repository structure

```text
fire_regime_cluster_analysis/
├── data/
│   └── 7_variables_EU_completed_table.csv
├── outputs/
├── scripts/
│   └── cluster_analysis.R
└── README.md
```

The `data/` and `outputs/` directories are included as placeholders. The input dataset should be placed in `data/` before running the script. Generated results are written to `outputs/`.

## Requirements

The analysis was developed in R and requires the following packages:

- `tidyverse`
- `factoextra`
- `cluster`
- `openxlsx`
- `GGally`

The script checks for missing packages and installs them automatically when necessary.

## Input data

The expected input file is:

```text
data/7_variables_EU_completed_table.csv
```

The dataset must contain the variables required by the script, including the fire-regime, land-cover, and severity columns described above.

Severity variables are identified automatically as columns whose names start with `sev`.

## Reproducibility settings

The clustering uses:

- Number of clusters: `k = 6`
- Random seed: `123`
- Number of random starts: `25`
- Maximum iterations: `100`

These settings are defined at the beginning of the script and can be modified if needed.

## Outputs

The script generates the following main outputs in `outputs/`:

### Diagnostic figures

- `cluster_diagnostic_elbow.png`
- `cluster_diagnostic_silhouette.png`

### Cluster visualizations

- `cluster_pca_visualization.png`
- `comparative_fire_regimes.png`
- `fire_regime_parallel_coordinates.png`

### Cluster statistics

- `Statistics_Clusters.csv`
- `Cluster_Rankings_By_Variable.xlsx`

### QGIS-ready dataset

- `Results_CA_final.csv`

This file contains the original spatial-unit attributes together with the assigned cluster ID, maximum severity, and descriptive fire-regime label.

### R objects

- `kmeans_k6_model.rds`
- `cluster_statistics.rds`

## Fire-regime labels

The selected six clusters are labelled in the script as:

1. High-Recurrence, High-Severity Regime
2. Dynamic Forest-Dominated Regime
3. High-Severity, Low-Recurrence Regime
4. Shrubland-Dominated Hotspots
5. Stable, Low-Fire Regime
6. Moderate-Intensity Transition Regime

These labels are descriptive summaries of the cluster characteristics and are assigned after the k-means solution is calculated.

## Running the analysis

From the project root, run:

```r
source("scripts/cluster_analysis.R")
```

Alternatively, open `scripts/cluster_analysis.R` in RStudio and run the script from top to bottom.

Before running the analysis, make sure that the input CSV is available at:

```text
data/7_variables_EU_completed_table.csv
```

## Notes on interpretation

The clustering is based on standardized variables because the input variables are measured on different scales. The resulting clusters represent multivariate fire-regime contexts rather than individual fire characteristics.

The cluster labels should therefore be interpreted as descriptive summaries of the multivariate profiles identified by the analysis.

## Author

This repository contains the analysis code developed for the thesis on wildfire and ecosystem functional spatial heterogeneity in mainland Portugal.
