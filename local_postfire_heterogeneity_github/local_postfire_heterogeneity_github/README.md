# Local Post-Fire Functional Spatial Heterogeneity – 2005 Wildfire

## Overview

This repository contains the R workflow used for the local-scale analysis of
post-fire functional spatial heterogeneity in selected areas affected by the
2005 wildfire in mainland Portugal.

The workflow is divided into two main scripts:

1. **Local focal-site selection and spatial data preparation**
2. **Post-fire spatial heterogeneity analysis**

The focal sites are intentionally selected to represent contrasting long-term
NDVI trajectories within fire-regime and burned-area strata. They are therefore
not a random sample of the 2005 burned landscape.

## Workflow

### 1. Local focal-site selection

`3.1_Local_Cells_Selection.R`:

- reconstructs contiguous 2005 fire scars;
- excludes areas that burned again after 2005;
- calculates total undisturbed burned area for each reconstructed fire scar;
- assigns the dominant fire-regime cluster;
- assigns burned-area size classes;
- extracts annual NDVI values for 2000–2024;
- calculates pre-fire (2000–2004) and post-fire (2020–2024) NDVI values;
- identifies contrasting NDVI trajectories within fire-regime and burned-area strata;
- summarizes contemporary COS2023 land-cover composition;
- exports focal-site spatial layers and tabular data for subsequent analysis.

The script explicitly performs fire-size classification at the reconstructed
fire-scar level so that fragments belonging to the same fire receive the same
size class.

### 2. Post-fire functional heterogeneity

`3.2_Postfire_Heterogeneity_Analysis.R`:

- assembles the focal-site analysis dataset;
- calculates changes in EFA spatial dispersion using MAD;
- calculates relative changes in MAD;
- summarizes pre-fire and post-fire heterogeneity by EFA;
- tests the relationship between total burned area and changes in spatial
  dispersion using Spearman rank correlations;
- evaluates the relationship between changes in EFA mean and changes in MAD;
- calculates the exploratory Structural Sensitivity Ratio (SSR);
- produces the main trajectory and land-cover contextual figures;
- exports site-level and summary tables.

The four ecosystem functional attributes are:

- NDVI
- Land Surface Temperature (LST)
- Tasseled Cap Brightness (TCTB)
- Tasseled Cap Wetness (TCTW)

## Key methodological definitions

### Annual EFA data

The analysis uses annual EFA raster datasets generated previously from the
remote-sensing workflow. Each annual raster represents the pixel-level median
of the available observations for that year.

### Spatial dispersion

Within each focal site, spatial dispersion is quantified using the
**median absolute deviation (MAD)** across pixels.

The main response variable is:

```text
ΔMAD = post-fire MAD − pre-fire MAD
```

Relative change is calculated as:

```text
%ΔMAD = 100 × ΔMAD / |pre-fire MAD|
```

### Functional-state change

Changes in the EFA central tendency are calculated as:

```text
ΔMean = post-fire mean − pre-fire mean
```

The analysis therefore distinguishes changes in the average functional state
from changes in the internal spatial organization of the site.

### Fire extent

`total_fire_area_ha` represents the total spatial extent of the undisturbed
2005 burned area for each focal site. It is used as a measure of wildfire
extent and **not as a measure of fire severity**.

### Structural Sensitivity Ratio

SSR is included as an exploratory descriptive metric:

```text
SSR = |ΔMAD| / |ΔMean|
```

A relative version is also calculated from the corresponding percentage
changes. SSR values can become unstable when the denominator approaches zero
and should therefore not be interpreted as a causal sensitivity measure.

## Land-cover context

Land-cover information is derived from COS2023.

COS2023 is used as **contemporary contextual information** and does not represent
the land-cover composition of the landscape in 2005. Land-cover class richness
and dominant-class composition are therefore interpreted as contextual
associations rather than direct pre-fire explanatory variables.

## Repository structure

```text
local_postfire_heterogeneity/
├── README.md
├── .gitignore
├── scripts/
│   ├── 3.1_Local_Cells_Selection.R
│   └── 3.2_Postfire_Heterogeneity_Analysis.R
├── data/
└── outputs/
    ├── Figures/
    └── Tables/
```

The `data/` and `outputs/` directories are included as placeholders. Large
input rasters, geospatial datasets and generated outputs should generally not
be committed to GitHub unless they are small enough and their redistribution
is permitted.

## Main input data

The scripts expect locally available datasets including:

- the fire-grid intersection dataset;
- annual NDVI rasters for 2000–2024;
- COS2023 land-cover data;
- the focal-site and heterogeneity objects required by the second script.

The paths to these datasets are defined in the configuration sections of the
scripts and should be adapted to the local analysis environment before running
the workflow.

## Required R packages

### Script 3.1

- `sf`
- `terra`
- `dplyr`
- `tidyr`
- `stringr`
- `ggplot2`

### Script 3.2

- `dplyr`
- `tidyr`
- `ggplot2`
- `readr`
- `ggrepel`
- `patchwork`

## Running the workflow

Run the scripts in sequence.

From the repository root:

```r
source("scripts/3.1_Local_Cells_Selection.R")
source("scripts/3.2_Postfire_Heterogeneity_Analysis.R")
```

The second script depends on objects generated during the preceding local EFA
processing workflow. In particular, it expects the required heterogeneity
summary and focal-site area information to be available in the R environment.

## Main outputs

### Focal-site selection

The first script produces:

- `Selected_Focal_Sites.csv`
- `Focal_Sites_COS2023_Composition.csv`
- `Focal_Sites_COS2023_Summary.csv`
- `Selected_Focal_Sites.gpkg`
- `Selected_Focal_Sites.kml`
- `Selected_Focal_Sites_GEE.shp`
- individual focal-site shapefiles for GEE
- `Focal_Site_NDVI_Selection.png`

### Post-fire heterogeneity analysis

The second script produces:

- `Site_Level_Postfire_Heterogeneity.csv`
- `Heterogeneity_Summary_By_EFA.csv`
- `Heterogeneity_Trajectory_Summary.csv`
- `Fire_Area_Heterogeneity_Relationship.csv`
- `Delta_Mean_Delta_MAD_Correlations.csv`
- `Structural_Sensitivity_Ratio_Summary.csv`

and the main figures:

- `Fire_Area_vs_Delta_MAD.png`
- `Fire_Area_vs_Percent_Delta_MAD.png`
- `Delta_Mean_vs_Delta_MAD.png`
- `Mean_MAD_Trajectory_Concordance.png`
- `Landcover_Class_Richness_MAD.png`

Optional land-cover tables are exported when the corresponding objects are
available in the R environment.

## Interpretation

The local analysis is designed to identify and compare contrasting focal-site
trajectories rather than estimate an unbiased regional average response.
Associations are therefore interpreted as observational relationships within
the selected sites.

In particular:

- burned area is treated as wildfire extent;
- ΔMAD is the primary measure of change in spatial dispersion;
- %ΔMAD provides a relative measure;
- ΔMean and ΔMAD represent distinct dimensions of functional response;
- SSR is exploratory;
- COS2023 land-cover variables provide contemporary context rather than a
  reconstruction of pre-fire land cover.

## Reproducibility

The scripts contain the configuration parameters and analytical steps required
to reproduce the workflow, provided that the required input datasets and the
preceding EFA-processing objects are available.

Paths to local data sources should be adapted to the user's analysis
environment before execution.

## Author

This repository contains the analysis code developed for the thesis on
wildfire and ecosystem functional spatial heterogeneity in mainland Portugal.
