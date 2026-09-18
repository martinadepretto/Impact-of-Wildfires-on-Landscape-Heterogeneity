# Impact-of-Wildfires-on-Landscape-Heterogeneity
Assessing the Impact of Wildfires on Landscape Heterogeneity — a Multi-Indicator Remote Sensing Analysis using Ecosystem Functioning Attributes for mainland Portugal

This repository contains the R scripts used to analyse the relationship between wildfire regimes and the spatial heterogeneity of ecosystem functioning across mainland Portugal.

The workflow corresponds to the national-scale analysis presented in Section 2 of the thesis and is based on four remotely sensed Ecosystem Functional Attributes (EFAs):

Normalized Difference Vegetation Index (NDVI)
Land Surface Temperature (LST)
Tasseled Cap Brightness (TCTB)
Tasseled Cap Wetness (TCTW)

Annual EFA raster datasets were generated in Google Earth Engine by calculating the temporal median of available observations for each pixel within each calendar year from 2000 to 2024. The resulting annual rasters were subsequently analysed in R.


# File 1: 1_cluster_analysis.R

Contains the central configuration used throughout the workflow.

This script defines:

input and output directories;
analysis years;
EFA names;
fire-regime variables;
spatial grid information;
statistical parameters.

It does not perform the main analyses and should be run or sourced before the other scripts.
