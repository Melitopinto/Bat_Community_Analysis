# Bat Community Analysis
Analysis of insectivorous bat species-specific activity and community diversity.

## Overview
This repository contains the analytical pipeline for modeling bat species-specific activity and community diversity in an urban landscape. The analytical pipeline is based off Cole et al. 2022 
EcoCountHelper: an R package and analytical pipeline for the analysis of ecological count data using GLMMs, and a case study of bats in Grand Teton National Park

## Taxonomic Note
**Important:** Please note that the taxonomic nomenclature used in this analysis reflects updated classifications, specifically regarding the genus *Alionoctula*. 

## Data
The analysis uses the following raw data files:
* `Raw data_bat passes.csv`: Contains the primary acoustic monitoring data.
* `Sampled_locations.csv`: Metadata regarding site characteristics.
* `Season wise bat community data.csv`: Aggregated community data.

## Analytical Pipeline
The code is written in R and follows the `EcoCountHelper` workflow for ecological count data using Generalized Linear Mixed Models (GLMMs).

## How to use
1. Ensure all data files are present in the working directory.
2. Open `Bat_Community_Analysis.Rproj` in RStudio.
3. Run the script `Urban_bats_analysis.R` to reproduce the statistical models and visualizations.
