# Spatial Extrapolation of Intercropping Performance Across Sub-Saharan Africa
Overview

This repository contains the data, scripts, and outputs for Chapter 2 of my PhD thesis. The chapter builds on the meta-analysis of Adam et al. (2025), which established how moisture variability, system composition, and agronomic practices drive intercropping performance (LER) across Africa. Chapter 2 takes those relationships and asks: where across Sub-Saharan Africa does intercropping consistently outperform sole cropping, and what environmental conditions govern this?

The analysis draws on the Evidence for Resilient Agriculture (ERA) database (Rosenstock et al., 2024) - 112,859 geolocated observations from 2,011 agricultural studies in Africa (1934–2018) - as the primary data source, and uses NEX-GDDP-CMIP6 climate data to compute SPEI and project LER suitability under current and future climate scenarios.

Research Objectives

Overall aim: To spatially predict and map the Land Equivalent Ratio (LER) of intercropping systems across Sub-Saharan Africa, identifying where intercropping consistently outperforms sole cropping and which environmental and agronomic factors govern this spatial variation.

Objective 1 — Characterise the ERA intercropping dataset and its spatial coverage
Extract and describe the intercropping subset from ERA (system types, LER distributions, geographic coverage, temporal range), assessing data sufficiency and representativeness across African agroecological zones. Identify gaps in the evidence base.

Objective 2 — Map environmental similarity and extrapolation confidence across Sub-Saharan Africa
Compute Multivariate Environmental Similarity Surfaces (MESS) to quantify, for every grid cell across Sub-Saharan Africa, how similar local environmental conditions are to the range of conditions represented in the ERA intercropping trial observations. Identify regions of high extrapolation risk (novel environments) where predictions should be interpreted with caution.

Objective 3 — Spatially extrapolate LER using analogue matching in environmental space
For each grid cell, identify the k most environmentally similar ERA intercropping trial observations using distance-weighted k-nearest neighbour matching (Gower distance across SPEI, soil, and temperature covariates) and predict LER as the similarity-weighted mean of analogue observations. Generate system-specific LER maps (cereal+legume, root/tuber+legume, cereal+cereal) delineating zones of consistent intercropping advantage (LER > 1) and underperformance, paired with prediction uncertainty maps derived from analogue variance.

Objective 4 — Identify the dominant environmental drivers of spatial LER variation
Using a complementary mixed-effects regression and variable importance analysis (random forest), quantify the relative contribution of aridity/SPEI, soil organic carbon, soil nitrogen, temperature, and agroecological zone to LER variation.

Objective 5 — Assess spatial shifts in intercropping suitability under projected climate change
Recompute environmental similarity and analogue-based LER predictions using NEX-GDDP-CMIP6 projections under SSP2-4.5 and SSP5-8.5 (2040–2069), identifying regions where climate change will expand or contract the advantage of intercropping over sole cropping.


Workflow
Stage 1 - Data Preparation

1.1 ERA data extraction

- Update the ERA data with more recent intercropping papers 

- Filter to intercropping observations only (practice codes for intercropping/mixed cropping)

- Output: clean intercropping subset with LER, coordinates, and metadata

1.2 Climate data - historical SPEI

- Download NEX-GDDP-CMIP6 historical daily precipitation (pr), maximum temperature (tasmax), and minimum temperature (tasmin) for 1950–2014

- For each ERA observation, extract seasonal climate values matching the reported study year and location

- Compute SPEI at the growing-season scale (consistent with Adam et al. methodology) using the SPEI R package or equivalent

- Assign SPEI value and SPEI class to each ERA observation

- Output: ERA intercropping dataset with matched historical SPEI values

1.3 Soil covariates

- Download iSDAsoil or SoilGrids layers for Sub-Saharan Africa: soil organic carbon (SOC), total nitrogen (TN), pH, texture

- Extract values at ERA observation coordinates

- Output: ERA dataset with soil covariate columns appended

1.4 Compile master modelling dataset

- Join all above outputs into a single analysis-ready dataframe

- Quality checks: remove observations with missing LER, implausible coordinates, or missing key covariates

- Descriptive summary: sample size by system type, region, AEZ, and SPEI class
  
- Output: final datasey
  
Stage 2 — Exploratory Analysis

2.1 Descriptive and spatial analysis

- Map ERA intercropping observation locations across Sub-Saharan Africa
- Plot LER distributions by system type, AEZ, SPEI class, and region
- Assess spatial clustering and coverage gaps
- Correlate LER with candidate predictors (SPEI, SOC, TN, AEZ, temperature, system type)

2.2 Environmental space characterisation

Define the multivariate environmental space occupied by ERA intercropping trial locations across the following dimensions:
- Mean growing-season SPEI
- Aridity index (AI)
- Soil organic carbon (SOC)
- Total nitrogen (TN)
- Mean growing-season temperature
- Rainfall seasonality (CV of monthly rainfall)
- AEZ class
- Visualise coverage of environmental space
- Identify underrepresented environmental conditions 

Stage 3 — Similarity-Based LER Prediction

3.1 Compute baseline similarity surface

- Compute a multivariate environmental similarity surface across the Sub-Saharan Africa prediction grid
- Reference point set: ERA intercropping trial locations with their extracted environmental covariate values
- Prediction grid: 0.25° × 0.25° grid covering Sub-Saharan Africa, masked to agricultural land
- Similarity is computed for each predictor variable and summarised as a minimum similarity score across all variables
- Positive values indicate conditions within the range of the ERA training data; negative values indicate novel environments beyond the training range
Output: continuous similarity raster and binary novel/non-novel mask

3.2 Map extrapolation confidence zones

- Classify similarity output into confidence tiers: high confidence (well-represented in ERA trial data), moderate confidence (within range but sparsely represented), and low confidence/novel (beyond the training range)
- Map confidence zones across Sub-Saharan Africa and overlay ERA trial locations to visualise spatial evidence density

3.3 Analogue-based LER prediction

- For each grid cell in the prediction grid, identify the k most environmentally similar ERA intercropping trial observations using a distance metric computed across all environmental covariates
- Predict LER as the distance-weighted mean of the k analogues, with closer analogues receiving higher weight
- Compute prediction uncertainty as the distance-weighted standard deviation across analogues
- Run separately for each intercropping system type (cereal+legume, root/tuber+legume, cereal+cereal) and for all systems combined
- Predictions in low-confidence/novel zones are explicitly flagged in all outputs
- Output: predicted mean LER raster, prediction uncertainty raster, and binary LER > 1 suitability layer per system type

3.4 Analogue diagnostics

- For a sample of grid cells, report the identity and characteristics of top analogues (which ERA studies, which countries, which environmental conditions) to validate that matches are agronomically plausible
- Map mean analogue distance across the prediction grid as a spatial diagnostic of prediction support

3.5 Future scenario prediction

- Repeat similarity computation and analogue-based LER prediction using NEX-GDDP-CMIP6 future climate covariates (SSP2-4.5 and SSP5-8.5, 2040–2069); soil and agroecological zone covariates held constant
- Identify regions where future climate conditions become increasingly novel relative to the current ERA evidence base
- Compute difference rasters (future LER − baseline LER) and map regions of increasing or decreasing intercropping advantage, including areas crossing the LER = 1 threshold in either direction
- Summarise changes by agroecological zone, country, and region

5.1 Mixed-effects model

- Fit a linear mixed-effects model-
- Use to derive predictor importance and partial dependence relationships
- Compare spatial predictions to analogue-based maps as a sensitivity check

5.2 Random forest

- Fit a random forest model on ERA intercropping observations
- Extract variable importance scores
- Generate partial dependence plots for key predictors (SPEI, SOC, AEZ)
- Use proximity matrix to validate environmental similarity structure used in Stage 4

