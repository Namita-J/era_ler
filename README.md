# Spatial Extrapolation of Intercropping Performance Across Sub-Saharan Africa
Overview

This repository contains the data, scripts, and outputs for Chapter 2 of my PhD thesis. The chapter builds on the meta-analysis of Adam et al. (2025), which established how moisture variability, system composition, and agronomic practices drive intercropping performance (LER) across Africa. Chapter 2 takes those relationships and asks: where across Sub-Saharan Africa does intercropping consistently outperform sole cropping, and what environmental conditions govern this?

The analysis draws on the Evidence for Resilient Agriculture (ERA) database (Rosenstock et al., 2024) — 112,859 geolocated observations from 2,011 agricultural studies in Africa (1934–2018) — as the primary data source, and uses NEX-GDDP-CMIP6 climate data to compute SPEI and project LER suitability under current and future climate scenarios.

Research Objectives

Overall aim: To spatially predict and map the Land Equivalent Ratio (LER) of intercropping systems across Sub-Saharan Africa, identifying where intercropping consistently outperforms sole cropping and which environmental and agronomic factors govern this spatial variation.

Objective 1 — Characterise the ERA intercropping dataset and its spatial coverage
Extract and describe the intercropping subset from ERA (system types, LER distributions, geographic coverage, temporal range), assessing data sufficiency and representativeness across African agroecological zones. Identify gaps in the evidence base.

Objective 2 — Model the relationship between LER and spatially explicit environmental and agronomic predictors
Develop a predictive model linking observed LER to gridded covariates — aridity/SPEI, soil organic carbon, soil nitrogen, agroecological zone, rainfall seasonality — using ERA's geolocated observations. Model selection will consider mixed-effects regression, random forest, and geographically weighted regression approaches.

Objective 3 — Spatially extrapolate predicted LER across Sub-Saharan Africa by intercropping system type
Apply the validated model to continental-scale gridded environmental data to generate maps of predicted LER for key intercropping system categories (cereal+legume, root/tuber+legume, cereal+cereal), delineating zones of consistent intercropping advantage (LER > 1) and underperformance.

Objective 4 — Assess spatial shifts in intercropping suitability under projected climate change
Compare current LER suitability zones against NEX-GDDP-CMIP6 projections under SSP2-4.5 and SSP5-8.5 (2040–2069) to identify regions where climate change will expand or contract the advantage of intercropping over sole cropping.


Workflow
Stage 1 — Data Preparation

1.1 ERA data extraction

Update the ERA data with more recent intercropping papers 
Filter to intercropping observations only (practice codes for intercropping/mixed cropping)
Extract outcome variable: Land Equivalent Ratio (LER)
Extract covariates already in ERA: country, coordinates, study year, crop system type, soil properties, management details
Output: clean intercropping subset with LER, coordinates, and metadata

1.2 Climate data — historical SPEI

Download NEX-GDDP-CMIP6 historical daily precipitation (pr), maximum temperature (tasmax), and minimum temperature (tasmin) for 1950–2014
For each ERA observation, extract seasonal climate values matching the reported study year and location
Compute SPEI at the growing-season scale (consistent with Adam et al. methodology) using the SPEI R package or equivalent
Assign SPEI value and SPEI class to each ERA observation
Output: ERA intercropping dataset with matched historical SPEI values

1.3 Soil covariates

Download iSDAsoil or SoilGrids layers for Sub-Saharan Africa: soil organic carbon (SOC), total nitrogen (TN), pH, texture
Extract values at ERA observation coordinates
Output: ERA dataset with soil covariate columns appended

1.4 Agroecological zones

Download FAO GAEZ agroecological zone raster for Africa
Extract AEZ class at each ERA observation location
Output: ERA dataset with AEZ class column appended

1.5 Compile master modelling dataset

Join all above outputs into a single analysis-ready dataframe
Quality checks: remove observations with missing LER, implausible coordinates, or missing key covariates
Descriptive summary: sample size by system type, region, AEZ, and SPEI class
Output: data/processed/ERA_intercrop_modelling_dataset.csv
Stage 2 — Exploratory Analysis and Model Development

2.1 Descriptive and spatial analysis

Map ERA intercropping observation locations across Sub-Saharan Africa
Plot LER distributions by system type, AEZ, SPEI class, and region
Assess spatial clustering and coverage gaps
Correlate LER with candidate predictors (SPEI, SOC, TN, AEZ, system type)

2.2 Model selection

Candidate models:
Linear mixed-effects model 
Random forest: non-parametric, handles interactions and non-linearity, provides variable importance
Geographically weighted regression (GWR): tests whether predictor-LER relationships vary spatially
Model comparison via cross-validation (leave-one-study-out or k-fold), RMSE, and R²
Select final model based on predictive performance and interpretability

2.3 Model validation

Validate against held-out observations (20% split or LOOCV at study level)
Assess spatial autocorrelation in residuals (Moran's I)
Check for extrapolation beyond training data range using environmental space analysis (convex hull or MESS)
Output: model performance metrics, residual maps, validation plots
Stage 3 — Spatial Extrapolation (Current Climate)

3.1 Prepare prediction grid

Create a 0.25° × 0.25° grid (matching NEX-GDDP-CMIP6 resolution) covering Sub-Saharan Africa
Mask to land area, excluding major water bodies and the Sahara (where no relevant agricultural systems exist)

3.2 Compute baseline climatological SPEI

Use NEX-GDDP-CMIP6 historical period (1985–2014, 30-year climatology) daily pr, tasmax, tasmin
Compute mean growing-season SPEI per grid cell across the climatological period
Output: baseline SPEI raster at 0.25° resolution

3.3 Extract soil and AEZ covariates at prediction grid

Resample iSDAsoil/SoilGrids SOC and TN to 0.25° grid
Resample FAO GAEZ AEZ raster to 0.25° grid
Stack all predictor rasters

3.4 Generate LER prediction maps

Apply final model to prediction grid, for each system type separately (cereal+legume, root/tuber+legume, cereal+cereal)
Output: predicted LER rasters per system type
Derive binary suitability layer: LER > 1 (intercropping advantage) vs. LER ≤ 1
Quantify area of intercropping advantage by AEZ, country, and region
Output: outputs/maps/LER_predicted_[systemtype]_baseline.tif
Stage 4 — Future Climate Scenarios

4.1 Compute future SPEI

Download NEX-GDDP-CMIP6 projected daily pr, tasmax, tasmin for SSP2-4.5 and SSP5-8.5 (2040–2069)
Use ensemble mean across available CMIP6 models (or subset of models with good African performance)
Compute growing-season SPEI per grid cell for future period
Output: future SPEI rasters for SSP2-4.5 and SSP5-8.5

4.2 Generate future LER prediction maps

Apply final model to future predictor stacks (SPEI updated; soil and AEZ held constant)
Generate predicted LER maps for each system type × scenario combination
Output: outputs/maps/LER_predicted_[systemtype]_[scenario].tif

4.3 Change analysis

Compute difference rasters: future LER − baseline LER
Identify and map:
Regions of increasing intercropping advantage (positive shift)
Regions of decreasing advantage or emerging underperformance (negative shift)
Regions crossing the LER = 1 threshold in either direction
Summarise changes by AEZ, country, and region
Output: change maps and summary tables
