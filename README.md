# Temp_Biodiv_Change

This repository contains code necessary to replicate data analysis, figures, and tables for the manuscript "Temperature-related biodiversity change across temperate marine and terrestrial systems".

It will be archived on Zenodo upon manuscript publication.

We quantify local biodiversity change in assemblage-level time series in relation to temperature change, across the temperate regions of the globe, and comparing marine and terrestrial realms. We use the BioTIME database [Dornelas et al. 2018](https://doi.org/10.1111/geb.12729) to estimate rates of change in species richness, total abundance, species gains and losses over time and across taxonomic groups. For the same locations and periods, we estimate temperature change trends. We then evaluate biodiversity responses related to temperature change using meta-analytical Bayesian models, accounting for data structure.

Specifically, we test two predictions for the effects of temperature change on assemblage-level diversity within temperate regions:
1. Species richness and total abundance will increase with warming, and such increases will be greatest across relatively warm regions that border the species-rich tropics;
2. The coupling of assemblage and temperature change will be tighter in the ocean than on land.



**Disclaimer:** The code in this repository represents one version of the code developed for the project, and may yet undergo changes and revisions.

**Authors:**  This work was developed through collaboration with Maria Dornelas, Aafke Schipper, Amanda Bates, Shane Blowes, Conor Waldock, Sarah Supp and Anne Magurran. Maria Dornelas and Aafke Schipper are senior authors on the paper.


**Contacts:** 
* Laura Antao - laura.antao@helsinki.fi


## Data 

Biodiversity data
The time series analysed were from references found in the BioTIME dataset and in other studies which were used with permission.

Approximately 92% (306 references) of the biodiversity studies analysed here are available as part of the published BioTIME Database28. The data are openly available, and can be accessed on Zenodo (https://doi.org/10.5281/zenodo.1211105) or through the BioTIME website (http://biotime.st-andrews.ac.uk/).

Dornelas, M., L.H. Antao, F. Moyes, A.E. Bates, A.E. Magurran, and BioTIME consortium (200+ authors). 2018. BioTIME: a database of biodiversity time-series for the Anthropocene. Global Ecology and Biogeography. 10.1111/geb.12729 

The remaining 8% (26 references) of biodiversity studies analysed were used with permission. Some of these studies are published and publicly available outside of the BioTIME database, and others are available with permission from the corresponding author on reasonable request. For more details regarding the original citations, corresponding authors, and availability of these datasets, please refer to Table S1 in Dornelas et al. (2018). 


The data analysed in this project originates from a series of steps described in detail in https://github.com/sChange-workshop/BioGeo-BioDiv-Change, and also in the manuscript [Blowes et al. 2019](https://doi.org/10.1126/science.aaw1620). First, the studies in BioTIME underwent a spatial harmonisation process, specifically being gridded into equal-extent grid cells. This means that each individual sample was allocated to a specific combination of study ID and grid cell based on its latitude and longitude. Each new time series contained samples from only one study, and thus the integrity of each sampling method within each study was maintained.
Then, for each time series samples were aggregated within years, and the coverage (or sample completeness) for each cell-year combination was calculated. All cell-year combinations with coverage < 0.85 were excluded. Additionally, we applied sample-based rarefaction to the time series, to equalise sampling effort over time, and calculated diversity metrics as the median of 199 bootstrap rarefactions. We then calculated linear slopes over time for species richness, total abundance, species gains and losses. We retained the estimated slope and standard error for each time series for use in our second-stage meta-analytic models.

For this project, we selected only temperate time series (absolute latitudinal cut-offs at 60° and 23.5°) with a minimum of five years sampled.


Temperature data
We extracted temperature records from HadCRUT4 (Jones et al. 2012; Harris et al. 2014), specifically the HadSST3 data for marine Sea Surface Temperature (SST) on a 1° resolution, and the CRUTEM4 data for air temperature on land on a 0.5° resolution. For the location of each biodiversity time series, we extracted monthly mean temperature records for the duration of the biodiversity monitoring period (Year_start:Year_end), and estimated mean temperature trends using generalized additive models (GAM).

We further extracted annual mean temperature data from the WorldClim (Fick & Hijmans 2017) database for terrestrial time series, and from the Bio-ORACLE database (Tyberghein et al. 2012; Assis et al. 2018) for marine time series (on a resolution of 0.01° for terrestrial and of 0.1° for marine systems). Finally, we extracted the variables "Mean Temperature of Warmest Quarter" from WorldClim and "Long-term maximum sea surface temperature" from Bio-ORACLE. For each realm, we standardized the long-term annual mean and maximum temperature across all the locations, by subtracting the mean and dividing by the standard deviation.


Jones, P.D. et al. Hemispheric and large-scale land-surface air temperature variations: An extensive revision and an update to 2010. J. Geophys. Res. Atmos. 117, D05127 (2012).
Harris, I., Jones, P.D., Osborn, T.J. & Lister, D.H. Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset. Int. J. Climatol. 34, 623-642 (2014).
Fick, S. E. & Hijmans, R. J. WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. Int. J. Climatol. 37, 4302-4315 (2017).
Tyberghein, L. et al. Bio-ORACLE: A global environmental dataset for marine species distribution modelling. Glob. Ecol. Biogeogr. 21, 272-281 (2012).
Assis, J. et al. Bio-ORACLE v2.0: Extending marine data layers for bioclimatic modelling. Glob. Ecol. Biogeogr. 27, 277-284 (2018).



## R Analysis Files 

These scripts produce plots summarising the data, prepare the data for analysis, fit the models for each metric, produce plots summarising the results, and illustrating the sensitivity analyses.

Script (01_) imports the data, produces maps and summary plots for the biodiversity and temperature trends, as well as plot summarising time series duration, number of years sampled, etc. (e.g. Figs. 2, S1 and S2).
Script (02_) fits the meta-analytical Bayesian models for each biodiversity metric and plots diagnostics of model fit.
Script (03_) contains code to produce figures summarising the main results (e.g. Figs 3, 4 and S3). 
script (04_) includes code to run the sensitivity analysis regarding the different baseline climate variables and produce figures (e.g. Fig. S4).
Script (05_) focuses on the sensitivity analysis regarding subsampling marine data (e.g. Fig. S5).
script (06_) produces figures for the remaining sensitivity analyses (e.g. Figs. S6 and S7).
script (07_) contains code to extract and plot the posterior densities of random parameters of interest per taxonomic group (e.g. Figs S8 and S9). 

Please note that some of the code in this repository was written to run on a HPC cluster.

* 01_data_import_summary.R
* 02_fit_models.R 
* 03_summarise_main_results.R
* 04_sensitivity_baseline_variables.R
* 05_sensitivity_marine_sampling.R
* 06_time_series_sensitivity.R
* 07_posterior_density_plots.R

This Rdata object contains all of the biodiversity and temperature trends estimates, as well as the baseline temperature values, for each time series. It further contains all the relevant metadata information.
* lm_slopes_meta_temperature.Rdata


* **model_fits_output**

This folder contains Rdata objects of all of the four biodiversity models.

    + NewRichTer_RS.Rdata
    + NewRichMar_RS.Rdata
    + NewAbundTer_RS.Rdata
    + NewAbundMar_RS.Rdata
    + NewGainsTer_RS.Rdata
    + NewGainsMar_RS.Rdata
    + NewLossTer_RS.Rdata
    + NewLossMar_RS.Rdata


* **sensitivity_analyses**

This folder contains four csv files with the results of running sensitivity analysis randomly subsetting marine data. It also contains an Rdata file with all the model estimates for sensitivity analysis for the different baseline climate variables.

    + RichMar_coeffsample.csv
    + AbundMar_coeffsample.csv
    + GainsMar_coeffsample.csv
    + LossesMar_coeffsample.csv
    
    + ALLcoeffs_baselines.Rdata

    

## Requirements
Program R version used 3.5.1

Necessary Packages: bayesplot, brms, broom, dplyr, ggplot2, ggridges, ggthemes, maps, maptools, purrr, readr, scales, tibble, tidyr, vegan, viridis



