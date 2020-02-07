##======================================================================
##	The biodiversity data used here originates from pre-processing steps
##	described in detail elsewhere.


# Primary Code Authors: Laura Antao
# Email: laura.antao@helsinki.fi
# DISCLAIMER: code is under development
# 31/01/2020


## Script (01_) imports the data, produces maps and summary plots for the biodiversity and temperature trends,
## as well as plot summarising time series duration, number of years sampled, etc. (e.g. Figs. 2, S1 and S2).
## Script (02_) fits the meta-analytical Bayesian models for each biodiversity metric and plots diagnostics of model fit.


# Input is object with biodiversity change estimates, temperature change estimates, and metadata for the time series

###the code assumes the libraries and objects from script 01_ are loaded in the R environment


##======================================================================

library(tidyverse)
library(brms)
library(ggpubr)
library(ggthemes)
library(viridis)
library(ggExtra)
library(cowplot)

# setwd("~/Temp_Biodiv_Change/R")

## Load data
#load("lm_slopes_meta_temperature3.Rdata")



##overall model: delta_Biodiversity ~ 0 + delta_Temperature * Baseline_Climate + (0+ delta_Temperature | Taxa) + (1 | Taxa/Study_ID)
##                                                    ||
##                                                    ||

##                                              baselines are

## (1) Long-term annual temperature (main results)
## (2) Long-term maximum temperature (e.g. Warmest-quarter)
## (3) Temperature in Year=1 from temporal database
## (4) Latitude

##the baseline temperature values were standardised as (x-mean/sd) per realm

##note the models can take several days to run

##the code below is for the main results using long-term annual temperature as baseline climate
##code for running models with the other variables can be easilly adapted
##(see also script "04_sensitivity_baseline_variables" for the different model formulas used)



###create formula for brms
metafull_formula <- bf(slope | se(std.error, sigma = TRUE) ~ 0 + TempGAMCoef * new_sTempYear + 
                         (0+ TempGAMCoef| taxa_mod1) + (1 | taxa_mod1/ STUDY_ID))



############################################################################
############################################################################

#####1. SPECIES RICHNESS #####
##create data subsets to run model for each realm separately

##Marine data
dataS_MarTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "logS_lm" & REALM == "Marine")


unique(dataS_MarTemp$model_id)
unique(dataS_MarTemp$REALM)
n_distinct(dataS_MarTemp$rarefyID)


##run model
NewRichMar_RS <- brm(formula= metafull_formula,
                     data= dataS_MarTemp,
                     control = list(adapt_delta = 0.99),
                     iter = 8000)

##the remaining arguments were used as default; models use non-informative flat priors

summary(NewRichMar_RS)
plot(NewRichMar_RS)

##save model output object
save(NewRichMar_RS, file="NewRichMar_RS.Rdata")


#####
##Terrestrial data
dataS_TerTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "logS_lm" & REALM == "Terrestrial") 


unique(dataS_TerTemp$model_id)
unique(dataS_TerTemp$REALM)
n_distinct(dataS_TerTemp$rarefyID)


##run model
NewRichTer_RS <- brm(formula= metafull_formula,
                     data= dataS_TerTemp,
                     control = list(adapt_delta = 0.99),
                     iter = 8000)

summary(NewRichTer_RS)
plot(NewRichTer_RS)

##save model output object
save(NewRichTer_RS, file="NewRichTer_RS.Rdata")




##======================================================================

#####2. ABUNDANCE #####

##Marine data
dataN_MarTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "logN_lm" & REALM == "Marine")

unique(dataN_MarTemp$model_id)
unique(dataN_MarTemp$REALM)
n_distinct(dataN_MarTemp$rarefyID)

NewAbundMar_RS <- brm(formula= metafull_formula,
                      data= dataN_MarTemp,
                      control = list(adapt_delta = 0.99),
                      iter = 8000)

summary(NewAbundMar_RS)
plot(NewAbundMar_RS)

save(NewAbundMar_RS, file="NewAbundMar_RS.Rdata")


#####
##Terrestrial data
dataN_TerTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "logN_lm" & REALM == "Terrestrial")

unique(dataN_TerTemp$model_id)
unique(dataN_TerTemp$REALM)
n_distinct(dataN_TerTemp$rarefyID)


NewAbundTer_RS <- brm(formula= metafull_formula,
                      data= dataN_TerTemp,
                      control = list(adapt_delta = 0.99),
                      iter = 8000)

summary(NewAbundTer_RS)
plot(NewAbundTer_RS)

save(NewAbundTer_RS, file="NewAbundTer_RS.Rdata")




##======================================================================

#####3. GAINS #####


##Marine data
dataG_MarTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "Gains_lm" & REALM == "Marine")

unique(dataG_MarTemp$model_id)
unique(dataG_MarTemp$REALM)
n_distinct(dataG_MarTemp$rarefyID)


NewGainsMar_RS <- brm(formula= metafull_formula,
                      data= dataG_MarTemp,
                      control = list(adapt_delta = 0.99),
                      iter = 8000)

summary(NewGainsMar_RS)
plot(NewGainsMar_RS)

save(NewGainsMar_RS, file="NewGainsMar_RS.Rdata")


#####
##Terrestrial data
dataG_TerTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "Gains_lm" & REALM == "Terrestrial")

unique(dataG_TerTemp$model_id)
unique(dataG_TerTemp$REALM)
n_distinct(dataG_TerTemp$rarefyID)


NewGainsTer_RS <- brm(formula= metafull_formula,
                      data= dataG_TerTemp,
                      control = list(adapt_delta = 0.99),
                      iter = 8000)

summary(NewGainsTer_RS)
plot(NewGainsTer_RS)

save(NewGainsTer_RS, file="NewGainsTer_RS.Rdata")




##======================================================================

#####4. LOSSES ####

##Marine data
dataL_MarTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "Losses_lm" & REALM == "Marine")

unique(dataL_MarTemp$model_id)
unique(dataL_MarTemp$REALM)
n_distinct(dataL_MarTemp$rarefyID)

NewLossMar_RS <- brm(formula= metafull_formula,
                     data= dataL_MarTemp,
                     control = list(adapt_delta = 0.99),
                     iter = 8000)

summary(NewLossMar_RS)
plot(NewLossMar_RS)

save(NewLossMar_RS, file="NewLossMar_RS.Rdata")


####
##Terrestrial data
dataL_TerTemp <- lm_slopes_meta_temperature3 %>%
  filter(model_id== "Losses_lm" & REALM == "Terrestrial")

unique(dataL_TerTemp$model_id)
unique(dataL_TerTemp$REALM)
n_distinct(dataL_TerTemp$rarefyID)


NewLossTer_RS <- brm(formula= metafull_formula,
                     data= dataL_TerTemp,
                     control = list(adapt_delta = 0.99),
                     iter = 8000)

summary(NewLossTer_RS)
plot(NewLossTer_RS)

save(NewLossTer_RS, file="NewLossTer_RS.Rdata")




##end models run##




##======================================================================

######to produce summary model tables with both fixed and random effects

##load model results objects
load("model_fits_output/NewRichTer_RS.Rdata")
load("model_fits_output/NewRichMar_RS.Rdata")
load("model_fits_output/NewAbundTer_RS.Rdata")
load("model_fits_output/NewAbundMar_RS.Rdata")
load("model_fits_output/NewGainsTer_RS.Rdata")
load("model_fits_output/NewGainsMar_RS.Rdata")
load("model_fits_output/NewLossTer_RS.Rdata")
load("model_fits_output/NewLossMar_RS.Rdata")


##for Richness
##marine
RichMar_postSummary<-data.frame(posterior_summary(NewRichMar_RS, robust = FALSE) %>% 
                                  round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Species Richness",
         Realm= "Marine") %>%
  select(Model, Realm, Term, everything()) %>%  #reorder cols
  filter(Term!= "lp__")


##terrestrial
RichTer_postSummary<-data.frame(posterior_summary(NewRichTer_RS, robust = FALSE) %>% 
                                  round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Species Richness",
         Realm= "Terrestrial") %>%
  select(Model, Realm, Term, everything()) %>%
  filter(Term!= "lp__")


##combine both
RichnessAll_postSummary<- rbind(RichMar_postSummary, RichTer_postSummary)

names(RichnessAll_postSummary)[c(1, 5:7)] <- c("Biodiversity Metric", "Std.error", "Lower 95% CI", "Upper 95% CI")


#write.csv(RichnessAll_postSummary,"RichnessAll_postSummary.csv", row.names = F)


########
##for Abundance
##marine
AbundMar_postSummary<-data.frame(posterior_summary(NewAbundMar_RS, robust = FALSE) %>% 
                                  round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Abundance",
         Realm= "Marine") %>%
  select(Model, Realm, Term, everything()) %>%  #reorder cols
  filter(Term!= "lp__")


##terrestrial
AbundTer_postSummary<-data.frame(posterior_summary(NewAbundTer_RS, robust = FALSE) %>% 
                                  round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Abundance",
         Realm= "Terrestrial") %>%
  select(Model, Realm, Term, everything()) %>%
  filter(Term!= "lp__")


##combine both
AbundanceAll_postSummary<- rbind(AbundMar_postSummary, AbundTer_postSummary)

names(AbundanceAll_postSummary)[c(1, 5:7)] <- c("Biodiversity Metric", "Std.error", "Lower 95% CI", "Upper 95% CI")


#write.csv(AbundanceAll_postSummary,"AbundanceAll_postSummary.csv", row.names = F)



########
##for Gains
##marine
GainMar_postSummary<-data.frame(posterior_summary(NewGainsMar_RS, robust = FALSE) %>% 
                                   round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Gains",
         Realm= "Marine") %>%
  select(Model, Realm, Term, everything()) %>%  #reorder cols
  filter(Term!= "lp__")


##terrestrial
GainTer_postSummary<-data.frame(posterior_summary(NewGainsTer_RS, robust = FALSE) %>% 
                                   round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Gains",
         Realm= "Terrestrial") %>%
  select(Model, Realm, Term, everything()) %>%
  filter(Term!= "lp__")


##combine both
GainAll_postSummary<- rbind(GainMar_postSummary, GainTer_postSummary)

names(GainAll_postSummary)[c(1, 5:7)] <- c("Biodiversity Metric", "Std.error", "Lower 95% CI", "Upper 95% CI")


#write.csv(GainAll_postSummary,"GainAll_postSummary.csv", row.names = F)




########
##for Losses
##marine
LossMar_postSummary<-data.frame(posterior_summary(NewLossMar_RS, robust = FALSE) %>% 
                                  round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Losses",
         Realm= "Marine") %>%
  select(Model, Realm, Term, everything()) %>%  #reorder cols
  filter(Term!= "lp__")


##terrestrial
LossTer_postSummary<-data.frame(posterior_summary(NewLossTer_RS, robust = FALSE) %>% 
                                  round(3)) %>%
  mutate(Term = make.names(rownames(.), unique=TRUE),
         Model = "Losses",
         Realm= "Terrestrial") %>%
  select(Model, Realm, Term, everything()) %>%
  filter(Term!= "lp__")


##combine both
LossAll_postSummary<- rbind(LossMar_postSummary, LossTer_postSummary)

names(LossAll_postSummary)[c(1, 5:7)] <- c("Biodiversity Metric", "Std.error", "Lower 95% CI", "Upper 95% CI")


#write.csv(LossAll_postSummary,"LossAll_postSummary.csv", row.names = F)



