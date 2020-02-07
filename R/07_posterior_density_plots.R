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
## Script (03_) contains code to produce figures summarising the main results (e.g. Figs 3, 4 and S3). 
## script (04_) includes code to run the sensitivity analysis regarding the different baseline climate variables and produce figures (e.g. Fig. S4).
## Script (05_) focuses on the sensitivity analysis regarding subsampling marine data (e.g. Fig. S5).
## script (06_) produces figures for the remaining sensitivity analyses (e.g. Figs. S6 and S7).
## script (07_) contains code to extract and plot the posterior densities of random parameters of interest per taxonomic group (e.g. Figs S8 and S9).


# Input are the model objects saved (see script 02_)

###the code assumes the libraries and objects from previous scripts are loaded in the R environment


##======================================================================

library(tidyverse)
library(brms)
library(ggpubr)
library(ggthemes)
library(viridis)
library(ggExtra)
library(cowplot)
library(ggridges)


# setwd("~/Temp_Biodiv_Change/R")


load("model_fits_output/NewRichTer_RS.Rdata")
load("model_fits_output/NewRichMar_RS.Rdata")
load("model_fits_output/NewAbundTer_RS.Rdata")
load("model_fits_output/NewAbundMar_RS.Rdata")
load("model_fits_output/NewGainsTer_RS.Rdata")
load("model_fits_output/NewGainsMar_RS.Rdata")
load("model_fits_output/NewLossTer_RS.Rdata")
load("model_fits_output/NewLossMar_RS.Rdata")




## this code extracts samples for the random slopes to produce density plots for the differente taxonomic groups
## and produces Figs. S8 and S9 in the manuscript

###for Richness Marine####

#get posterior draws for overall trend
RichMar_global_posterior <- tibble(
  S_global = posterior_samples(NewRichMar_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>%
    unlist() %>%
    as.numeric())

###get the levels for taxa
taxa_levels <- NewRichMar_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)


##get posterior distribution of the taxa slopes
richMar_sample_posterior <- taxa_levels %>%
  mutate(S_Mar_postSamp = purrr::map(data, ~posterior_samples(NewRichMar_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>%
                                       unlist() %>%
                                       as.numeric()))


richMar_sample_posterior <- richMar_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(S_global_slope = rep(RichMar_global_posterior$S_global, times = n_distinct(taxa_mod1)))



##to have nicer lables
richMar_sample_posterior$taxa2 <- factor(richMar_sample_posterior$taxa_mod1,
                                         levels = c("All", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"),
                                         labels = c("Multiple taxa", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"))


##plotting slopes
(richMar_slopes<- ggplot() +
    geom_density_ridges(data = richMar_sample_posterior,
                        aes(x = S_Mar_postSamp + S_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6)) +
    geom_vline(aes(xintercept = fixef(NewRichMar_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = '',
         x = 'Taxon-level slopes') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    xlim(-NA,1) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Richness", Term=="Temp change", Realm=="Marine"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=richMar_slopes, "richMar_slopes.tiff", width = 14, height = 12, units = c("cm"))



###for Richness Terrestrial####

#get posterior draws for overall trend
RichTer_global_posteriorTer <- tibble(
  S_global = posterior_samples(NewRichTer_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
    unlist() %>% 
    as.numeric())


######get the levels for taxa
taxa_levelsT <- NewRichTer_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)


##get posterior distribution of the taxa slopes
richTer_sample_posterior <- taxa_levelsT %>%
  mutate(S_Ter_postSamp = purrr::map(data, ~posterior_samples(NewRichTer_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
                                       unlist() %>%
                                       as.numeric()))

richTer_sample_posterior <- richTer_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(S_global_slope = rep(RichTer_global_posteriorTer$S_global, times = n_distinct(taxa_mod1)))


##to have nice lables
richTer_sample_posterior$taxa2 <- factor(richTer_sample_posterior$taxa_mod1,
                                         levels = c("All", "Terrestrial plants", "Birds", "Mammals", "Terrestrial invertebrates", "Amphibians"),
                                         labels = c("Multiple taxa", "Plants", "Birds", "Mammals", "Invertebrates", "Amphibians"))



##plotting slopes
(richTer_slopes<- ggplot() +
    geom_density_ridges(data = richTer_sample_posterior,
                        aes(x = S_Ter_postSamp + S_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6, option = "inferno")) +
    geom_vline(aes(xintercept = fixef(NewRichTer_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = "",
         x = "Taxon-level slopes") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    xlim(-0.5,NA) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Richness", Term=="Temp change", Realm=="Terrestrial"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=richTer_slopes, "richTer_slopes.tiff", width = 14, height = 12, units = c("cm"))



##======================================================================

###for Abundance Marine####

#get posterior draws for overall trend
AbundMar_global_posterior <- tibble(
  S_global = posterior_samples(NewAbundMar_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
    unlist() %>% 
    as.numeric())


######get the levels for taxa
taxa_levelsA <- NewAbundMar_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)


##get posterior distribution of the taxa slopes
abundMar_sample_posterior <- taxa_levelsA %>%
  mutate(N_Mar_postSamp = purrr::map(data, ~posterior_samples(NewAbundMar_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
                                       unlist() %>% 
                                       as.numeric()))


abundMar_sample_posterior <- abundMar_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(N_global_slope = rep(AbundMar_global_posterior$S_global, times = n_distinct(taxa_mod1)))


##to have nice lables
abundMar_sample_posterior$taxa2 <- factor(abundMar_sample_posterior$taxa_mod1,
                                          levels = c("All", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"),
                                          labels = c("Multiple taxa", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"))


##plotting slopes
(abunMar_slopes<- ggplot() +
    geom_density_ridges(data = abundMar_sample_posterior,
                        aes(x = N_Mar_postSamp + N_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6)) +
    geom_vline(aes(xintercept = fixef(NewAbundMar_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = '',
         x = 'Taxon-level slopes') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    xlim(NA,4.5) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Abundance", Term=="Temp change", Realm=="Marine"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=abunMar_slopes, "abunMar_slopes.tiff", width = 14, height = 12, units = c("cm"))



###for Abundance Terrestrial####

#get posterior draws for overall trend
AbundTer_global_posterior <- tibble(
  S_global = posterior_samples(NewAbundTer_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
    unlist() %>% 
    as.numeric())

######get the levels for taxa
taxa_levelsAt <- NewAbundTer_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)


##get posterior distribution of the taxa slopes
abundTer_sample_posterior <- taxa_levelsAt %>%
  mutate(N_Ter_postSamp = purrr::map(data, ~posterior_samples(NewAbundTer_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
                                       unlist() %>% 
                                       as.numeric()))


abundTer_sample_posterior <- abundTer_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(N_global_slope = rep(AbundTer_global_posterior$S_global, times = n_distinct(taxa_mod1)))


##to have nice lables
abundTer_sample_posterior$taxa2 <- factor(abundTer_sample_posterior$taxa_mod1,
                                          levels = c("All", "Terrestrial plants", "Birds", "Mammals", "Terrestrial invertebrates", "Amphibians"),
                                          labels = c("Multiple taxa", "Plants", "Birds", "Mammals", "Invertebrates", "Amphibians"))


##plotting slopes
(abunTer_slopes<- ggplot() +
    geom_density_ridges(data = abundTer_sample_posterior,
                        aes(x = N_Ter_postSamp + N_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6, option = "inferno")) +
    geom_vline(aes(xintercept = fixef(NewAbundTer_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = '',
         x = 'Taxon-level slopes') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    xlim(-0.9,1) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Abundance", Term=="Temp change", Realm=="Terrestrial"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=abunTer_slopes, "abunTer_slopes.tiff", width = 14, height = 12, units = c("cm"))



##======================================================================

###for Gains Marine####

#get posterior draws for overall trend
GainMar_global_posterior <- tibble(
  S_global = posterior_samples(NewGainsMar_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
    unlist() %>% 
    as.numeric())

######get the levels for taxa
taxa_levelsG <- NewGainsMar_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)


##get posterior distribution of the taxa slopes
GainMar_sample_posterior <- taxa_levelsG %>%
  mutate(G_Mar_postSamp = purrr::map(data, ~posterior_samples(NewGainsMar_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>%
                                       unlist() %>% 
                                       as.numeric()))


GainMar_sample_posterior <- GainMar_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(N_global_slope = rep(GainMar_global_posterior$S_global, times = n_distinct(taxa_mod1)))



##to have nice lables
GainMar_sample_posterior$taxa2 <- factor(GainMar_sample_posterior$taxa_mod1,
                                         levels = c("All", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"),
                                         labels = c("Multiple taxa", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"))


##plotting slopes
(gainMar_slopes<- ggplot() +
    geom_density_ridges(data = GainMar_sample_posterior,
                        aes(x = G_Mar_postSamp + N_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6)) +
    geom_vline(aes(xintercept = fixef(NewGainsMar_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = '',
         x = 'Taxon-level slopes') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    xlim(NA,3.5) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Gains", Term=="Temp change", Realm=="Marine"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=gainMar_slopes, "gainMar_slopes.tiff", width = 14, height = 12, units = c("cm"))



###for Gains Terrestrial####

#get posterior draws for overall trend
GainTer_global_posterior <- tibble(
  S_global = posterior_samples(NewGainsTer_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>%
    unlist() %>% 
    as.numeric())


######get the levels for taxa
taxa_levelsGt <- NewGainsTer_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)


##get posterior distribution of the taxa slopes
GainTer_sample_posterior <- taxa_levelsGt %>%
  mutate(G_Ter_postSamp = purrr::map(data, ~posterior_samples(NewGainsTer_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
                                       unlist() %>% 
                                       as.numeric()))


GainTer_sample_posterior <- GainTer_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(N_global_slope = rep(GainTer_global_posterior$S_global, times = n_distinct(taxa_mod1)))



##to have nice lables
GainTer_sample_posterior$taxa2 <- factor(GainTer_sample_posterior$taxa_mod1,
                                         levels = c("All", "Terrestrial plants", "Birds", "Mammals", "Terrestrial invertebrates", "Amphibians"),
                                         labels = c("Multiple taxa", "Plants", "Birds", "Mammals", "Invertebrates", "Amphibians"))


##plotting slopes
(gainTer_slopes<- ggplot() +
    geom_density_ridges(data = GainTer_sample_posterior,
                        aes(x = G_Ter_postSamp + N_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6, option = "inferno")) +
    geom_vline(aes(xintercept = fixef(NewGainsTer_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = '',
         x = 'Taxon-level slopes') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    xlim(-3.5,3.5) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Gains", Term=="Temp change", Realm=="Terrestrial"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=gainTer_slopes, "gainTer_slopes.tiff", width = 14, height = 12, units = c("cm"))



##======================================================================

###for Losses Marine####

#get posterior draws for overall trend
LossMar_global_posterior <- tibble(
  S_global = posterior_samples(NewLossMar_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>%
    unlist() %>% 
    as.numeric())


######get the levels for taxa
taxa_levelsL <- NewLossMar_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)



##get posterior distribution of the taxa slopes
LossMar_sample_posterior <- taxa_levelsL %>%
  mutate(L_Mar_postSamp = purrr::map(data, ~posterior_samples(NewLossMar_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>%
                                       unlist() %>% 
                                       as.numeric()))


LossMar_sample_posterior <- LossMar_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(N_global_slope = rep(LossMar_global_posterior$S_global, times = n_distinct(taxa_mod1)))


##to have nice lables
LossMar_sample_posterior$taxa2 <- factor(LossMar_sample_posterior$taxa_mod1,
                                         levels = c("All", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"),
                                         labels = c("Multiple taxa", "Fish", "Benthos", "Birds", "Marine invertebrates", "Mammals"))


##plotting slopes
(lossMar_slopes<- ggplot() +
    geom_density_ridges(data = LossMar_sample_posterior,
                        aes(x = L_Mar_postSamp + N_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6)) +
    geom_vline(aes(xintercept = fixef(NewLossMar_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = '',
         x = 'Taxon-level slopes') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    #xlim(-1.5,1.2) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Losses", Term=="Temp change", Realm=="Marine"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=lossMar_slopes, "lossMar_slopes.tiff", width = 14, height = 12, units = c("cm"))



###for Losses Terrestrial####

#get posterior draws for overall trend
LossTer_global_posterior <- tibble(
  S_global = posterior_samples(NewLossTer_RS,
                               pars = 'b_TempGAMCoef',
                               exact_match = TRUE,
                               subset = floor(runif(n = 1000, min = 1, max = 2000))) %>%
    unlist() %>% 
    as.numeric())

######get the levels for taxa
taxa_levelsLt <- NewLossTer_RS$data %>% 
  as_tibble() %>% 
  distinct(taxa_mod1) %>% 
  mutate(level = make.names(taxa_mod1, unique=TRUE)) %>%  ###this removes the space when creating the levels
  nest(level)


##get posterior distribution of the taxa slopes
LossTer_sample_posterior <- taxa_levelsLt %>%
  mutate(L_Ter_postSamp = purrr::map(data, ~posterior_samples(NewLossTer_RS, 
                                                              pars = paste('r_taxa_mod1[', as.character(.x$level), ',TempGAMCoef]', sep=''),
                                                              exact = TRUE,
                                                              subset = floor(runif(n = 1000, min = 1, max = 2000))) %>% 
                                       unlist() %>% 
                                       as.numeric()))


LossTer_sample_posterior <- LossTer_sample_posterior %>% 
  select(-data) %>% 
  unnest() %>%
  mutate(N_global_slope = rep(LossTer_global_posterior$S_global, times = n_distinct(taxa_mod1)))



##to have nice lables
LossTer_sample_posterior$taxa2 <- factor(LossTer_sample_posterior$taxa_mod1,
                                         levels = c("All", "Terrestrial plants", "Birds", "Mammals", "Terrestrial invertebrates", "Amphibians"),
                                         labels = c("Multiple taxa", "Plants", "Birds", "Mammals", "Invertebrates", "Amphibians"))


##plotting slopes
(lossTer_slopes<- ggplot() +
    geom_density_ridges(data = LossTer_sample_posterior,
                        aes(x = L_Ter_postSamp + N_global_slope, y = taxa2,
                            fill = taxa2),
                        scale = 2, alpha = 0.6,
                        linetype = 0) +
    scale_fill_manual(name = 'Taxon group', values = viridis(6, option = "inferno")) +
    geom_vline(aes(xintercept = fixef(NewLossTer_RS)[1]), lwd=1.2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    labs(y = '',
         x = 'Taxon-level slopes') +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size= 12),
          axis.title = element_text(size= 15),
          axis.text.y =  element_text(size= 10),
          legend.key = element_blank(),
          legend.position = 'none',
          legend.direction = 'horizontal',
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    xlim(-4,4.5) +
    
    ##add grey shading for the CI of main effect
    geom_rect(data=newcoeffALL %>%
                filter(Metric=="Losses", Term=="Temp change", Realm=="Terrestrial"),
              aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.4, inherit.aes = FALSE))

ggsave(plot=lossTer_slopes, "lossTer_slopes.tiff", width = 14, height = 12, units = c("cm"))




