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


# Input is Rdata object with biodiversity change estimates, temperature change estimates, and metadata for the time series

###the code assumes the libraries and objects from previous scripts are loaded in the R environment


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



## this code produces plots to show the patterns in biodiversity responses as a funtion of
## the number of years sampled, temporal duration or start year of the time series
## and produces Figs. S6 and S7 in the manuscript

lm_slopes_meta_temperature3$model_id = factor(lm_slopes_meta_temperature3$model_id, levels=c("logS_lm","logN_lm","Gains_lm","Losses_lm"))

##to plot different trends
col_change = c("up" = "blue", "down" = "red", "neutral" = "grey")

###to rename the panels
labels <- c(logS_lm = " Richness", logN_lm = "Abundance", Gains_lm = "Gains", Losses_lm = "Losses")


##to produce Fig. S6 == variation in estimated slopes

###slopes ~ number of years sampled
(slopes_numYears <- lm_slopes_meta_temperature3 %>%
    ggplot() +
    facet_wrap(~model_id, scales="free", nrow = 1, labeller=labeller(model_id = labels)) +
    geom_point(aes(num_years, slope, colour=change_mu, size = duration), alpha = 0.2) +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(name="Time series\nchange", values = col_change) +
    scale_size_area(name="Time series\nduration", breaks = c(2,4,8,16,32,64,128)) +
    labs(x = "Number of years sampled", y = "Slope estimate") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size= 16),
          strip.text = element_text(size = 16),
          plot.subtitle = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          strip.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size = 3))))


###slopes ~ duration of time-series
(slopes_duration <- lm_slopes_meta_temperature3 %>%
    ggplot() +
    facet_wrap(~model_id, scales="free", nrow = 1, labeller=labeller(model_id = labels)) +
    geom_point(aes(duration, slope, colour=change_mu, size = duration), alpha = 0.2) +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(guide = F,
                       name="Time series\nchange", values = col_change) +
    scale_size_area(name = "Time series\nduration", breaks = c(2,4,8,16,32,64)) +
    labs(x = "Duration of Study", y = "Slope estimate") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size= 16),
          strip.text = element_text(size = 16),
          plot.subtitle = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          strip.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size = 3))))


###slopes ~ start years sampled
(slopes_startYear <- lm_slopes_meta_temperature3 %>%
    ggplot() +
    facet_wrap(~model_id, scales="free", nrow = 1, labeller=labeller(model_id = labels)) +
    geom_point(aes(startYear, slope, colour=change_mu, size = num_years), alpha = 0.2) +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(guide = FALSE, 
                       name="Time series\nchange", values = col_change) +
    scale_size_area(name = "Number of\nyears sampled", breaks = c(2,4,8,16,32,64)) +
    labs(x = "Start year", y = "Slope estimate") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size= 16),
          strip.text = element_text(size = 16),
          plot.subtitle = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          strip.background = element_blank()))


##combine the plots
plot_grid(slopes_numYears, slopes_duration, slopes_startYear, nrow = 3, align = "hv", labels = c("a)","b)","c)"))
ggsave("FigS6_slope_numyears_sensitivity.png", width = 400, height = 300, units = "mm")





##======================================================================

##to produce Fig. S7 == variation in estimated standard errors

###std.error ~ number of years sampled
(stderror_numYears <- lm_slopes_meta_temperature3 %>%
   ggplot() +
   facet_wrap(~model_id, scales="free", nrow = 1, labeller=labeller(model_id = labels)) +
   geom_point(aes(num_years, std.error, colour=change_mu, size = duration), alpha = 0.2) +
   geom_hline(yintercept = 0, lty = 2) +
   scale_color_manual(name="Time series\nchange", values = col_change) +
   scale_size_area(name="Time series\nduration", breaks = c(2,4,8,16,32,64,128)) +
   labs(x = "Number of years sampled", y = "Std error estimate") +
   theme_bw() +
   theme(axis.text = element_text(size = 16),
         axis.title = element_text(size= 16),
         strip.text = element_text(size = 16),
         plot.subtitle = element_text(size = 14, face = "bold"),
         legend.text = element_text(size = 14),
         legend.title = element_text(size = 16),
         strip.background = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha=1, size = 3))))


###std.error ~ duration of time-series
(stderror_duration <- lm_slopes_meta_temperature3 %>%
    ggplot() +
    facet_wrap(~model_id, scales="free", nrow = 1, labeller=labeller(model_id = labels)) +
    geom_point(aes(duration, std.error, colour=change_mu, size = duration), alpha = 0.2) +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(guide = F,
                       name="Time series\nchange", values = col_change) +
    scale_size_area(name = "Time series\nduration", breaks = c(2,4,8,16,32,64)) +
    labs(x = "Duration of Study", y = "Std error estimate") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size= 16),
          strip.text = element_text(size = 16),
          plot.subtitle = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          strip.background = element_blank()))



###std.error ~ start year
(stderror_startYear <- lm_slopes_meta_temperature3 %>%
    ggplot() +
    facet_wrap(~model_id, scales="free", nrow = 1, labeller=labeller(model_id = labels)) +
    geom_point(aes(startYear, std.error, colour=change_mu, size = num_years), alpha = 0.2) +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(guide = FALSE, 
                       name="Time series\nchange", values = col_change) +
    scale_size_area(name = "Number of\nyears sampled", breaks = c(2,4,8,16,32,64)) +
    labs(x = "Start year", y = "Std error estimate") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size= 16),
          strip.text = element_text(size = 16),
          plot.subtitle = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          strip.background = element_blank()))


##combine the plots
plot_grid(stderror_numYears, stderror_duration, stderror_startYear, nrow = 3, align = "hv", labels = c("a)","b)","c)"))
ggsave("FigS7_stderror_numyears_sensitivity.png", width = 400, height = 300, units = "mm")







