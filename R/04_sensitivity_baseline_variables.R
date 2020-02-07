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


# Input is Rdata file containing all the model estimates for sensitivity analysis for the different baseline climate variables

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



## the code shows how the summaries from the different models were combined to procude the plot in Fig. S4
## we provide an Rdata file with the model estimates, rather than the model objects


###1. LONG-TERM MEAN TEMPERATURE (static databases) ====main results

# metafull_formula <- bf(slope | se(std.error, sigma = TRUE) ~ 0 + TempGAMCoef * new_sTempYear + 
#                          (0+ TempGAMCoef| taxa_mod1) + (1 | taxa_mod1/ STUDY_ID))

# load("model_fits_output/NewRichTer_RS.Rdata")
# load("model_fits_output/NewRichMar_RS.Rdata")
# load("model_fits_output/NewAbundTer_RS.Rdata")
# load("model_fits_output/NewAbundMar_RS.Rdata")
# load("model_fits_output/NewGainsTer_RS.Rdata")
# load("model_fits_output/NewGainsMar_RS.Rdata")
# load("model_fits_output/NewLossTer_RS.Rdata")
# load("model_fits_output/NewLossMar_RS.Rdata") ##named "AnnualTemp"



###2. LONG-TERM MAXIMUM TEMPERATURE (static databases)
# metafull_formulaMAX <- bf(slope | se(std.error, sigma = TRUE) ~ 0 + TempGAMCoef * MAX_sTempYear + 
#                             (0+ TempGAMCoef| taxa_mod1) + (1 | taxa_mod1/ STUDY_ID))

# load("~/NewRichTer_RS_max.Rdata")
# load("~/NewRichMar_RS_max.Rdata")
# load("~/NewAbundTer_RS_max.Rdata")
# load("~/NewAbundMar_RS_max.Rdata")
# load("~/NewGainsTer_RS_max.Rdata")
# load("~/NewGainsMar_RS_max.Rdata")
# load("~/NewLossTer_RS_max.Rdata")
# load("~/NewLossMar_RS_max.Rdata")   ##named "MaxTemp"


###3. YEAR=1 ANNUAL TEMPERATURE (temporal data)
# metafull_formulaY1 <- bf(slope | se(std.error, sigma = TRUE) ~ 0 + TempGAMCoef * sTemp_Year1 + 
#                            (0+ TempGAMCoef| taxa_mod1) + (1 | taxa_mod1/ STUDY_ID))

# load("~/NewRichTer_RS_Y1.Rdata")
# load("~/NewRichMar_RS_Y1.Rdata")
# load("~/NewAbundTer_RS_Y1.Rdata")
# load("~/NewAbundMar_RS_Y1.Rdata")
# load("~/NewGainsTer_RS_Y1.Rdata")
# load("~/NewGainsMar_RS_Y1.Rdata")
# load("~/NewLossTer_RS_Y1.Rdata")
# load("~/NewLossMar_RS_Y1.Rdata")   ##named "Year1Temp"



###4. LATITUDE as baseline
# metafull_formulaLAT <- bf(slope | se(std.error, sigma = TRUE) ~ 0 + TempGAMCoef * abs(rarefyID_y) +
# (0+ TempGAMCoef| taxa_mod1) + (1 | taxa_mod1/ STUDY_ID))

# load("~/NewRichTer_RS_LAT.Rdata")
# load("~/NewRichMar_RS_LAT.Rdata")
# load("~/NewAbundTer_RS_LAT.Rdata")
# load("~/NewAbundMar_RS_LAT.Rdata")
# load("~/NewGainsTer_RS_LAT.Rdata")
# load("~/NewGainsMar_RS_LAT.Rdata")
# load("~/NewLossTer_RS_LAT.Rdata")
# load("~/NewLossMar_RS_LAT.Rdata")    ##named "Latitude"



##======================================================================

##import the Rdata file from the folder
load("sensitivity_analyses/ALLcoeffs_baselines.Rdata")

###the code below produces the four panels in Fig. S4

##for richness #####

ALLcoeffs_SR <- rbind(data.frame(rbind(fixef(NewRichMar_RS), fixef(NewRichTer_RS)),
                                baseline =rep("AnnualTemp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####annual temp
                     
                     data.frame(rbind(fixef(NewRichMar_RS_max), fixef(NewRichTer_RS_max)),
                                baseline =rep("MaxTemp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####max temp
                     
                     data.frame(rbind(fixef(NewRichMar_RS_Y1), fixef(NewRichTer_RS_Y1)),
                                baseline =rep("Year1Temp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####YEAR=1 temp
                     
                     data.frame(rbind(fixef(NewRichMar_RS_LAT), fixef(NewRichTer_RS_LAT)),
                                baseline =rep("Latitude",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3))) %>%   ###Latitude
  
  mutate(Term= rep(c("Temp change", "Climatology", "Interaction"), 8))

ALLcoeffs_SR$baseline <- factor(ALLcoeffs_SR$baseline, levels=c("AnnualTemp", "MaxTemp", "Year1Temp", "Latitude"))


###different baseline variable comparisons
pal <- viridisLite::viridis(4)

(SR_compare<- ggplot(data = ALLcoeffs_SR, 
                     aes(x = Term, y = Estimate, group = Term, label = Term, color=baseline)) +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), position=position_jitter(), size=.9, alpha=.8) +
    scale_x_discrete(limits=c("Interaction", "Climatology", "Temp change"),  ##order of plot for x axis 
                     labels = c("Interaction" = "Interaction",
                                "Climatology" = "Baseline Climate\n/Latitude",
                                "Temp change" = "Temperature change")) +   ##rename parameters 
    coord_flip() +
    labs(y="Estimated Coefficient", x= "", title = "Richness", color="") +
    geom_hline(yintercept=0, lty=2) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = c(.79,.2),
          strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(face = "bold", size=10),
          plot.margin = unit(c(0.1, 0.5, 0.1, -0.5), "lines")) +
    facet_grid(~ Realm) +
    scale_color_manual(values = c(pal[-4], "gold2")))

#ggsave("SR_compare.tiff")



##for abundance #####

ALLcoeffs_N <- rbind(data.frame(rbind(fixef(NewAbundMar_RS), fixef(NewAbundTer_RS)),
                                baseline =rep("AnnualTemp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####annual temp
                     
                     data.frame(rbind(fixef(NewAbundMar_RS_max), fixef(NewAbundTer_RS_max)),
                                baseline =rep("MaxTemp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####max temp
                     
                     data.frame(rbind(fixef(NewAbundMar_RS_Y1), fixef(NewAbundTer_RS_Y1)),
                                baseline =rep("Year1Temp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####YEAR=1 temp
                     
                     data.frame(rbind(fixef(NewAbundMar_RS_LAT), fixef(NewAbundTer_RS_LAT)),
                                baseline =rep("Latitude",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3))) %>%   ###Latitude
  
  mutate(Term= rep(c("Temp change", "Climatology", "Interaction"), 8))

ALLcoeffs_N$baseline <- factor(ALLcoeffs_N$baseline, levels=c("AnnualTemp", "MaxTemp", "Year1Temp", "Latitude"))


(Abund_compare<- ggplot(data = ALLcoeffs_N, 
                        aes(x = Term, y = Estimate, group = Term, label = Term, color=baseline)) +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), position=position_jitter(), size=.9, alpha=.8) +
    scale_x_discrete(limits=c("Interaction", "Climatology", "Temp change")) +   ##order of plot for x axis 
    coord_flip() +
    labs(y="Estimated Coefficient", x= "", title = "Abundance", color="") +
    geom_hline(yintercept=0, lty=2) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(face = "bold", size=10),
          plot.margin = unit(c(0.1, 0.5, 0.1, -0.5), "lines")) +
    facet_grid(~ Realm) +
    scale_color_manual(values = c(pal[-4], "gold2")))

#ggsave("Abund_compare.tiff")



## for gains #####

ALLcoeffs_G <- rbind(data.frame(rbind(fixef(NewGainsMar_RS), fixef(NewGainsTer_RS)),
                                baseline =rep("AnnualTemp",6),
                               Realm =rep(c("Marine","Terrestrial"),each=3)),  ####annual temp
                    
                    data.frame(rbind(fixef(NewGainsMar_RS_max), fixef(NewGainsTer_RS_max)),
                               baseline =rep("MaxTemp",6),
                               Realm =rep(c("Marine","Terrestrial"),each=3)),  ####max temp
                    
                    data.frame(rbind(fixef(NewGainsMar_RS_Y1), fixef(NewGainsTer_RS_Y1)),
                               baseline =rep("Year1Temp",6),
                               Realm =rep(c("Marine","Terrestrial"),each=3)),  ####YEAR=1 temp
                    
                    data.frame(rbind(fixef(NewGainsMar_RS_LAT), fixef(NewGainsTer_RS_LAT)),
                               baseline =rep("Latitude",6),
                               Realm =rep(c("Marine","Terrestrial"),each=3))) %>%   ###Latitude
  
  mutate(Term= rep(c("Temp change", "Climatology", "Interaction"), 8))

ALLcoeffs_G$baseline <- factor(ALLcoeffs_G$baseline, levels=c("oldERR", "AnnualTemp", "MaxTemp", "Year1Temp", "Latitude"))


(Gains_compare<- ggplot(data = ALLcoeffs_G,
                        aes(x = Term, y = Estimate, group = Term, label = Term, color=baseline)) +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), position=position_jitter(), size=.9, alpha=.8) +
    scale_x_discrete(limits=c("Interaction", "Climatology", "Temp change")) +   ##order of plot for x axis 
    coord_flip() +
    labs(y="Estimated Coefficient", x= "", title = "Gains", color="") +
    geom_hline(yintercept=0, lty=2) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(face = "bold", size=10),
          plot.margin = unit(c(0.1, 0.5, 0.1, -0.5), "lines")) +
    facet_grid(~ Realm) +
    scale_color_manual(values = c(pal[-4], "gold2")))

#ggsave("Gains_compare.tiff")



## for losses #####

ALLcoeffs_L <- rbind(data.frame(rbind(fixef(NewLossMar_RS), fixef(NewLossTer_RS)),
                                baseline =rep("AnnualTemp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####annual temp
                     
                     data.frame(rbind(fixef(NewLossMar_RS_max), fixef(NewLossTer_RS_max)),
                                baseline =rep("MaxTemp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####max temp
                     
                     data.frame(rbind(fixef(NewLossMar_RS_Y1), fixef(NewLossTer_RS_Y1)),
                                baseline =rep("Year1Temp",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3)),  ####YEAR=1 temp
                     
                     data.frame(rbind(fixef(NewLossMar_RS_LAT), fixef(NewLossTer_RS_LAT)),
                                baseline =rep("Latitude",6),
                                Realm =rep(c("Marine","Terrestrial"),each=3))) %>%   ###Latitude
  
  mutate(Term= rep(c("Temp change", "Climatology", "Interaction"), 8))

ALLcoeffs_L$baseline <- factor(ALLcoeffs_L$baseline, levels=c("AnnualTemp", "MaxTemp", "Year1Temp", "Latitude"))


(Losses_compare<- ggplot(data = ALLcoeffs_L,
                         aes(x = Term, y = Estimate, group = Term, label = Term, color=baseline)) +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), position=position_jitter(), size=.9, alpha=.8) +
    scale_x_discrete(limits=c("Interaction", "Climatology", "Temp change")) +   ##order of plot for x axis 
    coord_flip() +
    labs(y="Estimated Coefficient", x= "", title = "Losses", color="") +
    geom_hline(yintercept=0, lty=2) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(face = "bold", size=10),
          plot.margin = unit(c(0.1, 0.5, 0.1, -0.5), "lines")) +
    facet_grid(~ Realm) +
    scale_color_manual(values = c(pal[-4], "gold2")))

#ggsave("Losses_compare.tiff")

##then combine the four plots




# save(ALLcoeffs_SR, ALLcoeffs_N, ALLcoeffs_G, ALLcoeffs_L, file="ALLcoeffs_baselines.Rdata")

##then import this object from the GitHub repo






