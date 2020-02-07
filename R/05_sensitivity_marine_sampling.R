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


# Input are csv files containing the results of running sensitivity analysis by randomly subsetting marine data

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


## this code subsamples the marine data to match the number of time series
## and the latitdinal span of the terrestrial data
## and then fits the meta-analytical models to 100 random subsamples



###create formula for brms
metafull_formula <- bf(slope | se(std.error, sigma = TRUE) ~ 0 + TempGAMCoef * new_sTempYear + 
                         (0+ TempGAMCoef| taxa_mod1) + (1 | taxa_mod1/ STUDY_ID))


########################
#####1. RICHNESS
dataS_Ter <- lm_slopes_meta_temperature %>%
  filter(model_id== "logS_lm" & REALM == "Terrestrial")

unique(dataS_Ter$model_id)
unique(dataS_Ter$REALM)
n_distinct(dataS_Ter$rarefyID)



##create list to save the parameters of each run
RichMar_sensitivity<- list()

for(u in 1:100){
  
  ####create the sample Marine dataset
  dataS_Mar <- lm_slopes_meta_temperature %>%
    filter(model_id== "logS_lm" & REALM == "Marine") %>%
    
    ##select latitudinal range only
    filter(rarefyID_y >= min(dataS_Ter$rarefyID_y) & rarefyID_y <= max(dataS_Ter$rarefyID_y)) %>%
    
    ##sample number of locations
    sample_n(., n_distinct(dataS_Ter$rarefyID))
  
  
  ##run the model for each subset
  modbrmssubset <- brm(formula= metafull_formula,
                       data= dataS_Mar,
                       control = list(adapt_delta = 0.99),
                       iter = 8000)
  
  
  ##save estimates into list element
  RichMar_sensitivity[[u]]<- data.frame(summary(modbrmssubset)[14]) %>%
    mutate(Metric = rep("Richness", 3),
           term1 = c("Temp change", "Average Temp", "Interaction"),
           run= u)
  
  print(u)
  
}

##combine all runs into a dataframe
RichMar_coeffsample <- do.call(rbind, RichMar_sensitivity)


write.csv(RichMar_coeffsample, "RichMar_coeffsample.csv")



########################################
######2. ABUNDANCE
##create list to save the parameters
AbundMar_sensitivity<- list()

for(u in 1:100){
  
  ####create the sample Marine dataset
  dataAbd_Mar <- lm_slopes_meta_temperature %>%
    filter(model_id== "logN_lm" & REALM == "Marine") %>%
    
    ##select latitudinal range only
    filter(rarefyID_y >= min(dataS_Ter$rarefyID_y) & rarefyID_y <= max(dataS_Ter$rarefyID_y)) %>%
    
    sample_n(., n_distinct(dataS_Ter$rarefyID))
  
  
  ##run the model for each data
  modbrmssubsetAbd <- brm(formula= metafull_formula,
                          data= dataAbd_Mar,
                          control = list(adapt_delta = 0.99),
                          iter = 8000)
  
  
  ##saves estimates into list element
  AbundMar_sensitivity[[u]]<- data.frame(summary(modbrmssubsetAbd)[14]) %>%
    mutate(Metric = rep("Abundance", 3),
           term1 = c("Temp change", "Average Temp", "Interaction"),
           run= u)
  
  print(u)
  
}

##combine all runs into a dataframe
AbundMar_coeffsample<- do.call(rbind, AbundMar_sensitivity) ##combine all runs into a dataframe

write.csv(AbundMar_coeffsample, "AbundMar_coeffsample.csv")





########################################
######3. GAINS
##create list to save the parameters
GainsMar_sensitivity<- list()

for(u in 1:100){
  
  ####create the sample Marine dataset
  dataG_Mar <- lm_slopes_meta_temperature %>%
    filter(model_id== "Gains_lm" & REALM == "Marine") %>%
    
    ##select latitudinal range only
    filter(rarefyID_y >= min(dataS_Ter$rarefyID_y) & rarefyID_y <= max(dataS_Ter$rarefyID_y)) %>%
    
    sample_n(., n_distinct(dataS_Ter$rarefyID))
  
  ##run the model for each data
  modbrmssubsetG <- brm(formula= metafull_formula,
                        data= dataG_Mar,
                        control = list(adapt_delta = 0.99),
                        iter = 8000)
  
  
  ##saves estimates into list element
  GainsMar_sensitivity[[u]]<- data.frame(summary(modbrmssubsetG)[14]) %>%
    mutate(Metric = rep("Gains",3),
           term1= c("Temp change", "Average Temp", "Interaction"),
           run= u)
  
}

##combine all runs into a dataframe
GainsMar_coeffsample<- do.call(rbind, GainsMar_sensitivity) ##combine all runs into a dataframe

write.csv(GainsMar_coeffsample, "GainsMar_coeffsample.csv")




########################################
######4. LOSSES
##create list to save the parameters
LossesMar_sensitivity<- list()

for(u in 1:100){
  
  ####create the sample Marine dataset
  dataL_Mar <- lm_slopes_meta_temperature %>%
    filter(model_id== "Losses_lm" & REALM == "Marine") %>%
    
    ##select latitudinal range only
    filter(rarefyID_y >= min(dataS_Ter$rarefyID_y) & rarefyID_y <= max(dataS_Ter$rarefyID_y)) %>%
    
    sample_n(., n_distinct(dataS_Ter$rarefyID))
  
  ##run the model for each data
  modbrmssubsetL <- brm(formula= metafull_formula,
                        data= dataL_Mar,
                        control = list(adapt_delta = 0.99),
                        iter = 8000)
  
  
  ##saves estimates into list element
  LossesMar_sensitivity[[u]]<- data.frame(summary(modbrmssubsetL)[14]) %>%
    mutate(Metric =rep("Losses", 3),
           term1= c("Temp change", "Average Temp", "Interaction"),
           run= u)
  
  print(u)
  
}

##combine all runs into a dataframe
LossesMar_coeffsample<- do.call(rbind, LossesMar_sensitivity) ##combine all runs into a dataframe


write.csv(LossesMar_coeffsample, "LossesMar_coeffsample.csv")







##======================================================================


####for producing the plots
##import the csv files from the folder

RichMar_coeffsample<- read.csv("sensitivity_analyses/RichMar_coeffsample.csv")
AbundMar_coeffsample<- read.csv("sensitivity_analyses/AbundMar_coeffsample.csv")
GainsMar_coeffsample<- read.csv("sensitivity_analyses/GainsMar_coeffsample.csv")
LossesMar_coeffsample<- read.csv("sensitivity_analyses/LossesMar_coeffsample.csv")

##renaming variable
RichMar_coeffsample$Term<- plyr::revalue(RichMar_coeffsample$Term, c("Average Temp"="Baseline Climate"))
AbundMar_coeffsample$Term<- plyr::revalue(AbundMar_coeffsample$Term, c("Average Temp"="Baseline Climate"))
GainsMar_coeffsample$Term<- plyr::revalue(GainsMar_coeffsample$Term, c("Average Temp"="Baseline Climate"))
LossesMar_coeffsample$Term<- plyr::revalue(LossesMar_coeffsample$Term, c("Average Temp"="Baseline Climate"))


###the code below produces the four panels in Fig. S5

###for RICHNESS
ggplot(data = RichMar_coeffsample, aes(x = Term, y = fixed.Estimate, group = Term, label = Term)) + 
  geom_pointrange(aes(ymin = fixed.l.95..CI, ymax = fixed.u.95..CI), position=position_jitter(w=.4), color="grey70") +
  scale_x_discrete(limits=c("Interaction", "Baseline Climate", "Temp change")) +   ##order of plot for x axis 
  coord_flip() +
  labs(y="Coefficient", x= "", title = "Richness") +
  
  geom_hline(yintercept=0, lty=2) +
  
  theme(text = element_text(size= 15),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.text.y=element_text(size = 1),
        axis.line.x = element_line(color="black", size=.8), 
        axis.line.y = element_line(color="black", size=.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  
  ###adding point for each "true" paramereter --getting values from newcoeffALL
  geom_pointrange(data= newcoeffALL %>%
                    filter(Metric == "Richness", Realm=="Marine"), aes(x=Term, y=Estimate, ymin= Q2.5, ymax=Q97.5, color=Realm),
                  lwd=1.8, fatten = 3.5) +
  scale_color_manual(values = c("#0073C2FF", "#5da93d"))



###for ABUNDANCE
ggplot(data = AbundMar_coeffsample, aes(x = Term, y = fixed.Estimate, group = Term, label = Term)) + 
  geom_pointrange(aes(ymin = fixed.l.95..CI, ymax = fixed.u.95..CI), position=position_jitter(w=.4), color="grey70") +
  scale_x_discrete(limits=c("Interaction", "Baseline Climate", "Temp change")) +   ##order of plot for x axis 
  coord_flip() +
  labs(y="Coefficient", x= "", title = "Abundance") +
  
  geom_hline(yintercept=0, lty=2) +
  
  theme(text = element_text(size= 15),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.text.y=element_text(size = 1),
        axis.line.x = element_line(color="black", size=.8), 
        axis.line.y = element_line(color="black", size=.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  
  ###adding point for each "true" paramereter --getting values from newcoeffALL
  geom_pointrange(data= newcoeffALL %>%
                    filter(Metric == "Abundance", Realm=="Marine"), aes(x=Term, y=Estimate, ymin= Q2.5, ymax=Q97.5, color=Realm),
                  lwd=1.8, fatten = 3.5) +
  scale_color_manual(values = c("#0073C2FF", "#5da93d")) +
  
  ##added to better plotting
  ylim(-2.5, NA)

AbundMar_coeffsample %>% 
  group_by(Term) %>% 
  summarise(mean(fixed.Estimate))



###for GAINS
ggplot(data = GainsMar_coeffsample, aes(x = Term, y = fixed.Estimate, group = Term, label = Term)) + 
  geom_pointrange(aes(ymin = fixed.l.95..CI, ymax = fixed.u.95..CI), position=position_jitter(w=.4), color="grey70") +
  scale_x_discrete(limits=c("Interaction", "Baseline Climate", "Temp change")) +   ##order of plot for x axis 
  coord_flip() +
  labs(y="Coefficient", x= "", title = "Gains") +
  
  geom_hline(yintercept=0, lty=2) +
  
  theme(text = element_text(size= 15),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.text.y=element_text(size = 1),
        axis.line.x = element_line(color="black", size=.8), 
        axis.line.y = element_line(color="black", size=.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  
  ###adding point for each "true" paramereter --getting values from newcoeffALL
  geom_pointrange(data= newcoeffALL %>%
                    filter(Metric == "Gains", Realm=="Marine"), aes(x=Term, y=Estimate, ymin= Q2.5, ymax=Q97.5, color=Realm),
                  lwd=1.8, fatten = 3.5) +
  scale_color_manual(values = c("#0073C2FF", "#5da93d"))

GainsMar_coeffsample %>% 
  group_by(Term) %>% 
  summarise(mean(fixed.Estimate))



###for LOSSES
ggplot(data = LossesMar_coeffsample, aes(x = Term, y = fixed.Estimate, group = Term, label = Term)) + 
  geom_pointrange(aes(ymin = fixed.l.95..CI, ymax = fixed.u.95..CI), position=position_jitter(w=.4), color="grey70") +
  scale_x_discrete(limits=c("Interaction", "Baseline Climate", "Temp change")) +   ##order of plot for x axis 
  coord_flip() +
  labs(y="Coefficient", x= "", title = "Losses") +
  
  geom_hline(yintercept=0, lty=2) +
  
  theme(text = element_text(size= 15),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.text.y=element_text(size = 1),
        axis.line.x = element_line(color="black", size=.8), 
        axis.line.y = element_line(color="black", size=.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  
  ###adding point for each "true" paramereter --getting values from newcoeffALL
  geom_pointrange(data= newcoeffALL %>%
                    filter(Metric == "Losses", Realm=="Marine"), aes(x=Term, y=Estimate, ymin= Q2.5, ymax=Q97.5, color=Realm),
                  lwd=1.8, fatten = 3.5) +
  scale_color_manual(values = c("#0073C2FF", "#5da93d"))

LossesMar_coeffsample %>% 
  group_by(Term) %>% 
  summarise(mean(fixed.Estimate))




