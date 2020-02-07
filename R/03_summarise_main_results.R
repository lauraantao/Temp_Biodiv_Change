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


# setwd("~/Temp_Biodiv_Change/R")

##load model results objects
load("model_fits_output/NewRichTer_RS.Rdata")
load("model_fits_output/NewRichMar_RS.Rdata")
load("model_fits_output/NewAbundTer_RS.Rdata")
load("model_fits_output/NewAbundMar_RS.Rdata")
load("model_fits_output/NewGainsTer_RS.Rdata")
load("model_fits_output/NewGainsMar_RS.Rdata")
load("model_fits_output/NewLossTer_RS.Rdata")
load("model_fits_output/NewLossMar_RS.Rdata")



#####Fig 3 ####

##combine model estimates from all the metrics into a single data.frame
newcoeffALL <- rbind(data.frame(rbind(fixef(NewRichMar_RS), fixef(NewRichTer_RS)),
                                Metric =rep("Richness",6),
                                Realm =rep(c("Marine","Terrestrial" ),each=3)),
                     
                     data.frame(rbind(fixef(NewAbundMar_RS), fixef(NewAbundTer_RS)),
                                Metric =rep("Abundance",6),
                                Realm = rep(c("Marine","Terrestrial" ),each=3)),
                     
                     data.frame(rbind(fixef(NewGainsMar_RS), fixef(NewGainsTer_RS)),
                                Metric =rep("Gains",6),
                                Realm = rep(c("Marine","Terrestrial" ),each=3)),
                     
                     data.frame(rbind(fixef(NewLossMar_RS), fixef(NewLossTer_RS)),
                                Metric =rep("Losses",6),
                                Realm = rep(c("Marine","Terrestrial" ),each=3))) %>%
  mutate(Term= rep(c("Temp change", "Baseline Climate", "Interaction"),8))


##to produce figure
ggplot(data= newcoeffALL, #%>%
       aes(x=Term, y=Estimate, fill=Realm)) + 
  facet_wrap(~ Metric, ncol=4, scales = "free_x")+ 
  scale_x_discrete(limits=c("Interaction", "Baseline Climate", "Temp change"),   ##order of plot for x axis
                   labels = c("Interaction" = "Interaction",
                              "Baseline Climate" = "Baseline\nClimate",
                              "Temp change" = "Temperature\nchange")) +       ##rename variables 
  coord_flip() +
  geom_bar(position=position_dodge(), stat="identity", alpha=0.8, width=.7) +
  geom_errorbar(aes(ymin= Q2.5, ymax= Q97.5, color=Realm),
                width=.15,                    # Width of the error bars
                position=position_dodge(.7)) +
  labs(x="", y="", title = "") +
  scale_fill_manual(values = c("#0073C2FF", "#5da93d"))+ 
  scale_color_manual(values = c("#2e3e74", "#3e7129"))+  ##darker colors for the CI lines
  theme(axis.text.x = element_text(size= 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        #axis.text.y = element_text(margin = unit(c(t = 0, r = -2.6, b = 0, l = 0), "cm"), hjust = 0, face = "bold"),
        axis.text.y = element_text(hjust = 1, size= 8),
        strip.background = element_rect(fill = "white"),
        #strip.text.x = element_blank(),
        panel.spacing = unit(2, "lines")) +
  geom_hline(yintercept=0, lty=2)


##or ommiting the main effect for baseline climate
ggplot(data= newcoeffALL %>%
         filter(Term != "Baseline Climate"),
       aes(x=Term, y=Estimate, fill=Realm)) + 
  facet_wrap(~ Metric, ncol=4, scales = "free_x")+ 
  scale_x_discrete(limits=c("Interaction",  "Temp change"),   ##order of plot for x axis
                   labels = c("Interaction" = "Interaction\nBaseline Climate",
                              #"Baseline Climate" = "Baseline\nClimate",
                              "Temp change" = "Temperature\nchange")) +       ##rename variables 
  coord_flip() +
  geom_bar(position=position_dodge(), stat="identity", alpha=0.8, width=.7) +
  geom_errorbar(aes(ymin= Q2.5, ymax= Q97.5, color=Realm),
                width=.15,                    # Width of the error bars
                position=position_dodge(.7)) +
  labs(x="", y="", title = "") +
  scale_fill_manual(values = c("#0073C2FF", "#5da93d"))+ 
  scale_color_manual(values = c("#2e3e74", "#3e7129"))+  ##darker colors for the CI lines
  theme(axis.text.x = element_text(size= 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        #axis.text.y = element_text(margin = unit(c(t = 0, r = -2.6, b = 0, l = 0), "cm"), hjust = 0, face = "bold"),
        axis.text.y = element_text(hjust = 1, size= 8),
        strip.background = element_rect(fill = "white"),
        #strip.text.x = element_blank(),
        panel.spacing = unit(2, "lines")) +
  geom_hline(yintercept=0, lty=2)



##======================================================================

#####Fig 4 ####

####plotting the interaction surfaces

###for RICHNESS

##MARINE
meSRichMar <- marginal_effects(NewRichMar_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, method = "fitted", stype="raster")

(pmeSRichMar<- plot(meSRichMar, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Marine Richness") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 12),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 5000, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= richnessdata %>%
                 filter(REALM=="Marine"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot=pmeSRichMar, "surface_richMar.tiff", width = 14, height = 12, units = c("cm"))




##TERRESTRIAL
meSRichTer <- marginal_effects(NewRichTer_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, stype="raster", method = "fitted")

(pmeSRichTer<- plot(meSRichTer, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Terrestrial Richness") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 12),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 50, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= richnessdata %>%
                 filter(REALM=="Terrestrial"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot = pmeSRichTer, "surface_richTer.tiff", width = 14, height = 12, units = c("cm"))



########
###for ABUNDANCE

##MARINE
meAbundMar <- marginal_effects(NewAbundMar_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, stype="raster", method = "fitted")

(pmeAbundMar<- plot(meAbundMar, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Marine Abundance") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 10),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 50, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= abundancedata %>%
                 filter(REALM=="Marine"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot = pmeAbundMar, "surface_abundMar.tiff", width = 14, height = 12, units = c("cm"))


##TERRESTRIAL
meAbundTer <- marginal_effects(NewAbundTer_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, stype="raster", method = "fitted")

(pmeAbundTer<- plot(meAbundTer, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Terrestrial Abundance") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 10),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 50, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= abundancedata %>%
                 filter(REALM=="Terrestrial"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot = pmeAbundTer, "surface_abundTer.tiff", width = 14, height = 12, units = c("cm"))




########
###for GAINS

##MARINE
meGainsMar <- marginal_effects(NewGainsMar_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, stype="raster", method = "fitted")

(pmeGainsMar<- plot(meGainsMar, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Marine Gains") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 10),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 500, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= gainsdata %>%
                 filter(REALM=="Marine"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot = pmeGainsMar, "surface_gainsMar.tiff", width = 14, height = 12, units = c("cm"))


##TERRESTRIAL
meGainsTer <- marginal_effects(NewGainsTer_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, stype="raster", method = "fitted")

(pmeGainsTer<- plot(meGainsTer, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Terrestrial Gains") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 10),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 50, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= gainsdata %>%
                 filter(REALM=="Terrestrial"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot = pmeGainsTer, "surface_gainsTer.tiff", width = 14, height = 12, units = c("cm"))




########
###for LOSSES

##MARINE
meLossesMar <- marginal_effects(NewLossMar_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, stype="raster", method = "fitted")

(pmeLossesMar<- plot(meLossesMar, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Marine Losses") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 10),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 500, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= lossesdata %>%
                 filter(REALM=="Marine"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot = pmeLossesMar, "surface_lossesMar.tiff", width = 14, height = 12, units = c("cm"))


##TERRESTRIAL
meLossesTer <- marginal_effects(NewLossTer_RS, "TempGAMCoef:new_sTempYear", surface = TRUE, stype="raster", method = "fitted")

(pmeLossesTer<- plot(meLossesTer, plot = FALSE, stype="raster")[[1]] + 
    labs(x = "", y="", colour = "", fill="Slope",
         subtitle = "Terrestrial Losses") +
    theme(axis.text = element_text(size = 18), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0.5, 0.7, -1, -1), "lines"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -20, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 10),
          legend.title=element_text(size= 18),
          legend.title.align = 1) +
    
    scale_color_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0) +
    guides(fill = guide_colourbar(nbin = 50, title.vjust=1)) +
    
    ##adding sampling points
    geom_point(data= lossesdata %>%
                 filter(REALM=="Terrestrial"),
               aes(x=TempGAMCoef, y=new_sTempYear), color="grey", alpha=0.5, shape=1))

ggsave(plot = pmeLossesTer, "surface_lossesTer.tiff", width = 14, height = 12, units = c("cm"))


##to combine plots for different metrics and realms
top <- plot_grid(pmeSRichMar, pmeAbundMar, pmeGainsMar, pmeLossesMar, nrow = 1)  ##marine
bottom <- plot_grid(pmeSRichTer, pmeAbundTer, pmeGainsTer, pmeLossesTer, nrow = 1)  ##terrestrial
plot_grid(top, bottom, nrow = 2)

# ggsave("Fig4_rev2.pdf", width = 120, height = 120, units = "mm")






##======================================================================


#####Fig S3 ####

####plotting marginal effects plots

###for RICHNESS

##MARINE
interRichMar<- plot(marginal_effects(NewRichMar_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interRichMar[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.78, 0.12),
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_richMar.tiff")


##TERRESTRIAL
interRichTer<- plot(marginal_effects(NewRichTer_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interRichTer[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.17, 0.17),
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_richTer.tiff")



########
###for ABUNDANCE

##MARINE
interAbndMar<- plot(marginal_effects(NewAbundMar_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interAbndMar[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_abundMar.tiff")


##TERRESTRIAL
interAbndTer<- plot(marginal_effects(NewAbundTer_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interAbndTer[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_abundTer.tiff")



########
###for GAINS

##MARINE
interGainMar<- plot(marginal_effects(NewGainsMar_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interGainMar[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_gainsMar.tiff")


##TERRESTRIAL
interGainTer<- plot(marginal_effects(NewGainsTer_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interGainTer[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_gainsTer.tiff")



########
###for LOSSES

##MARINE
interLossMar<- plot(marginal_effects(NewLossMar_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interLossMar[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_lossesMar.tiff")


##TERRESTRIAL
interLossTer<- plot(marginal_effects(NewLossTer_RS, "TempGAMCoef:new_sTempYear", method = "fitted")) 

plot(interLossTer[[1]]) +
  labs(x = "", y= "", color="Long-term\nTemperature", fill="Long-term\nTemperature") +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, -1, -1), "lines")) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2)

ggsave("inter_lossesTer.tiff")







