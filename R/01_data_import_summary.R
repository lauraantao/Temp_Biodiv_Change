##======================================================================
##	The biodiversity data used here originates from pre-processing steps
##	described in detail elsewhere.

##  Biodiversity data are from the BioTIME gridded data (each time series =="rarefy_ID"), 
##  selecting temperate time series with at least 5 years sampled

##  Linear slopes were calculated for different biodiversity metrics for each time series:
##  model_id== "logS_lm", "logN_lm", "Gains_lm", "Losses_lm"

##  Temperature values were extracted for each rarefy_ID (monthly records for the duration of each time series)
##  GAM models were used to calculate the trends for temperature change

##  The data imported here is a dataframe with a slope estimate for each biodiversity metric 
##  && a slope estimate for temperature change per time series
##  Each rarefy_ID is nested within Study_ID and REALM

# Primary Code Authors: Laura Antao
# Code to process biodiversity and temperature data was built on extensive previous code developed by Shane Blowes, Sarah Supp and Conor Waldock
# Email: laura.antao@helsinki.fi
# DISCLAIMER: code is under development
# 31/01/2020


## Script (01_) imports the data, produces maps and summary plots for the biodiversity and temperature trends, 
## as well as plot summarising time series duration, number of years sampled, etc. (e.g. Figs. 2, S1 and S2).


# Input is Rdata object with biodiversity change estimates, temperature change estimates, and metadata for the time series


##======================================================================

library(tidyverse)
library(brms)
library(ggpubr)
library(ggthemes)
library(viridis)
library(ggExtra)
library(cowplot)

#setwd("~/Temp_Biodiv_Change/R")

## Load data
load("lm_slopes_meta_temperature3.Rdata")


n_distinct(lm_slopes_meta_temperature3$rarefyID) #21500 time series
n_distinct(lm_slopes_meta_temperature3$STUDY_ID) #156 original Study IDs
n_distinct(lm_slopes_meta_temperature3$taxa_mod1)  #9 taxonomic groups

# names(lm_slopes_meta_temperature3)
###explaining the variables in the data:
## STUDY_ID == the unique datasets in BioTIME (refers to Table S1 in Dornelas et al. (2018))
## REALM == Marine or Terrestrial
## TAXA == original taxonomic groups as described in BioTIME metadata
## TITLE == title of the original studies
## model_id == indicates the metric for which the trends were estimated (Richness="logS_lm", Abundance="logN_lm", Gains="Gains_lm", Losses="Losses_lm")
## rarefyID == corresponds to the unique combination of study ID and grid cell, i.e. standardised time series that will be analysed
## slope == estimated biodiversity trend over time (for each biodiversity metric, via lm)
## std.error == estimated standard error (for each slope estimate)
## p.value == estimated p.value of the trend
## change_mu == indicates "neutral", "up" or "down" trends according to the slope and p.values
## rarefyID_x == longitude of each grid cell
## rarefyID_y == latitude of each grid cell
## climate_mod == latitudinal band based on latitudinal cut-offs (in this case all data are "Temperate")
## startYear == start year of the time series
## endYear == end year of the time series
## num_years == number of years sampled in the time series
## duration == time series duration (i.e. Yend-Ystart+1)
## StudyMethod == classification based on samples locations and study extent
## TempGAMCoef == estimated temperature trend (for each time series (i.e. rarefy_ID), via GAM)
## TempGAMCoefSE == estimated standard error (for each slope estimate)
## changegam == indicates "up" or "down" trends according to the slope
## taxa_mod1 == simplified taxonomic groups
## newtempvalues == long-term annual average temperature values (i.e. baseline climate, per rarefy_ID, extracted from global databases)
## MAXtempvalues == long-term annual maximum temperature values (i.e. baseline climate, per rarefy_ID, extracted from global databases)
## new_sTempYear == standardised long-term annual average temperature values
## MAX_sTempYear == standardised long-term annual maximum temperature values
## model_clean == better model names





# create subsets for each biodiv metric
##for species richness
richnessdata <- lm_slopes_meta_temperature3 %>%
  filter(model_id=="logS_lm")

n_distinct(richnessdata$rarefyID) #21500 time series
n_distinct(richnessdata$STUDY_ID) #156 original Study IDs
n_distinct(richnessdata$taxa_mod1)  #9 taxonomic groups

##for total abundance
abundancedata <- lm_slopes_meta_temperature3 %>%
  filter(model_id=="logN_lm")

##for species gains
gainsdata <- lm_slopes_meta_temperature3 %>%
  filter(model_id=="Gains_lm")

##for species losses
lossesdata <- lm_slopes_meta_temperature3 %>%
  filter(model_id=="Losses_lm")

##number of time series per taxa and realm
richnessdata %>% 
  group_by(taxa_mod1,REALM) %>%
  summarise(ntimeseries=n()) %>% 
  arrange(REALM)


##======================================================================

###illustrating the estimated changes in temperature and in biodiversity across the time series
##this code produces Fig 2 in the manuscript


##temperature change ~ latitude panel
(temp_lat_change<- ggplot(richnessdata, aes(abs(rarefyID_y), TempGAMCoef)) +
   facet_wrap(~REALM, scales = "free_x") +
   geom_point(aes(color=newtempvalues), size=2) +
   coord_flip() +
   theme(text = element_text(size= 15),
         axis.text = element_text(size = 15), 
         axis.title = element_text(size = 15),
         axis.line.x = element_line(color="black", size=.8), 
         axis.line.y = element_line(color="black", size=.8),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black"),
         legend.text = element_text(size = 12),
         strip.background = element_rect( fill = "white"),
         strip.text.x = element_text( face = "bold"),
         legend.position = c(0.9,0.15),
         legend.direction = "horizontal",
         legend.title=element_text(size=12)) +
   geom_hline(aes(yintercept = 0), color = "grey20", linetype = "dashed") +  ##zero line
   
   ##adding meridian lines
   geom_vline(aes(xintercept = 23.5), color = "grey75") +
   geom_vline(aes(xintercept = 60), color = "grey75") +
   
   labs(y = "", x ="Absolute Latitude", colour="Long-term\nTemperature") +
   scale_x_continuous(breaks = c(0, 23.5, 60))+
   scale_color_viridis(option = "plasma", guide = "colourbar") +
  guides(colour = guide_colourbar(title.position = "top")))


##temperature density plot
auxtemp<- richnessdata %>%
  group_by(REALM) %>%
  summarise(mean=mean(TempGAMCoef))

(tempdens<- ggplot(richnessdata, aes(x=TempGAMCoef)) +
    geom_density(aes(colour = REALM, fill = REALM), alpha=0.5) + ###slopes of temperature change for each realm
    scale_fill_manual(values = c("#0073C2FF", "#5da93d")) + 
    scale_color_manual(values = c("#0073C2FF", "#5da93d")) +
    theme(text = element_text(size= 15),
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 15),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, "cm")) +
    labs(x = "Temperature change", y="", fill= "Realm", color="Realm") +
    geom_vline(xintercept=0, lty=2) +
    geom_vline(aes(xintercept = mean), data= auxtemp, size=1.1, lty=c(2, 6), colour= c("blue", "forestgreen")) +  ##add mean line per realm
    xlim(c(-0.2,0.2)))  ##cuttig off the extreme values from the density plot



####creating density plots for each biodiv metric

##richness density plot
auxrich<- richnessdata%>%
  group_by(REALM) %>%
  summarise(mean=mean(slope))

(richdens<- ggplot(richnessdata, aes(x=slope)) +
    geom_density(aes(colour = REALM, fill = REALM), alpha=0.5) +
    scale_fill_manual(values = c("#0073C2FF", "#5da93d")) + 
    scale_color_manual(values = c("#0073C2FF", "#5da93d")) +
    theme(text = element_text(size= 15),
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 15),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 1)) +
    labs(title = "Richness", x="", y="", fill= "Realm", color="Realm") +
    geom_vline(xintercept=0, lty=2) + ##add zero line
    geom_vline(aes(xintercept = mean), data = auxrich, size=1.1, lty=c(2, 6), colour= c("blue", "forestgreen")) +  ##add mean value per realm
    xlim(c(-0.2,0.2)))   ##cuttig off the extreme values from the density plot



##abundance density plot
auxabund<- abundancedata%>%
  group_by(REALM) %>%
  summarise(mean=mean(slope))

(abundens<- ggplot(abundancedata, aes(x=slope)) +
    geom_density(aes(colour = REALM, fill = REALM), alpha=0.5) +
    scale_fill_manual(values = c("#0073C2FF", "#5da93d")) + 
    scale_color_manual(values = c("#0073C2FF", "#5da93d")) +
    theme(text = element_text(size= 15),
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 15),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 1)) +
    labs(title = "Abundance", x="", y="", fill= "Realm", color="Realm") +
    geom_vline(xintercept=0, lty=2) +
    geom_vline(aes(xintercept = mean), data = auxabund, size=1.1, lty=c(2, 6), colour= c("blue", "forestgreen")) +  ##add mean line per realm
    xlim(c(-.5,.5)))   ##cuttig off the extreme values from the density plot



##gains density plot
auxgains<- gainsdata%>%
  group_by(REALM) %>%
  summarise(mean=mean(slope))

(gaindens<- ggplot(gainsdata, aes(x=slope)) +
  geom_density(aes(colour = REALM, fill = REALM), alpha=0.5) +
  scale_fill_manual(values = c("#0073C2FF", "#5da93d")) + 
  scale_color_manual(values = c("#0073C2FF", "#5da93d")) +
  theme(text = element_text(size= 15),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.line.x = element_line(color="black", size=.8), 
        axis.line.y = element_line(color="black", size=.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 1)) +
  labs(title = "Gains", x="", y="", fill= "Realm", color="Realm") +
  geom_vline(xintercept=0, lty=2) +
  geom_vline(aes(xintercept = mean), data = auxgains, size=1.1, lty=c(2, 6), colour= c("blue", "forestgreen")) +  ##add mean line per realm
  xlim(c(-2,2)))   ##cuttig off the extreme values from the density plot



##losses density plot
auxlosses<- lossesdata%>%
  group_by(REALM) %>%
  summarise(mean=mean(slope))

(lossdens<- ggplot(lossesdata, aes(x=slope)) +
    geom_density(aes(colour = REALM, fill = REALM), alpha=0.5) +
    scale_fill_manual(values = c("#0073C2FF", "#5da93d")) + 
    scale_color_manual(values = c("#0073C2FF", "#5da93d")) +
    theme(text = element_text(size= 15),
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 15),
          axis.line.x = element_line(color="black", size=.8), 
          axis.line.y = element_line(color="black", size=.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 1)) +
    labs(title = "Losses", x="Rates of change", y="", fill= "Realm", color="Realm") +
    geom_vline(xintercept=0, lty=2) +
    geom_vline(aes(xintercept = mean), data = auxlosses, size=1.1, lty=c(2, 6), colour= c("blue", "forestgreen")) +  ##add mean line per realm
    xlim(c(-1,1)))   ##cuttig off the extreme values from the density plot


##combine them
a<- plot_grid(temp_lat_change, tempdens, ncol = 1, rel_heights = c(3, 1))
b<- plot_grid(richdens, abundens, gaindens, lossdens, ncol = 1)


##to produce Fig 2
plot_grid(a, b, ncol = 2, rel_widths = c(2, 1))




##======================================================================

#####create maps with time series showing where the data come from
###also illustrating the estimated change in temperature and then in species richness
##this code produces Fig. S1 in the manuscript


#####Fig  S1####

World <- map_data("world")

# pdf(width = 7.25, height = 8, file = "FigS1.pdf")


##for temperature change
##Marine
Tchange_map_MAR <- ggplot() +
  geom_polygon(data=World, aes(long, lat, group = group), fill="#7f7f7f", size=0, alpha=0.6) +
  geom_point(data= arrange(richnessdata, TempGAMCoef) %>%
               filter(REALM== "Marine"),
             aes(x=rarefyID_x, y=rarefyID_y, colour= TempGAMCoef), size=1.8, alpha= 0.7) +  #shape=REALM
  scale_color_viridis(option = "inferno", guide = "colourbar") +
  labs(x= "", y="Temperature change", colour="Temperature\nchange", #shape="Realm",
       #tag = "Temperature change",
       subtitle = "Marine") +
  theme_bw() +
  scale_size_area() +
  coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme(#panel.grid.minor = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.position = c(0.1,0.3),
    legend.direction = "horizontal", 
    legend.key.size = unit(6, "mm"),
    legend.spacing.x=unit(0, "mm"),
    plot.margin = margin(0,0,0,0, "mm"),
    legend.background = element_rect(fill = "transparent"),
    legend.text.align = 0.9,
    plot.subtitle = element_text(size = 10, face = "bold"))  +
  guides(colour = guide_colourbar(title.position = "top")) #shape =  guide_legend(title.position = "top", override.aes = list(size = 3))


##Terrestrial
Tchange_map_TER <- ggplot() +
  geom_polygon(data=World, aes(long, lat, group = group), fill="#7f7f7f", size=0, alpha=0.6) +
  geom_point(data= arrange(richnessdata, TempGAMCoef) %>%
               filter(REALM== "Terrestrial"),
             aes(x=rarefyID_x, y=rarefyID_y, colour= TempGAMCoef), size=1.8, alpha= 0.7) +  ##shape=REALM
  scale_color_viridis(option = "inferno", guide = "colourbar") +
  labs(x= "", y="", colour="Temperature\nchange", #shape="Realm",
       subtitle = "Terrestial") +
  theme_bw() +
  scale_size_area() +
  coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme(#panel.grid.minor = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.position = c(0.1,0.3),
    legend.direction = "horizontal", 
    legend.key.size = unit(6, "mm"),
    legend.spacing.x=unit(0, "mm"),
    plot.margin = margin(0,0,0,0, "mm"),
    legend.background = element_rect(fill = "transparent"),
    legend.text.align = 0.9,
    plot.subtitle = element_text(size = 10, face = "bold"))  +
  guides(colour = guide_colourbar(title.position = "top")) #shape =  guide_legend(title.position = "top", override.aes = list(size = 3))


# pdf(width = 7.25, height = 8, file = "FigS1mar_ter.pdf")
# png(width = 11, height = 12, file = "Fig1.png", units = "in", res = 150)


##for species richness
##Marine
SRchange_map_MAR <- ggplot() +
  geom_polygon(data=World, aes(long, lat, group = group), fill="#7f7f7f", size=0, alpha=0.6) +
  geom_point(data= arrange(richnessdata, slope) %>%
               filter(REALM== "Marine"),
             aes(x=rarefyID_x, y=rarefyID_y, colour=slope), size=1.8, alpha= 0.7) +
  scale_color_viridis(option = "magma") +
  labs(x= "", y="Species richness change", colour="Species richness\nchange", #shape="",
       subtitle = "") +
  theme_bw() +
  scale_size_area() +
  coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme(#panel.grid.minor = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.position = c(0.1,0.3),
    legend.direction = "horizontal", 
    legend.key.size = unit(6, "mm"),
    legend.spacing.x=unit(3, "mm"),
    plot.margin = margin(0,0,0,0, "mm"),
    legend.background = element_rect(fill = "transparent"),
    legend.text.align = 0.9,
    plot.subtitle = element_text(size = 10, face = "bold"))  +
  guides(shape = guide_legend(title.position = "top", override.aes = list(size = 3)),
         colour = guide_colourbar(title.position = "top"))


##Terrestrial
SRchange_map_TER <- ggplot() +
  geom_polygon(data=World, aes(long, lat, group = group), fill="#7f7f7f", size=0, alpha=0.6) +
  geom_point(data= arrange(richnessdata, slope) %>%
               filter(REALM== "Terrestrial"),
             aes(x=rarefyID_x, y=rarefyID_y, colour=slope), size=1.8, alpha= 0.7) +
  scale_color_viridis(option = "magma") +
  labs(x= "", y="", colour="Species richness\nchange", #shape="",
       subtitle = "") +
  theme_bw() +
  scale_size_area() +
  coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme(#panel.grid.minor = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.position = c(0.1,0.3),
    legend.direction = "horizontal", 
    legend.key.size = unit(6, "mm"),
    legend.spacing.x=unit(3, "mm"),
    plot.margin = margin(0,0,0,0, "mm"),
    legend.background = element_rect(fill = "transparent"),
    legend.text.align = 0.9,
    plot.subtitle = element_text(size = 10, face = "bold"))  +
  guides(shape = guide_legend(title.position = "top", override.aes = list(size = 3)),
         colour = guide_colourbar(title.position = "top"))


#pdf(width = 7.25, height = 10, file = "FigS1mar_ter.pdf")
plot_grid(Tchange_map_MAR, Tchange_map_TER,
          SRchange_map_MAR, SRchange_map, ncol = 2)
#dev.off()



##======================================================================

##to show the temporal heterogeneity across time series
##this code produces Fig. S2 in the manuscript


#####Fig  S2####

(hist_startYear <- richnessdata %>%
   ggplot() +
   geom_histogram(aes(x = startYear), binwidth = 2) +
   scale_y_continuous(breaks = c(0,2500,5000,7500)) +
   labs(tag = "a)",
        y = "Number of time series",
        x = "First year sampled") +
   theme_bw() +
   theme(panel.grid.minor.x = element_blank()))



(hist_numYears <- richnessdata %>%
    ggplot() +
    geom_histogram(aes(x = num_years), binwidth = 2) +
    scale_x_continuous(name = "Number of years sampled",
                       breaks = c(5,10,15,20,30,50,95)) +
    scale_y_continuous(breaks = c(0,2500,5000,7500,10000,12500, 15000, 17500, 20000,22500)) +
    labs(tag = "b)",
         y = "Number of time series") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank()))


(hist_duration <- richnessdata %>%
    ggplot() +
    geom_histogram(aes(x = duration), binwidth = 1) +
    scale_x_continuous(name = "Duration of study",
                       breaks = c(5,10,15,20,30,50,95),
                       labels = c(5,10,15,20,30,50,95)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) + 
    labs(tag = "c)",
         y = "Number of time series"))


(xy_density <- richnessdata %>%
    ggplot() + 
    geom_bin2d(aes(x = num_years, y = duration)) +
    scale_x_continuous("Number of years sampled",
                       breaks = c(5,10,15,20,30,50,95)) +
    scale_y_continuous("Duration of study",
                       breaks = c(5,10,15,20,30,50,95),
                       labels = c(5,10,15,20,30,50,95)) +
    scale_fill_gradientn(name = "Number of time series", 
                         trans = "log2",
                         # breaks = c(2,16,128,1024,4112),
                         colours = c("#fdd49e",
                                     "#fdbb84", "#fc8d59",
                                     "#e34a33", "#b30000")) +
    labs(tag = "d)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = c(0.79, 0.3),
          legend.box.background = element_rect(linetype = 1, size = 0.25, fill = "white"),
          legend.background = element_blank()))

plot_grid(hist_startYear, hist_numYears, hist_duration, xy_density, nrow = 2, align = "hv")
#ggsave("FigS2.png", width = 240, height = 180, units = "mm")






