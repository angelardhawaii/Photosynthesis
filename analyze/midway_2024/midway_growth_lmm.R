#This script brings in growth data for Midway Chondria tumulosa analysis
#By Angela Richards Donà
#Date: 9/13/24

# open packages for mixed model effects analysis
library(lme4)
library(lmerTest)
library(effects)
library(car)
library(MuMIn)
library (dplyr)
library(emmeans)
library(DHARMa)
library(performance)
library(patchwork)
library(rstatix)
#for plots and tables
library(ggplot2)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(strengejacke)
#library(mmtable2)
library(gt)
library(purrr)
library(stringr)
library(tidyr)
library(piecewiseSEM)
library(easystats)
library(magrittr)

#open weight dataset and make columns for growth rate from initial and final weights
mid_growth <- read.csv("../../../data/limu/midway/input/growth/midway_growth.csv")

#make a new column for weight change (difference final from initial)
mid_growth$d9_growth_percent <- round(((mid_growth$final_weight - mid_growth$initial_weight) / mid_growth$initial_weight * 100), digits = 2)
mid_growth$d5_growth_percent <- round(((mid_growth$d5_weight - mid_growth$initial_weight) / mid_growth$initial_weight * 100), digits = 2)

#add a new column that gets rid of characters in ID column
mid_growth$plant_ID <- as.factor(substr(mid_growth$ID, 3, 5))

#make a new column for run using the third value in the ID
mid_growth$run <- as.factor(substr(mid_growth$ID, 3, 3))

#assigns treatment as characters from integers then to factors
mid_growth$treatment <- as.factor(as.character(mid_growth$treatment))

#remove ppt from salinity data cells
mid_growth$salinity <- substr(mid_growth$salinity, 1, 2)

#create subsets for the plots
#toggle between the plant_part for output. Use Day 9 for final analysis
canopy <- subset(mid_growth, plant_part == "canopy")
canopy$treatment_graph[canopy$treatment == 1] <- "0.5umol"
canopy$treatment_graph[canopy$treatment == 2] <- "2umol" 
canopy$treatment_graph[canopy$treatment == 3] <- "4umol"
canopy$treatment_graph[canopy$treatment == 4] <- "8umol" 

under <- subset(mid_growth, plant_part == "under" & d9_growth_percent != -62.82)
under$treatment_graph[under$treatment == 1] <- "0.5umol"
under$treatment_graph[under$treatment == 2] <- "2umol" 
under$treatment_graph[under$treatment == 3] <- "4umol"
under$treatment_graph[under$treatment == 4] <- "8umol"


#-----------make a histogram of the growth rate data
under %>% ggplot(aes(d5_growth_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()
under %>% ggplot(aes(d9_growth_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

canopy %>% ggplot(aes(d5_growth_percent)) +
  geom_histogram(binwidth=5, fill = "goldenrod1", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()
canopy %>% ggplot(aes(d9_growth_percent)) +
  geom_histogram(binwidth=5, fill = "goldenrod1", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run model without RLC_order as this has little effect
#add number of initial axes to the random effects
growth_model_canopy <- lmer(formula = d9_growth_percent ~ treatment +
                            (1 | plant_ID) + (1 | run), data = canopy, REML = TRUE)

hist(resid(growth_model_acan))
plot(resid(growth_model_acan) ~ fitted(growth_model_acan))
qqnorm(resid(growth_model_acan))
qqline(resid(growth_model_acan))

#check the performance of the model for dataset: acan
performance::check_model(growth_model_acan)
rsquared(growth_model_acan)
summary(growth_model_acan)
#view random effects levels
ranef(growth_model_acan)
tab_model(growth_model_acan, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_acan))

#construct null model to perform likelihood ratio test REML must be FALSE
mid_growth_treatment_null <- lmer(formula = growth_rate_percent ~ temperature + (1 | plant_ID) + (1 | run), data = mid_growth, REML = FALSE)
mid_growth_model2 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant_ID) + (1 | run), data = mid_growth, REML = FALSE)
anova(mid_growth_treatment_null, mid_growth_model2)
mid_growth_temperature_null <- lmer(formula = growth_rate_percent ~ treatment + (1 | plant_ID) + (1 | run), data = mid_growth, REML = FALSE)
mid_growth_model3 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant_ID) + (1 | run), data = mid_growth, REML = FALSE)
anova(mid_growth_temperature_null, mid_growth_model3)

#plots
mid_growth %>% ggplot(aes(treatment_graph, growth_rate_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = temperature), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "A", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("T0) 35 ppt/0.5 μmol N", "T1) 28 ppt/14 μmol N")) + 
  ylim(-60, 60) + stat_mean() + 
  scale_color_manual(values = c("azure4", "darkgoldenrod1")) +
  geom_hline(yintercept=0, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means
mid_growth %>% group_by(treatment) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))
mid_growth %>% group_by(temperature) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))

#plot temperature
mid_growth %>% ggplot(aes(temperature_graph, growth_rate_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="Temperature", y= "9-day Growth (%)", title= "B", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("24°C", "28°C")) + 
  ylim(-60, 50) + stat_mean() + 
  geom_hline(yintercept=0, color = "red", linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("azure3", "darkgoldenrod")) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))





#RUN MODEL FOR SECONDARY APICES
mid_growth %>% ggplot(aes(secondary_apices)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()
#run model without RLC_order as this has little effect
apices_model_acan <- lmer(formula = secondary_apices ~ treatment + temperature +
                            (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = mid_growth, REML = TRUE)

hist(resid(apices_model_acan))
plot(resid(apices_model_acan) ~ fitted(apices_model_acan))
qqnorm(resid(apices_model_acan))
qqline(resid(apices_model_acan))

#check the performance of the model for dataset: acan
performance::check_model(apices_model_acan)
rsquared(apices_model_acan)
summary(apices_model_acan)
#view random effects levels
ranef(apices_model_acan)
tab_model(apices_model_acan, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(apices_model_acan))

#construct null model to perform likelihood ratio test REML must be FALSE
acan_apices_treatment_null <- lmer(formula = secondary_apices~ temperature + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = mid_growth, REML = FALSE)
acan_apices_model2 <- lmer(formula = secondary_apices ~ treatment + temperature + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = mid_growth, REML = FALSE)
anova(acan_apices_treatment_null, acan_apices_model2)
acan_apices_temperature_null <- lmer(formula = secondary_apices ~ treatment + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = mid_growth, REML = FALSE)
acan_apices_model3 <- lmer(formula = secondary_apices ~ treatment + temperature + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = mid_growth, REML = FALSE)
anova(acan_apices_temperature_null, acan_apices_model3)

#plots
mid_growth %>% ggplot(aes(treatment_graph, secondary_apices)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = temperature), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "Number of Secondary Apices", title= "A", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("T0) 35 ppt/0.5 μmol N", "T1) 28 ppt/14 μmol N")) + 
  ylim(-1, 130) + stat_mean() + 
  scale_color_manual(values = c("azure4", "darkgoldenrod1")) +
  geom_hline(yintercept=0, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means
mid_growth %>% group_by(treatment) %>% summarise_at(vars(secondary_apices), list(mean = mean))
mid_growth %>% group_by(temperature) %>% summarise_at(vars(secondary_apices), list(mean = mean))


#plot temperature
mid_growth %>% ggplot(aes(temperature_graph, secondary_apices)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="Temperature", y= "Number of Secondary Apices", title= "B", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("24°C", "28°C")) + 
  ylim(-1, 130) + stat_mean() + 
  geom_hline(yintercept=0, color = "red", linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("azure3", "darkgoldenrod")) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

