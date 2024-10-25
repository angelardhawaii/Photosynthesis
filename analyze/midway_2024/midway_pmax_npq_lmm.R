#Script to run model for photosynthesis data for Midway data on Chondria tumulosa
#By Angela Richards Donà
#Date: 9/13/24

library(lme4)
library(lmerTest)
library(afex)
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
#library(mmtable2)
library(gt)
library(purrr)
library(stringr)
library(tidyr)
library(hms)


#load this file for normalized to quantum efficiency of photosynthesis per Silsbe and Kromkamp
mid_ps <- read.csv("../../../data/limu/midway/transformed/midway_ek_alpha_normalized_2024.csv")

# make sure time is hms
mid_ps$rlc_end_time <- as_hms(mid_ps$rlc_end_time)

# assign run as a factor
mid_ps$run <- as.factor(mid_ps$run)

#assign temperature as a factor
mid_ps$salinity <- as.factor(mid_ps$salinity)

#assigns treatment as characters from integers then to factors
mid_ps$treatment <- as.factor(as.character(mid_ps$treatment))

# assign deltaNPQ as a numeric
mid_ps$deltaNPQ <- as.numeric(mid_ps$deltaNPQ)

#assign new column for chronological ranking of individuals in each day_group
mid_ps <- mid_ps %>%
  group_by(day_group) %>%
  mutate(rank = rank(rlc_end_time)) %>%
  ungroup()

#use ranking to make smaller groups of ~15 minutes for rlc_order
mid_ps <- mid_ps %>%
  group_by(day_group) %>%
  mutate(rlc_order = floor((rank -1)/3) +1)

#toggle between the plant_part for output. Use Day 9 for final analysis
canopy <- subset(mid_ps, rlc_day == 9 & plant_part == "canopy")
canopy$treatment_graph[mid_day9$treatment == 1] <- "1) 0.5umol"
canopy$treatment_graph[mid_day9$treatment == 2] <- "2) 2umol" 
canopy$treatment_graph[mid_day9$treatment == 3] <- "3) 4umol"
canopy$treatment_graph[mid_day9$treatment == 4] <- "4) 8umol" 

under <- subset(mid_ps, rlc_day == 9 & plant_part == "under")
under$treatment_graph[mid_day9$treatment == 1] <- "1) 0.5umol"
under$treatment_graph[mid_day9$treatment == 2] <- "2) 2umol" 
under$treatment_graph[mid_day9$treatment == 3] <- "3) 4umol"
under$treatment_graph[mid_day9$treatment == 4] <- "4) 8umol"


mid_day9$plantpart_graph[mid_day9$plant_part == "canopy"] <- "canopy"
mid_day9$plantpart_graph[mid_day9$plant_part == "under"] <- "under"


#------------Pmax----------------
#plot histogram
canopy %>% ggplot(aes(pmax)) +
  geom_histogram(binwidth=5, fill = "goldenrod1", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
under %>% ggplot(aes(pmax)) +
  geom_histogram(binwidth=5, fill = "maroon", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
canopy_pmax_model <- lmer(formula = pmax ~ treatment + salinity + (1 | run) + (1 | plantID) + (1 | rlc_order), data = canopy)
under_pmax_model <- lmer(formula = pmax ~ treatment + salinity + (1 | run) + (1 | plantID), data = under)

#construct null model to perform likelihood ratio test REML must be FALSE
canopy_pmax_treatment_null <- lmer(formula = pmax ~ salinity + (1 | run) + (1 | plantID) + (1 | rlc_order), data = canopy, REML = FALSE)
canopy_pmax_model2 <- lmer(formula = pmax ~ treatment + salinity + (1 | run) + (1 | plantID) + (1 | rlc_order), data = canopy, REML = FALSE)
anova(canopy_pmax_treatment_null, canopy_pmax_model2)
canopy_pmax_salinity_null <- lmer(formula = pmax ~ treatment + (1 | run) + (1 | plantID) + (1 | rlc_order), data = canopy, REML = FALSE)
canopy_pmax_model3 <- lmer(formula = pmax ~ treatment + salinity + (1 | run) + (1 | plantID) + (1 | rlc_order), data = canopy, REML = FALSE)
anova(canopy_pmax_salinity_null, canopy_pmax_model3)

under_pmax_treatment_null <- lmer(formula = pmax ~ salinity + (1 | run) + (1 | plantID), data = under, REML = FALSE)
under_pmax_model2 <- lmer(formula = pmax ~ treatment + salinity + (1 | run) + (1 | plantID), data = under, REML = FALSE)
anova(under_pmax_treatment_null, under_pmax_model2)
under_pmax_salinity_null <- lmer(formula = pmax ~ treatment + (1 | run) + (1 | plantID), data = under, REML = FALSE)
under_pmax_model3 <- lmer(formula = pmax ~ treatment + salinity + (1 | run) + (1 | plantID), data = under, REML = FALSE)
anova(under_pmax_salinity_null, under_pmax_model3)

#make residual plots of the data for Acanthophora
hist(resid(canopy_pmax_model))
plot(resid(canopy_pmax_model) ~ fitted(canopy_pmax_model))
qqnorm(resid(canopy_pmax_model))
qqline(resid(canopy_pmax_model))

hist(resid(under_pmax_model))
plot(resid(under_pmax_model) ~ fitted(canopy_pmax_model))
qqnorm(resid(under_pmax_model))
qqline(resid(under_pmax_model))

#check the performance of the model
performance ::check_model(canopy_pmax_model)
r.squaredGLMM(canopy_pmax_model)

performance ::check_model(under_pmax_model)
r.squaredGLMM(under_pmax_model)

#make plots and tables for the data
tab_model(canopy_pmax_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(canopy_pmax_model))

tab_model(under_pmax_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(under_pmax_model))

canopy %>% ggplot(aes(treatment_graph, pmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = salinity), position = "jitter", show.legend = TRUE) + 
  labs(x="nitrate", y= "Day 9 Pmax (μmols electrons m-2 s-1)", title= "A", subtitle = "Chondria tumulosa -- Canopy") + 
  scale_x_discrete(labels = c("0.5umolN", "2umolN", "4umolN", "8umolN")) + 
  ylim(25, 175) + stat_mean() + 
  scale_color_manual(values = c("gold", "goldenrod")) +
  geom_hline(yintercept=100, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

under %>% ggplot(aes(treatment_graph, pmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = salinity), position = "jitter", show.legend = TRUE) + 
  labs(x="nitrate", y= "Day 9 Pmax (μmols electrons m-2 s-1)", title= "B", subtitle = "Chondria tumulosa understory") + 
  scale_x_discrete(labels = c("0.5umolN", "2umolN", "4umolN", "8umolN")) + 
  ylim(0, 150) + stat_mean() + 
  scale_color_manual(values = c("maroon1", "maroon4")) +
  geom_hline(yintercept=100, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))


#summarize the means for pmax
canopy %>% group_by(treatment) %>% summarise_at(vars(pmax), list(mean = mean))
canopy %>% group_by(salinity) %>% summarise_at(vars(pmax), list(mean = mean))

under %>% group_by(treatment) %>% summarise_at(vars(pmax), list(mean = mean))
under %>% group_by(salinity) %>% summarise_at(vars(pmax), list(mean = mean))


#Linear regression Pmax vs Growth--------------------------------------------------------------------------------

#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("../../../data/limu/acan_2024/input/acanthophora_growth_r1.csv")
growth_rate$treatment <- as.factor(growth_rate$treatment)

#make a new column for weight change (difference final from initial)
growth_rate$growth_rate_percent <- (growth_rate$final_weight - growth_rate$initial_weight) / growth_rate$initial_weight * 100
mid_day9$growth_rate_percent <- growth_rate$growth_rate_percent
mid_day9$secondary_apices <- growth_rate$secondary_apices
mid_day9_sub <- subset(mid_day9, secondary_apices != "200")

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_growth_pmax_graph <- ggplot(mid_day9_sub, aes(x=pmax, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Pmax vs 9-Day Growth (%)", x = "Pmax (μmols electrons m-2 s-1)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
acan_growth_pmax_graph

#plot a regression between the photosynthetic independent variables of interest and apices/axis
acan_apices_pmax_graph <- ggplot(mid_day9_sub, aes(x=pmax, y=secondary_apices)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Pmax vs # Secondary Apices", x = "Pmax (μmols electrons m-2 s-1)", 
       y = "Secondary Apices (#)") + stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
acan_apices_pmax_graph





#--------------------NPQmax--------------------------

#plot histogram
mid_day9 %>% ggplot(aes(maxNPQ_Ypoint1)) +
  geom_histogram(binwidth=0.05, fill = "violet", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
mid_day9_maxNPQ_model <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + temp + (1 | plantID), data = mid_day9)

#construct null model to perform likelihood ratio test REML must be FALSE
acan_maxNPQ_treatment_null <- lmer(formula = maxNPQ_Ypoint1 ~ temp + (1 | plantID), data = mid_day9, REML = FALSE)
acan_maxNPQ_model2 <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + temp + (1 | plantID), data = mid_day9, REML = FALSE)
anova(acan_maxNPQ_treatment_null, acan_maxNPQ_model2)
acan_maxNPQ_temperature_null <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + (1 | plantID), data = mid_day9, REML = FALSE)
acan_maxNPQ_model3 <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + temp + (1 | plantID), data = mid_day9, REML = FALSE)
anova(acan_maxNPQ_temperature_null, acan_maxNPQ_model3)

#make residual plots of the data for Acanthophora
hist(resid(mid_day9_maxNPQ_model))
plot(resid(mid_day9_maxNPQ_model) ~ fitted(mid_day9_maxNPQ_model))
qqnorm(resid(mid_day9_maxNPQ_model))
qqline(resid(mid_day9_maxNPQ_model))

#check the performance of the model
performance ::check_model(mid_day9_maxNPQ_model)
r.squaredGLMM(mid_day9_maxNPQ_model)

#make plots and tables for the data
tab_model(mid_day9_maxNPQ_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(mid_day9_maxNPQ_model))

mid_day9 %>% ggplot(aes(treatment_graph, maxNPQ_Ypoint1)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = temp), position = "jitter", show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "Day 9 Maximum NPQ (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "28ppt/14umolN")) + 
  ylim(0, 1) + stat_mean() + 
  scale_color_manual(values = c("olivedrab2", "maroon2")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

mid_day9 %>% ggplot(aes(temperature_graph, maxNPQ_Ypoint1)) +
  geom_boxplot(size=0.5) +
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) +
  labs(x="temperature", y= "Day 9 Maximum NPQ (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") +
  ylim(0, 1) + stat_mean() + 
  scale_color_manual(values = c("olivedrab3", "maroon4")) +
  geom_hline(yintercept=0.5, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for pmax
mid_day9 %>% group_by(treatment) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))
mid_day9 %>% group_by(temp) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))
#mid_day9 %>% group_by(treatment, rlc_day) %>% summarise_at(vars(pmax), list(mean = mean))

#Linear regression NPQ vs Growth and Apices

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_growth_NPQmax_graph <- ggplot(mid_day9_sub, aes(x=maxNPQ_Ypoint1, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera NPQmax vs 9-Day Growth (%)", x = " NPQmax (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
acan_growth_NPQmax_graph

#plot a regression between the photosynthetic independent variables of interest and apices/axis
acan_apices_NPQmax_graph <- ggplot(mid_day9_sub, aes(x=maxNPQ_Ypoint1, y=secondary_apices)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera NPQmax vs # Secondary Apices", x = "NPQmax (rel. units)", 
       y = "Secondary Apices (#)") + stat_regline_equation(label.x = 0.25, label.y = 125) + stat_cor(label.x = 0.25, label.y = 128)
acan_apices_NPQmax_graph

#-----------------delta NPQ------------------
#plot histogram
mid_day9 %>% ggplot(aes(deltaNPQ)) +
  geom_histogram(binwidth=.05, fill = "violet", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
mid_day9_deltaNPQ_model <- lmer(formula = deltaNPQ ~ treatment + temp + (1 | plantID), data = mid_day9)

#construct null model to perform likelihood ratio test REML must be FALSE
acan_deltaNPQ_treatment_null <- lmer(formula = deltaNPQ ~ temp + (1 | plantID), data = mid_day9, REML = FALSE)
acan_deltaNPQ_model2 <- lmer(formula = deltaNPQ ~ treatment + temp + (1 | plantID), data = mid_day9, REML = FALSE)
anova(acan_deltaNPQ_treatment_null, acan_deltaNPQ_model2)
acan_deltaNPQ_temperature_null <- lmer(formula = deltaNPQ ~ treatment + (1 | plantID), data = mid_day9, REML = FALSE)
acan_deltaNPQ_model3 <- lmer(formula = deltaNPQ ~ treatment + temp + (1 | plantID), data = mid_day9, REML = FALSE)
anova(acan_deltaNPQ_temperature_null, acan_deltaNPQ_model3)

#make residual plots of the data for Acanthophora
hist(resid(mid_day9_deltaNPQ_model))
plot(resid(mid_day9_deltaNPQ_model) ~ fitted(mid_day9_deltaNPQ_model))
qqnorm(resid(mid_day9_deltaNPQ_model))
qqline(resid(mid_day9_deltaNPQ_model))

#check the performance of the model
performance ::check_model(mid_day9_deltaNPQ_model)
r.squaredGLMM(mid_day9_deltaNPQ_model)

#make plots and tables for the data
tab_model(mid_day9_deltaNPQ_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(mid_day9_deltaNPQ_model))

mid_day9 %>% ggplot(aes(treatment_graph, deltaNPQ)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = temp), position = "jitter", show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "Day 9 Delta NPQ (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "28ppt/14umolN")) + 
  ylim(0, 0.75) + stat_mean() + 
  scale_color_manual(values = c("olivedrab2", "maroon2")) +
  geom_hline(yintercept=0.25, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

mid_day9 %>% ggplot(aes(temperature_graph, deltaNPQ)) +
  geom_boxplot(size=0.5) +
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) +
  labs(x="temperature", y= "Day 9 Delta NPQ (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") +
  ylim(0, 0.75) + stat_mean() + 
  scale_color_manual(values = c("olivedrab3", "maroon4")) +
  geom_hline(yintercept=0.25, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for deltaNPQ
mid_day9 %>% group_by(treatment) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
mid_day9 %>% group_by(temp) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
#mid_day9 %>% group_by(treatment, rlc_day) %>% summarise_at(vars(pmax), list(mean = mean))

#Linear regression deltaNPQ vs Growth and Apices

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_growth_deltaNPQ_graph <- ggplot(mid_day9_sub, aes(x=deltaNPQ, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Delta NPQ vs 9-Day Growth (%)", x = " Delta NPQ (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
acan_growth_deltaNPQ_graph

#plot a regression between the photosynthetic independent variables of interest and apices/axis
acan_apices_deltaNPQ_graph <- ggplot(mid_day9_sub, aes(x=deltaNPQ, y=secondary_apices)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Delta NPQ vs # Secondary Apices", x = "Delta NPQ (rel. units)", 
       y = "Secondary Apices (#)") + stat_regline_equation(label.x = 0.25, label.y = 125) + stat_cor(label.x = 0.25, label.y = 128)
acan_apices_deltaNPQ_graph
