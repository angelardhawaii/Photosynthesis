#script to run model for Ek -- similar to rETRmax and growth_rate scripts
#ek.est = Ek estimate - terminology from WebFitt model that produced the values on Phytotools package

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


#load this file for normalized to quantum efficiency of photosynthesis per Silsbe and Kromkamp
acan_ps <- read.csv("../../../data/limu/acan_2024/transformed/acan_ek_alpha_normalized_2024.csv")

# assign run as a factor
acan_ps$run <- as.factor(acan_ps$run)

#assign temperature as a factor
acan_ps$temp <- as.factor(acan_ps$temp)

#assigns treatment as characters from integers then to factors
acan_ps$treatment <- as.factor(as.character(acan_ps$treatment))


#toggle between the species for output. Use Day 9 for final analysis
acan_day9 <- subset(acan_ps, rlc_day == 9)
acan_day9$treatment_graph[acan_day9$treatment == 0] <- "1) 35ppt/0.5umol"
acan_day9$treatment_graph[acan_day9$treatment == 1] <- "2) 28ppt/14umol" 

acan_day9$temperature_graph[acan_day9$temp == 24] <- "24°C"
acan_day9$temperature_graph[acan_day9$temp == 28] <- "28°C"


acan_day9 %>% ggplot(aes(ek.est)) +
  geom_histogram(binwidth=5, fill = "orange", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
acan_day9_ek.est_model <- lmer(formula = ek.est ~ treatment + temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9)

#construct null model to perform likelihood ratio test REML must be FALSE
acan_ek.est_treatment_null <- lmer(formula = ek.est ~ temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
acan_ek.est_model2 <- lmer(formula = ek.est ~ treatment + temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
anova(acan_ek.est_treatment_null, acan_ek.est_model2)
acan_ek.est_temperature_null <- lmer(formula = ek.est ~ treatment + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
acan_ek.est_model3 <- lmer(formula = ek.est ~ treatment + temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
anova(acan_ek.est_temperature_null, acan_ek.est_model3)

#make residual plots of the data for Acanthophora
hist(resid(acan_day9_ek.est_model))
plot(resid(acan_day9_ek.est_model) ~ fitted(acan_day9_ek.est_model))
qqnorm(resid(acan_day9_ek.est_model))
qqline(resid(acan_day9_ek.est_model))

#check the performance of the model
performance::check_model(acan_day9_ek.est_model)
rsquared(acan_day9_ek.est_model)

#make plots and tables for the data
tab_model(acan_day9_ek.est_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(acan_day9_ek.est_model))

acan_day9 %>% ggplot(aes(treatment_graph, ek.est)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = temp), position = "jitter", show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "Day 9 Ek (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "28ppt/14umolN")) + 
  ylim(-1, 150) + stat_mean() + 
  scale_color_manual(values = c("brown", "brown2")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

acan_day9 %>% ggplot(aes(temperature_graph, ek.est)) +
  geom_boxplot(size=0.5) +
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) +
  labs(x="temperature", y= "Day 9 Ek (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") +
  ylim(-1, 150) + stat_mean() + 
  scale_color_manual(values = c("brown4", "brown1")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for ek.est
acan_day9 %>% group_by(treatment) %>% summarise_at(vars(ek.est), list(mean = mean))
acan_day9 %>% group_by(temp) %>% summarise_at(vars(ek.est), list(mean = mean))

#Linear regression Ek vs Growth and Ek vs Apices/Axis--------------------------------------------------------

#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("../../../data/limu/acan_2024/input/acanthophora_growth_r1.csv")
growth_rate$treatment <- as.factor(growth_rate$treatment)

#make a new column for weight change (difference final from initial)
growth_rate$growth_rate_percent <- (growth_rate$final_weight - growth_rate$initial_weight) / growth_rate$initial_weight * 100
acan_day9$growth_rate <- growth_rate$growth_rate_percent
acan_day9$mean_apices_per_axis <- growth_rate$mean_apices_per_axis

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_growth_ek.est_graph <- ggplot(acan_day9, aes(x=ek.est, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
acan_growth_ek.est_graph

#plot a regression between the photosynthetic independent variables of interest and apices/axis
acan_apices_ek_graph <- ggplot(acan_day9, aes(x=ek.est, y=mean_apices_per_axis)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Ek vs Mean Apices/Axis", x = "Ek (μmols photons m-2 s-1)", 
       y = "Mean Apices per Axis") + stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
acan_apices_ek_graph

