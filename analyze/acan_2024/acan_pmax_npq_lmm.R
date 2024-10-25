#script to run model for photosynthesis data -- similar to growth_rate script

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

# assign deltaNPQ as a numeric
acan_ps$deltaNPQ <- as.numeric(acan_ps$deltaNPQ)


#toggle between the species for output. Use Day 9 for final analysis
acan_day9 <- subset(acan_ps, rlc_day == 9)
acan_day9$treatment_graph[acan_day9$treatment == 0] <- "1) 35ppt/0.5umol"
acan_day9$treatment_graph[acan_day9$treatment == 1] <- "2) 28ppt/14umol" 

acan_day9$temperature_graph[acan_day9$temp == 24] <- "24°C"
acan_day9$temperature_graph[acan_day9$temp == 28] <- "28°C"


#------------Pmax----------------
#plot histogram
acan_day9 %>% ggplot(aes(pmax)) +
  geom_histogram(binwidth=5, fill = "goldenrod1", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
acan_day9_pmax_model <- lmer(formula = pmax ~ treatment + temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9)

#construct null model to perform likelihood ratio test REML must be FALSE
acan_pmax_treatment_null <- lmer(formula = pmax ~ temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
acan_pmax_model2 <- lmer(formula = pmax ~ treatment + temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
anova(acan_pmax_treatment_null, acan_pmax_model2)
acan_pmax_temperature_null <- lmer(formula = pmax ~ treatment + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
acan_pmax_model3 <- lmer(formula = pmax ~ treatment + temp + (1 | run) + (1 | plantID) + (1 | rlc_order), data = acan_day9, REML = FALSE)
anova(acan_pmax_temperature_null, acan_pmax_model3)

#make residual plots of the data for Acanthophora
hist(resid(acan_day9_pmax_model))
plot(resid(acan_day9_pmax_model) ~ fitted(acan_day9_pmax_model))
qqnorm(resid(acan_day9_pmax_model))
qqline(resid(acan_day9_pmax_model))

#check the performance of the model
performance ::check_model(acan_day9_pmax_model)
r.squaredGLMM(acan_day9_pmax_model)

#make plots and tables for the data
tab_model(acan_day9_pmax_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(acan_day9_pmax_model))

acan_day9 %>% ggplot(aes(treatment_graph, pmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = temp), position = "jitter", show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "Day 9 Pmax (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "28ppt/14umolN")) + 
  ylim(-1, 90) + stat_mean() + 
  scale_color_manual(values = c("brown", "brown2")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

acan_day9 %>% ggplot(aes(temperature_graph, pmax)) +
  geom_boxplot(size=0.5) +
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) +
  labs(x="temperature", y= "Day 9 Pmax (μmols electrons m-2 s-1)", title= "A", subtitle = "Acanthophora spicifera") +
  ylim(-1, 90) + stat_mean() + 
  scale_color_manual(values = c("brown4", "brown1")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for pmax
acan_day9 %>% group_by(treatment) %>% summarise_at(vars(pmax), list(mean = mean))
acan_day9 %>% group_by(temp) %>% summarise_at(vars(pmax), list(mean = mean))
#acan_day9 %>% group_by(treatment, rlc_day) %>% summarise_at(vars(pmax), list(mean = mean))


#Linear regression Pmax vs Growth--------------------------------------------------------------------------------

#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("../../../data/limu/acan_2024/input/acanthophora_growth_r1.csv")
growth_rate$treatment <- as.factor(growth_rate$treatment)

#make a new column for weight change (difference final from initial)
growth_rate$growth_rate_percent <- (growth_rate$final_weight - growth_rate$initial_weight) / growth_rate$initial_weight * 100
acan_day9$growth_rate_percent <- growth_rate$growth_rate_percent
acan_day9$secondary_apices <- growth_rate$secondary_apices
acan_day9_sub <- subset(acan_day9, secondary_apices != "200")

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_growth_pmax_graph <- ggplot(acan_day9_sub, aes(x=pmax, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Pmax vs 9-Day Growth (%)", x = "Pmax (μmols electrons m-2 s-1)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
acan_growth_pmax_graph

#plot a regression between the photosynthetic independent variables of interest and apices/axis
acan_apices_pmax_graph <- ggplot(acan_day9_sub, aes(x=pmax, y=secondary_apices)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Pmax vs # Secondary Apices", x = "Pmax (μmols electrons m-2 s-1)", 
       y = "Secondary Apices (#)") + stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
acan_apices_pmax_graph





#--------------------NPQmax--------------------------

#plot histogram
acan_day9 %>% ggplot(aes(maxNPQ_Ypoint1)) +
  geom_histogram(binwidth=0.05, fill = "violet", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
acan_day9_maxNPQ_model <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + temp + (1 | plantID), data = acan_day9)

#construct null model to perform likelihood ratio test REML must be FALSE
acan_maxNPQ_treatment_null <- lmer(formula = maxNPQ_Ypoint1 ~ temp + (1 | plantID), data = acan_day9, REML = FALSE)
acan_maxNPQ_model2 <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + temp + (1 | plantID), data = acan_day9, REML = FALSE)
anova(acan_maxNPQ_treatment_null, acan_maxNPQ_model2)
acan_maxNPQ_temperature_null <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + (1 | plantID), data = acan_day9, REML = FALSE)
acan_maxNPQ_model3 <- lmer(formula = maxNPQ_Ypoint1 ~ treatment + temp + (1 | plantID), data = acan_day9, REML = FALSE)
anova(acan_maxNPQ_temperature_null, acan_maxNPQ_model3)

#make residual plots of the data for Acanthophora
hist(resid(acan_day9_maxNPQ_model))
plot(resid(acan_day9_maxNPQ_model) ~ fitted(acan_day9_maxNPQ_model))
qqnorm(resid(acan_day9_maxNPQ_model))
qqline(resid(acan_day9_maxNPQ_model))

#check the performance of the model
performance ::check_model(acan_day9_maxNPQ_model)
r.squaredGLMM(acan_day9_maxNPQ_model)

#make plots and tables for the data
tab_model(acan_day9_maxNPQ_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(acan_day9_maxNPQ_model))

acan_day9 %>% ggplot(aes(treatment_graph, maxNPQ_Ypoint1)) + 
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

acan_day9 %>% ggplot(aes(temperature_graph, maxNPQ_Ypoint1)) +
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
acan_day9 %>% group_by(treatment) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))
acan_day9 %>% group_by(temp) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))
#acan_day9 %>% group_by(treatment, rlc_day) %>% summarise_at(vars(pmax), list(mean = mean))

#Linear regression NPQ vs Growth and Apices

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_growth_NPQmax_graph <- ggplot(acan_day9_sub, aes(x=maxNPQ_Ypoint1, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera NPQmax vs 9-Day Growth (%)", x = " NPQmax (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
acan_growth_NPQmax_graph

#plot a regression between the photosynthetic independent variables of interest and apices/axis
acan_apices_NPQmax_graph <- ggplot(acan_day9_sub, aes(x=maxNPQ_Ypoint1, y=secondary_apices)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera NPQmax vs # Secondary Apices", x = "NPQmax (rel. units)", 
       y = "Secondary Apices (#)") + stat_regline_equation(label.x = 0.25, label.y = 125) + stat_cor(label.x = 0.25, label.y = 128)
acan_apices_NPQmax_graph

#-----------------delta NPQ------------------
#plot histogram
acan_day9 %>% ggplot(aes(deltaNPQ)) +
  geom_histogram(binwidth=.05, fill = "violet", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
acan_day9_deltaNPQ_model <- lmer(formula = deltaNPQ ~ treatment + temp + (1 | plantID), data = acan_day9)

#construct null model to perform likelihood ratio test REML must be FALSE
acan_deltaNPQ_treatment_null <- lmer(formula = deltaNPQ ~ temp + (1 | plantID), data = acan_day9, REML = FALSE)
acan_deltaNPQ_model2 <- lmer(formula = deltaNPQ ~ treatment + temp + (1 | plantID), data = acan_day9, REML = FALSE)
anova(acan_deltaNPQ_treatment_null, acan_deltaNPQ_model2)
acan_deltaNPQ_temperature_null <- lmer(formula = deltaNPQ ~ treatment + (1 | plantID), data = acan_day9, REML = FALSE)
acan_deltaNPQ_model3 <- lmer(formula = deltaNPQ ~ treatment + temp + (1 | plantID), data = acan_day9, REML = FALSE)
anova(acan_deltaNPQ_temperature_null, acan_deltaNPQ_model3)

#make residual plots of the data for Acanthophora
hist(resid(acan_day9_deltaNPQ_model))
plot(resid(acan_day9_deltaNPQ_model) ~ fitted(acan_day9_deltaNPQ_model))
qqnorm(resid(acan_day9_deltaNPQ_model))
qqline(resid(acan_day9_deltaNPQ_model))

#check the performance of the model
performance ::check_model(acan_day9_deltaNPQ_model)
r.squaredGLMM(acan_day9_deltaNPQ_model)

#make plots and tables for the data
tab_model(acan_day9_deltaNPQ_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(acan_day9_deltaNPQ_model))

acan_day9 %>% ggplot(aes(treatment_graph, deltaNPQ)) + 
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

acan_day9 %>% ggplot(aes(temperature_graph, deltaNPQ)) +
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
acan_day9 %>% group_by(treatment) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
acan_day9 %>% group_by(temp) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
#acan_day9 %>% group_by(treatment, rlc_day) %>% summarise_at(vars(pmax), list(mean = mean))

#Linear regression deltaNPQ vs Growth and Apices

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_growth_deltaNPQ_graph <- ggplot(acan_day9_sub, aes(x=deltaNPQ, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Delta NPQ vs 9-Day Growth (%)", x = " Delta NPQ (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
acan_growth_deltaNPQ_graph

#plot a regression between the photosynthetic independent variables of interest and apices/axis
acan_apices_deltaNPQ_graph <- ggplot(acan_day9_sub, aes(x=deltaNPQ, y=secondary_apices)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera Delta NPQ vs # Secondary Apices", x = "Delta NPQ (rel. units)", 
       y = "Secondary Apices (#)") + stat_regline_equation(label.x = 0.25, label.y = 125) + stat_cor(label.x = 0.25, label.y = 128)
acan_apices_deltaNPQ_graph
