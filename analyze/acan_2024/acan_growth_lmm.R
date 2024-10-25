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
acan_growth <- read.csv("../../../data/limu/acan_2024/input/acanthophora_growth_r1.csv")

#make a new column for weight change (difference final from initial)
acan_growth$growth_rate_percent <- (acan_growth$final_weight - acan_growth$initial_weight) / acan_growth$initial_weight * 100

#make a new column for daily growth rate from 8 day study (steady growth rate assumed rather than exponential)
acan_growth$steady_growth_daily <- acan_growth$growth_rate_percent / 8

#make a new column that keeps only the numerical values (effectively removes the C in temperatures for consistency)
#acan_growth$temp_clean <- as.factor(substr(acan_growth$temperature, 1, 2))

#assigns temperature as a factor
acan_growth$temperature <- as.factor(acan_growth$temperature)

#assigns treatment as characters from integers then to factors
acan_growth$treatment <- as.factor(as.character(acan_growth$treatment))

#assign run as a factor
acan_growth$run <- as.factor(acan_growth$run)

#assign plant ID as a factor
acan_growth$plant_ID <- as.factor(acan_growth$plant_ID)

#assign RLC order as a factor
acan_growth$rlc_order <- as.factor(acan_growth$rlc_order)


#create subsets for the plots
acan_growth <- subset(acan_growth, secondary_apices != "200")

acan_growth$treatment_graph[acan_growth$treatment == 0] <- "1) 35ppt/0.5umol"
acan_growth$treatment_graph[acan_growth$treatment == 1] <- "2) 28ppt/14umol" 

acan_growth$temperature_graph[acan_growth$temperature == 24] <- "24°C"
acan_growth$temperature_graph[acan_growth$temperature == 28] <- "28°C"



#-----------make a histogram of the growth rate data
hist(acan_growth$growth_rate_percent, main = paste("Acanthophora spicifera 9-Day Growth (%)"), col = "maroon", labels = TRUE)
#or
acan_growth %>% ggplot(aes(growth_rate_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run model without RLC_order as this has little effect
#add number of initial axes to the random effects
growth_model_acan <- lmer(formula = growth_rate_percent ~ treatment + temperature +
                            (1 | plant_ID) + (1 | run) + (1 | initial_axes), data = acan_growth, REML = TRUE)

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
acan_growth_treatment_null <- lmer(formula = growth_rate_percent ~ temperature + (1 | plant_ID) + (1 | run), data = acan_growth, REML = FALSE)
acan_growth_model2 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant_ID) + (1 | run), data = acan_growth, REML = FALSE)
anova(acan_growth_treatment_null, acan_growth_model2)
acan_growth_temperature_null <- lmer(formula = growth_rate_percent ~ treatment + (1 | plant_ID) + (1 | run), data = acan_growth, REML = FALSE)
acan_growth_model3 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant_ID) + (1 | run), data = acan_growth, REML = FALSE)
anova(acan_growth_temperature_null, acan_growth_model3)

#plots
acan_growth %>% ggplot(aes(treatment_graph, growth_rate_percent)) + 
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
acan_growth %>% group_by(treatment) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))
acan_growth %>% group_by(temperature) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))

#plot temperature
acan_growth %>% ggplot(aes(temperature_graph, growth_rate_percent)) + 
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
acan_growth %>% ggplot(aes(secondary_apices)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()
#run model without RLC_order as this has little effect
apices_model_acan <- lmer(formula = secondary_apices ~ treatment + temperature +
                            (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = acan_growth, REML = TRUE)

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
acan_apices_treatment_null <- lmer(formula = secondary_apices~ temperature + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = acan_growth, REML = FALSE)
acan_apices_model2 <- lmer(formula = secondary_apices ~ treatment + temperature + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = acan_growth, REML = FALSE)
anova(acan_apices_treatment_null, acan_apices_model2)
acan_apices_temperature_null <- lmer(formula = secondary_apices ~ treatment + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = acan_growth, REML = FALSE)
acan_apices_model3 <- lmer(formula = secondary_apices ~ treatment + temperature + (1 | plant_ID) + (1 | run) + (1 | rlc_order), data = acan_growth, REML = FALSE)
anova(acan_apices_temperature_null, acan_apices_model3)

#plots
acan_growth %>% ggplot(aes(treatment_graph, secondary_apices)) + 
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
acan_growth %>% group_by(treatment) %>% summarise_at(vars(secondary_apices), list(mean = mean))
acan_growth %>% group_by(temperature) %>% summarise_at(vars(secondary_apices), list(mean = mean))


#plot temperature
acan_growth %>% ggplot(aes(temperature_graph, secondary_apices)) + 
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

