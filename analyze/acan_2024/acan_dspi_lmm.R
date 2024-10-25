#script for mixed model analysis of Diurnal Saturated Photosynthesis Index


#load the various libraries
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
#library(mmtable2)
library(gt)
library(purrr)
library(stringr)
library(tidyr)

ek_irrad_data <- read.csv("../../../data/limu/acan_2024/transformed/acan_hsat_dspi.csv")

# assign run as a factor
ek_irrad_data$run <- as.factor(ek_irrad_data$run)

#assign temperature as a factor
ek_irrad_data$temperature <- as.factor(ek_irrad_data$temperature)

#assigns treatment as characters from integers then to factors
ek_irrad_data$treatment <- as.factor(as.character(ek_irrad_data$treatment))

ek_irrad_data$rlc_order <- as.factor(ek_irrad_data$rlc_order)

#convert dspi_pmax_total from micromol to mol
ek_irrad_data$dspi_mol <- ek_irrad_data$dspi_pmax_total / 1e+6

#for plots-match up number of treatment with detailed
ek_irrad_data$treatment_graph[ek_irrad_data$treatment == 0] <- "1) 35ppt/0.5umol"
ek_irrad_data$treatment_graph[ek_irrad_data$treatment == 1] <- "2) 28ppt/14umol" 

#for plots - same but for temp
ek_irrad_data$temperature_graph[ek_irrad_data$temperature == 24] <- "1) 24 °C"
ek_irrad_data$temperature_graph[ek_irrad_data$temperature == 28] <- "2) 28 °C" 

#make a histogram and residual plots of the data
hist(ek_irrad_data$dspi_mol, main = paste("Acanthophora spicifera"), col = "brown", labels = TRUE)
ek_irrad_data %>% ggplot(aes(dspi_mol)) +
        geom_histogram(binwidth=0.15, fill = "brown3", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

#run model
dspi_model_acan <- lmer(formula = dspi_mol ~ treatment + temperature + (1 |run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data)

#plot residuals for Pmax * Hsat
hist(resid(dspi_model_acan))
plot(resid(dspi_model_acan) ~ fitted(dspi_model_acan))
qqnorm(resid(dspi_model_acan))
qqline(resid(dspi_model_acan))

#check the performance of the model with Pmax
performance::check_model(dspi_model_acan)
r.squaredGLMM(dspi_model_acan)
summary(dspi_model_acan)
plot(allEffects(dspi_model_acan))
tab_model(dspi_model_acan, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

#Construct null model to perform likelihood ratio test REML must be FALSE
acan_dspi_treatment_null <- lmer(formula = dspi_mol ~ temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
acan_dspi_model2 <- lmer(formula = dspi_mol ~ treatment + temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
anova(acan_dspi_treatment_null, acan_dspi_model2)
acan_dspi_temperature_null <- lmer(formula = dspi_mol ~ treatment + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
acan_dspi_model3 <- lmer(formula = dspi_mol ~ treatment + temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
anova(acan_dspi_temperature_null, acan_dspi_model3)

#plot DSPI
ek_irrad_data %>% ggplot(aes(treatment_graph, dspi_mol)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.75, size = 3, aes(color = temperature), position = "jitter", show.legend = TRUE) + 
        labs(x="Treatment", y= "DSPI (mol m-2)", title= "C", subtitle = "Acanthophora spicifera") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "28ppt/14umolN")) + 
        ylim(0, 10) + stat_mean() + 
        geom_hline(yintercept=5.0, color = "red", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#494949", "#898989")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

ek_irrad_data %>% ggplot(aes(temperature_graph, dspi_mol)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="Temperature", y= "DSPI (mol m-2)", title= "D", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("24°C", "28°C")) + 
  ylim(0, 10) + stat_mean() + 
  geom_hline(yintercept=5.0, color = "red", size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("#292929", "#999999")) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for DSPI using Pmax
ek_irrad_data %>% group_by(treatment) %>% summarise_at(vars(dspi_mol), list(mean = mean))
ek_irrad_data %>% group_by(temperature) %>% summarise_at(vars(dspi_mol), list(mean = mean))
ek_irrad_data %>% group_by(run) %>% summarise_at(vars(day_length_avg), list(mean = mean))

#add growth rate from other dataset to this one and subset by species for regression
growth_rate <- read.csv("../../../data/limu/acan_2024/input/acanthophora_growth_r1.csv")
growth_rate$treatment <- as.factor(growth_rate$treatment)
growth_rate$growth_rate_percent <- (growth_rate$final_weight - growth_rate$initial_weight) / growth_rate$initial_weight * 100
ek_irrad_data$growth_rate_percent <- growth_rate$growth_rate_percent
ek_irrad_data$apices <- growth_rate$mean_apices_per_axis

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_apices_dspi_graph <- ggplot(ek_irrad_data, aes(x=dspi_mol, y=apices)) + 
  geom_point(alpha = 0.5, size = 1, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera DSPI vs Mean # Apices per Axis", x = "DSPI (mol)", 
       y = "mean apices/axis") + stat_regline_equation(label.x = 1, label.y = 47) + stat_cor()
acan_apices_dspi_graph

acan_growth_dspi_graph <- ggplot(ek_irrad_data, aes(x=dspi_mol, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 1, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera DSPI vs Growth Rate", x = "DSPI (mol)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 1, label.y = 30) + stat_cor()
acan_growth_dspi_graph
