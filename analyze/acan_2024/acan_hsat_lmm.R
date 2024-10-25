#script for mixed model analysis of Ek and irradiance (Hsat)
#Supersat is the time above Ek for each individual in minutes, these are binned into three time periods from day 1 to day 9
#Supersat avg is the mean of the three time periods

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
library(magrittr)


ek_irrad_data <- read.csv("../../../data/limu/acan_2024/transformed/acan_hsat_dspi.csv")

# assign run as a factor
ek_irrad_data$run <- as.factor(ek_irrad_data$run)

#assign temperature as a factor
ek_irrad_data$temperature <- as.factor(as.character(ek_irrad_data$temperature))

#assigns treatment as characters from integers then to factors
ek_irrad_data$treatment <- as.factor(as.character(ek_irrad_data$treatment))

#for plots-match up number of treatment with detailed
ek_irrad_data$treatment_graph[ek_irrad_data$treatment == 0] <- "1) 35ppt/0.5umol"
ek_irrad_data$treatment_graph[ek_irrad_data$treatment == 1] <- "2) 28ppt/14umol" 

#for plots - same but for temp
ek_irrad_data$temperature_graph[ek_irrad_data$temperature == 24] <- "1) 24 °C"
ek_irrad_data$temperature_graph[ek_irrad_data$temperature == 28] <- "2) 28 °C" 


#make a histogram of the data 
hist(ek_irrad_data$supersat_avg, main = paste("Acanthophora spicifera"), col = "maroon", labels = TRUE)

ek_irrad_data %>% ggplot(aes(supersat_avg)) +
        geom_histogram(binwidth=25, fill = "maroon", color = "black", linewidth = 0.25, alpha = 0.85) +
        theme_bw()

#run model without interaction between the treatments and temperature - supersat_avg is in minutes
hsat_model <- lmer(formula = supersat_avg ~ treatment + temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data)

#check residual plots
hist(resid(hsat_model))
plot(resid(hsat_model) ~ fitted(hsat_model))
qqnorm(resid(hsat_model))
qqline(resid(hsat_model))

#check the performance of the model
performance::check_model(hsat_model)
r.squaredGLMM(hsat_model)
summary(hsat_model)
plot(allEffects(hsat_model))
tab_model(hsat_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

#construct null model to perform likelihood ratio test REML must be FALSE
hsat_treatment_null <- lmer(formula = supersat_avg ~ temperature + (1 | run) + (1 |plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
hsat_model2 <- lmer(formula = supersat_avg ~ treatment + temperature + (1 | run) + (1 |plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
anova(hsat_treatment_null, hsat_model2)
hsat_temperature_null <- lmer(formula = supersat_avg ~ treatment + (1 | run) + (1 |plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
hsat_model3 <- lmer(formula = supersat_avg ~ treatment + temperature + (1 | run) + (1 |plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
anova(hsat_temperature_null, hsat_model3)

ek_irrad_data %>% ggplot(aes(treatment_graph, supersat_avg)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.75, size = 3, aes(color = temperature), position = "jitter", show.legend = TRUE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Hsat Time (min)", title= "A", subtitle = "Acanthophora spicifera") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "28ppt/14umolN")) + 
        ylim(0, 650) + stat_mean() + 
        geom_hline(yintercept=200, color = "red", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#9C0627", "#BB589F")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))


ek_irrad_data %>% ggplot(aes(temperature_graph, supersat_avg)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="Temperature (°C)", y= "Hsat Time (min)", title= "B", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("24 °C", "28 °C")) + 
  ylim(0, 650) + stat_mean() + 
  geom_hline(yintercept=200, color = "red", size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("#BB589F", "#F4B4E2")) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for variables of interest
ek_irrad_data %>% group_by(treatment) %>% summarise_at(vars(supersat_avg), list(mean = mean))
ek_irrad_data %>% group_by(temperature) %>% summarise_at(vars(supersat_avg), list(mean = mean))
ek_irrad_data %>% group_by(run) %>% summarise_at(vars(day_length_avg), list(mean = mean))

#add growth rate from other dataset to this one and subset by species for regression
growth_rate <- read.csv("../../../data/limu/acan_2024/input/acanthophora_growth_r1.csv")
growth_rate$treatment <- as.factor(growth_rate$treatment)
growth_rate$growth_rate_percent <- (growth_rate$final_weight - growth_rate$initial_weight) / growth_rate$initial_weight * 100
ek_irrad_data$growth_rate_percent <- growth_rate$growth_rate_percent
ek_irrad_data$apices <- growth_rate$mean_apices_per_axis

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_apices_hsat_graph <- ggplot(ek_irrad_data, aes(x=supersat_avg, y=apices)) + 
        geom_point(alpha = 0.5, size = 1, show.legend = TRUE, aes(color = treatment)) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "Acanthophora spicifera Hsat vs Mean # Apices per Axis", x = "Hsat time (minutes)", 
             y = "mean apices/axis") + stat_regline_equation(label.x = 100, label.y = 30) + stat_cor()
acan_apices_hsat_graph

acan_growth_hsat_graph <- ggplot(ek_irrad_data, aes(x=supersat_avg, y=growth_rate_percent)) + 
        geom_point(alpha = 0.5, size = 1, show.legend = TRUE, aes(color = treatment)) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "Acanthophora spicifera Hsat vs Growth Rate", x = "Hsat time (minutes)", 
             y = "growth rate (%)") + stat_regline_equation(label.x = 100, label.y = 30) + stat_cor()
acan_growth_hsat_graph

