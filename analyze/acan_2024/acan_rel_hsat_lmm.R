#This script is similar to the ek_irrad_model but uses relative Hsat for analysis


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
library(mmtable2)
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



#for plots-match up number of treatment with detailed
ek_irrad_data$treatment_graph[ek_irrad_data$treatment == 0] <- "1) 35ppt/0.5umol"
ek_irrad_data$treatment_graph[ek_irrad_data$treatment == 1] <- "2) 28ppt/14umol" 

#for plots - same but for temp
ek_irrad_data$temperature_graph[ek_irrad_data$temperature == 24] <- "1) 24 °C"
ek_irrad_data$temperature_graph[ek_irrad_data$temperature == 28] <- "2) 28 °C" 

#make a histogram
hist(ek_irrad_data$supersat_rel, main = paste("Acanthophora spicifera"), col = "olivedrab3", labels = TRUE)

ek_irrad_data %>% ggplot(aes(supersat_rel)) +
        geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()
#run model without interaction between the treatments and temperature
#take RLC.Order out of the random effects because causing problems of singularity. R2 is same with or without (+ (1 | RLC.Order))
acan_rel_hsat_model <- lmer(formula = supersat_rel ~ treatment + temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data)

#residual plots of the data
hist(resid(acan_rel_hsat_model))
plot(resid(acan_rel_hsat_model) ~ fitted(acan_rel_hsat_model))
qqnorm(resid(acan_rel_hsat_model))
qqline(resid(acan_rel_hsat_model))


#check the performance of the model
performance::check_model(acan_rel_hsat_model)
r.squaredGLMM(acan_rel_hsat_model)
summary(acan_rel_hsat_model)
plot(allEffects(acan_rel_hsat_model))
tab_model(acan_rel_hsat_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

#construct null model to perform likelihood ratio test REML must be FALSE
acan_relhsat_treatment_null <- lmer(formula = supersat_rel ~ temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
acan_relhsat_model2 <- lmer(formula = supersat_rel ~ treatment + temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
anova(acan_relhsat_treatment_null, acan_relhsat_model2)
acan_relhsat_temperature_null <- lmer(formula = supersat_rel ~ treatment + (1 | run) + (1 | plantID) + (1 | rlc_order), 
                                      data = ek_irrad_data, REML = FALSE,
                                      control = lmerControl(optimizer = "bobyqa"))
acan_relhsat_model3 <- lmer(formula = supersat_rel ~ treatment + temperature + (1 | run) + (1 | plantID) + (1 | rlc_order), data = ek_irrad_data, REML = FALSE)
anova(acan_relhsat_temperature_null, acan_relhsat_model3)

#plots
ek_irrad_data %>% ggplot(aes(treatment_graph, supersat_rel)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.75, size = 3, aes(color = temperature), position = "jitter", show.legend = TRUE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Relative Hsat (%)", title= "A", subtitle = "Acanthophora spicifera") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "28ppt/14umolN")) + 
        ylim(10, 110) + stat_mean() + 
        geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#295102", "#7CB950")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", size = 14, vjust = -18, hjust = 0.05))

ek_irrad_data %>% ggplot(aes(temperature_graph, supersat_rel)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="Temperature (°C)", y= "Relative Hsat (%)", title= "B", subtitle = "Acanthophora spicifera") + 
  scale_x_discrete(labels = c("24°C", "28°C")) + 
  ylim(10, 110) + stat_mean() + 
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("#7CB950", "#BDE269")) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -18, hjust = 0.05))

#summarize the means for relative Hsat
ek_irrad_data %>% group_by(treatment) %>% summarise_at(vars(supersat_rel), list(mean = mean))
ek_irrad_data %>% group_by(temperature) %>% summarise_at(vars(supersat_rel), list(mean = mean))
ek_irrad_data %>% group_by(run) %>% summarise_at(vars(day_length_avg), list(mean = mean))

#add growth rate from other dataset to this one and subset by species for regression
growth_rate <- read.csv("../../../data/limu/acan_2024/input/acanthophora_growth_r1.csv")
growth_rate$treatment <- as.factor(growth_rate$treatment)
growth_rate$growth_rate_percent <- (growth_rate$final_weight - growth_rate$initial_weight) / growth_rate$initial_weight * 100
ek_irrad_data$growth_rate_percent <- growth_rate$growth_rate_percent
ek_irrad_data$apices <- growth_rate$mean_apices_per_axis

#plot a regression between the photosynthetic independent variables of interest and growth rate
acan_apices_relhsat_graph <- ggplot(ek_irrad_data, aes(x=supersat_rel, y=apices)) + 
  geom_point(alpha = 0.5, size = 1, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera relHsat vs Mean # Apices per Axis", x = "relHsat (%)", 
       y = "mean apices/axis") + stat_regline_equation(label.x = 25, label.y = 45) + stat_cor()
acan_apices_relhsat_graph

acan_growth_relhsat_graph <- ggplot(ek_irrad_data, aes(x=supersat_rel, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 1, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Acanthophora spicifera relHsat vs Growth Rate", x = "relHsat (%)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 30) + stat_cor()
acan_growth_relhsat_graph




#NO LONGER IN USE________________________________________________________________________
#check for equal variance
bartlett.test(supersat_rel ~ Treatment, data = ulva)
#run Welch's ANOVA if not equal variance
welch_anova_treatment_ulva <- oneway.test(supersat_rel ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment_ulva
welch_anova_temp_ulva <- oneway.test(supersat_rel ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp_ulva
games_howell_test(ulva, supersat_rel ~ Treatment, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(ek_irrad_model_noint, type = c("III"), ddf = "Satterthwaite")
#ulva_ek_irrad_model_aov <- aov(supersat_total ~ Treatment + Temperature, data = ulva)
#TukeyHSD(ulva_ek_irrad_model_aov, "Treatment", ordered = FALSE)
#TukeyHSD(ulva_ek_irrad_model_aov, "Temperature", ordered = FALSE)

#check for equal variance
bartlett.test(supersat_rel ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(supersat_rel ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(supersat_rel ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, supersat_rel ~ Treatment, conf.level = 0.95, detailed = TRUE)
games_howell_test(hypnea, supersat_rel ~ Temperature, conf.level = 0.95, detailed = TRUE)