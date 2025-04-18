#This script is still very messy and was used only to do preliminary plotting of the data for growth of the wastewater phase I experiment
#By Angela Richards Donà
#Date: 1/16/25

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
waste_p1_growth <- read.csv("../data/wastewater/input/growth_2025/waste_phase1_growth_all.csv")

#assigns treatment as characters from integers then to factors
waste_p1_growth$treatment <- as.factor(as.character(waste_p1_growth$treatment))

#make a column that determines species
waste_p1_growth <- waste_p1_growth %>% 
  mutate(species = substr(ID, 1, 1))

#create subsets for the plots
#toggle between the plant_part for output. Use Day 9 for final analysis
ulva_g <- subset(waste_p1_growth, species == "u")
ulva_g$treatment_graph[ulva_g$treatment == 1] <- "1) 80 umol"
ulva_g$treatment_graph[ulva_g$treatment == 2] <- "2) 140 umol"
ulva_g$treatment_graph[ulva_g$treatment == 3] <- "3) 245 umol" 
ulva_g$treatment_graph[ulva_g$treatment == 4] <- "4) 428 umol"
ulva_g$treatment_graph[ulva_g$treatment == 5] <- "5) 748umol"
ulva_g$treatment_graph[ulva_g$treatment == 6] <- "6) 1308 umol"
ulva_g$treatment_graph[ulva_g$treatment == 7] <- "7) 2287 umol"
ulva_g$treatment_graph[ulva_g$treatment == 8] <- "8) 4000 umol"


hypnea_g <- subset(waste_p1_growth, species == "h")
hypnea_g$treatment_graph[hypnea_g$treatment == 1] <- "1) 80 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 2] <- "2) 140 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 3] <- "3)  245 umol" 
hypnea_g$treatment_graph[hypnea_g$treatment == 4] <- "4) 428 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 5] <- "5) 748 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 6] <- "6) 1308 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 7] <- "7) 2287 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 8] <- "8) 4000 umol"

#make histograms
canopy %>% ggplot(aes(d5_growth_percent)) +
  geom_histogram(binwidth=5, fill = "goldenrod1", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()
canopy %>% ggplot(aes(d9_growth_percent)) +
  geom_histogram(binwidth=5, fill = "goldenrod1", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run model without RLC_order as this has little effect
#mean_min28 used as random effect in lieu of run
growth_model_canopy <- lmer(formula = d9_growth_percent ~ salinity + nitrate +
                              (1 | plant_ID) + (1 | mean_mins28), data = canopy_g, REML = TRUE)

isSingular(growth_model_canopy)
hist(resid(growth_model_canopy))
plot(resid(growth_model_canopy) ~ fitted(growth_model_canopy))
qqnorm(resid(growth_model_canopy))
qqline(resid(growth_model_canopy))

#check the performance of the model for dataset: canopy
performance::check_model(growth_model_canopy)
r.squaredGLMM(growth_model_canopy)
summary(growth_model_canopy)
#view random effects levels
ranef(growth_model_canopy)
tab_model(growth_model_canopy, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_canopy))

#construct null model to perform likelihood ratio test REML must be FALSE
canopy_nitrate_null <- lmer(formula = d9_growth_percent ~ salinity + 
                              (1 | plant_ID) + (1 | mean_mins28), data = canopy_g, REML = FALSE)
canopy_model2 <- lmer(formula = d9_growth_percent ~ nitrate + salinity +
                        (1 | plant_ID) + (1 | mean_mins28), data = canopy_g, REML = FALSE)
anova(canopy_nitrate_null, canopy_model2)
canopy_salinity_null <- lmer(formula = d9_growth_percent ~ nitrate + (1 | plant_ID) + (1 | mean_mins28), data = canopy_g, REML = FALSE)
canopy_model3 <- lmer(formula = d9_growth_percent ~ nitrate + salinity + (1 | plant_ID) + (1 | mean_mins28), data = canopy_g, REML = FALSE)
anova(canopy_salinity_null, canopy_model3)

#plot canopy growth
ulva_p1_growth_plot <- ulva_g %>% 
  ggplot(aes(treatment_graph, growth)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "A - Growth", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("80 μmol N", "140 μmol N", "245 μmol N", "428 μmol N", 
                              "748 μmol N", "1308 μmol N", "2287 μmol N", "4000 μmol N")) + 
  ylim(0, 90) + stat_mean() + 
  #scale_color_manual(values = c("chartreuse4", "darkolivegreen2")) +
  geom_hline(yintercept=34, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position.inside = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
ulva_p1_growth_plot
ggsave("ulva_p1_growth_plot.png", path = "../analyze/wastewater/plots/")

#summarize the means for canopy
waste_p1_growth %>% group_by(species, treatment) %>% summarise_at(vars(growth), list(mean = mean))
waste_p1_growth %>% group_by(species) %>% summarise_at(vars(growth), list(mean = mean))



#UNDERSTORY

#-----------make a histogram of the growth rate data
hypnea_g %>% ggplot(aes(d5_growth_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()
hypnea_g %>% ggplot(aes(d9_growth_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run model without RLC_order as this has little effect
#mean_min28 used as random effect in lieu of run
growth_model_under <- lmer(formula = d9_growth_percent ~ salinity + nitrate +
                             (1 | plant_ID) + (1 | mean_mins28), data = hypnea_g, REML = TRUE)

isSingular(growth_model_under)
hist(resid(growth_model_under))
plot(resid(growth_model_under) ~ fitted(growth_model_under))
qqnorm(resid(growth_model_under))
qqline(resid(growth_model_under))

#check the performance of the model for dataset: under
performance::check_model(growth_model_under)
r.squaredGLMM(growth_model_under)
summary(growth_model_under)
#view random effects levels
ranef(growth_model_under)
tab_model(growth_model_under, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_under))

#construct null model to perform likelihood ratio test REML must be FALSE
under_nitrate_null <- lmer(formula = d9_growth_percent ~ salinity + 
                             (1 | plant_ID) + (1 | mean_mins28), data = hypnea_g, REML = FALSE)
under_model2 <- lmer(formula = d9_growth_percent ~ nitrate + salinity +
                       (1 | plant_ID) + (1 | mean_mins28), data = hypnea_g, REML = FALSE)
anova(under_nitrate_null, under_model2)
under_salinity_null <- lmer(formula = d9_growth_percent ~ nitrate + (1 | plant_ID) + (1 | mean_mins28), data = hypnea_g, REML = FALSE)
under_model3 <- lmer(formula = d9_growth_percent ~ nitrate + salinity + (1 | plant_ID) + (1 | mean_mins28), data = hypnea_g, REML = FALSE)
anova(under_salinity_null, under_model3)

#plot under growth
hypnea_p1_growth_plot <- hypnea_g %>% 
  ggplot(aes(treatment_graph, growth)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "B - Growth", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("80 μmol N", "140 μmol N", "245 μmol N", "428 μmol N", 
                              "748 μmol N", "1308 μmol N", "2287 μmol N", "4000 μmol N")) + 
  ylim(0, 90) + stat_mean() + 
  #scale_color_manual(values = c("maroon", "maroon2")) +
  geom_hline(yintercept=26, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position.inside = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
hypnea_p1_growth_plot
ggsave("hypnea_p1_growth_plot.png", path = "../analyze/wastewater/plots/")

#summarize the means for under
waste_p1_growth %>% group_by(nitrate, plant_part) %>% summarise_at(vars(d9_growth_percent), list(mean = mean))
waste_p1_growth %>% group_by(salinity, plant_part) %>% summarise_at(vars(d9_growth_percent), list(mean = mean))







