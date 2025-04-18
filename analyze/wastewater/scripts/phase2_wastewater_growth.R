#This script is for growth of wastewater treated ulva and hypnea
#By Angela Richards Donà
#Date: 04/08/25

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
waste_p2_growth <- read.csv("../data/wastewater/input/growth_2025/waste_p2_growth.csv")

#assigns treatment as characters from integers then to factors
waste_p2_growth$treatment <- as.factor(as.character(waste_p2_growth$treatment))
waste_p2_growth$salinity <- as.factor(as.character(waste_p2_growth$salinity))


#create subsets for the plots
#toggle between the plant_part for output. Use Day 9 for final analysis
ulva_g <- subset(waste_p2_growth, species == "u")
ulva_g$treatment_graph[ulva_g$treatment == 1] <- "1) 80 umol"
ulva_g$treatment_graph[ulva_g$treatment == 3] <- "2) 245 umol" 
ulva_g$treatment_graph[ulva_g$treatment == 5] <- "3) 748umol"
ulva_g$treatment_graph[ulva_g$treatment == 7] <- "4) 2287 umol"


hypnea_g <- subset(waste_p2_growth, species == "h")
hypnea_g$treatment_graph[hypnea_g$treatment == 1] <- "1) 80 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 3] <- "2)  245 umol" 
hypnea_g$treatment_graph[hypnea_g$treatment == 5] <- "3) 748 umol"
hypnea_g$treatment_graph[hypnea_g$treatment == 7] <- "4) 2287 umol"

#make histograms
ulva_g %>% ggplot(aes(growth_d9)) +
  geom_histogram(binwidth=5, fill = "darkolivegreen3", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run the model (removed lunar_phase since this had no effect and caused issues of Singularity)
ulva_growth_model <- lmer(formula = growth_d9 ~ salinity + treatment +
                              (1 | plant_id) + (1 | run), data = ulva_g, REML = TRUE)

#Check performance
isSingular(ulva_growth_model)
hist(resid(ulva_growth_model))
plot(resid(ulva_growth_model) ~ fitted(ulva_growth_model))
qqnorm(resid(ulva_growth_model))
qqline(resid(ulva_growth_model))

#check the performance of the model for dataset: ulva
performance::check_model(ulva_growth_model)
r.squaredGLMM(ulva_growth_model)
summary(ulva_growth_model)

#view random effects levels
tab_model(ulva_growth_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(ulva_growth_model))

#construct null model to perform likelihood ratio test REML must be FALSE
#treatment = nitrate/phosphate in this case
ulva_treatment_null <- lmer(formula = growth_d9 ~ salinity + 
                              (1 | plant_id) + (1 | run) + (1 | lunar_phase), data = ulva_g, REML = FALSE)
ulva_with_treatment <- lmer(formula = growth_d9 ~ treatment + salinity +
                        (1 | plant_id) + (1 | run) + (1 | lunar_phase), data = ulva_g, REML = FALSE)
anova(ulva_nitrate_null, ulva_with_treatment)
ulva_salinity_null <- lmer(formula = growth_d9 ~ treatment + (1 | plant_id) + (1 | run) + (1 | lunar_phase), data = ulva_g, REML = FALSE)
ulva_with_salinity <- lmer(formula = growth_d9 ~ treatment + salinity + (1 | plant_id) + (1 | run) + (1 | lunar_phase), data = ulva_g, REML = FALSE)
anova(ulva_salinity_null, ulva_with_salinity)

#plot ulva growth
ulva_p2_growth_plot <- ulva_g %>% 
  ggplot(aes(treatment_graph, growth_d9)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = salinity), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "A - Growth", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("80 μmol N", "245 μmol N", 
                              "748 μmol N", "2287 μmol N")) + 
  ylim(-30, 370) + stat_mean() + 
  scale_color_manual(values = c("18" = "chartreuse4","22" = "darkolivegreen2")) +
  geom_hline(yintercept=125, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.95,0.95), legend.justification = c(1,1), legend.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
ulva_p2_growth_plot
ggsave("ulva_p2_growth_plot.png", path = "../analyze/wastewater/plots/")

#plot for salinity
ulva_p2_growth_plot2 <- ulva_g %>% 
  ggplot(aes(salinity, growth_d9)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "B - Growth", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("18 ppt", "22 ppt")) + 
  ylim(-30, 370) + stat_mean() + 
  #scale_color_manual(values = c("deeppink", "maroon4", "darkred", "firebrick1")) +
  geom_hline(yintercept=64, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.95,0.95), legend.justification = c(1,1), legend.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
ulva_p2_growth_plot2
ggsave("ulva_p2_growth_plot.png", path = "../analyze/wastewater/plots/")

#summarize the means for ulva
waste_p2_growth %>% group_by(species, treatment) %>% summarise_at(vars(growth_d9), list(mean = mean))
waste_p2_growth %>% group_by(species) %>% summarise_at(vars(growth_d9), list(mean = mean))



#Hypnea

#-----------make a histogram of the growth rate data
hypnea_g %>% ggplot(aes(growth_d9)) +
  geom_histogram(binwidth=5, fill = "deeppink4", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run model without RLC_order as this has little effect
#mean_min28 used as random effect in lieu of run
growth_model_hypnea <- lmer(formula = growth_d9 ~ salinity + treatment +
                             (1 | plant_id) + (1 | run), data = hypnea_g, REML = TRUE)

isSingular(growth_model_hypnea)
hist(resid(growth_model_hypnea))
plot(resid(growth_model_hypnea) ~ fitted(growth_model_hypnea))
qqnorm(resid(growth_model_hypnea))
qqline(resid(growth_model_hypnea))

#check the performance of the model for dataset: hypnea
performance::check_model(growth_model_hypnea)
r.squaredGLMM(growth_model_hypnea)
summary(growth_model_hypnea)
#view random effects levels
#ranef(growth_model_hypnea)
tab_model(growth_model_hypnea, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_hypnea))

#construct null model to perform likelihood ratio test REML must be FALSE
hypnea_nitrate_null <- lmer(formula = growth_d9 ~ salinity + 
                             (1 | plant_id) + (1 | run), data = hypnea_g, REML = FALSE)
hypnea_model2 <- lmer(formula = growth_d9 ~ treatment + salinity +
                       (1 | plant_id) + (1 | run), data = hypnea_g, REML = FALSE)
anova(hypnea_nitrate_null, hypnea_model2)
hypnea_salinity_null <- lmer(formula = growth_d9 ~ treatment + (1 | plant_id) + (1 | run), data = hypnea_g, REML = FALSE)
hypnea_model3 <- lmer(formula = growth_d9 ~ treatment + salinity + (1 | plant_id) + (1 | run), data = hypnea_g, REML = FALSE)
anova(hypnea_salinity_null, hypnea_model3)

#plot hypnea growth
hypnea_p2_growth_plot <- hypnea_g %>% 
  ggplot(aes(treatment_graph, growth_d9)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = salinity), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "B - Growth", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("80 μmol N", "245 μmol N",
                              "748 μmol N", "2287 μmol N")) + 
  ylim(-30, 370) + stat_mean() + 
  scale_color_manual(values = c("deeppink", "maroon4")) +
  geom_hline(yintercept=64, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.95,0.95), legend.justification = c(1,1), legend.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
hypnea_p2_growth_plot
ggsave("hypnea_p2_growth_plot.png", path = "../analyze/wastewater/plots/")

#plot for salinity
hypnea_p2_growth_plot2 <- hypnea_g %>% 
  ggplot(aes(salinity, growth_d9)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "B - Growth", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("18 ppt", "22 ppt")) + 
  ylim(-0, 160) + stat_mean() + 
  #scale_color_manual(values = c("deeppink", "maroon4", "darkred", "firebrick1")) +
  geom_hline(yintercept=64, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.95,0.95), legend.justification = c(1,1), legend.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
hypnea_p2_growth_plot2
ggsave("hypnea_p2_growth_plot2.png", path = "../analyze/wastewater/plots/")






