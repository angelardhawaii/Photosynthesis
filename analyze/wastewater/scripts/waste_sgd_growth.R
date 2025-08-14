#Run analysis on the combined dataset using 2021, 2022, 2023, and 2025 data from
#Three separate experimental periods
#2021-2022 SGD, 2023 Continuation of SGD, 2025 Wastewater Scenarios
#Subsetted for analysis of the individual species
#By Angela Richards Dona
#August 6, 2025

#load libraries
#libraries for tidying up
library(gt)
library(purrr)
library(stringr)
library(tidyr)
library(hms)
library(lubridate)
library(ggeffects)

# libraries for analysis
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

#libraries for plots and tables
library(ggplot2)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(lubridate)

#read in combined dataset
combined_growth_all <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/combined_data_all.csv")

#change variables to factors
combined_growth_all <- combined_growth_all %>%
  mutate(date = ymd(date)) %>%
  mutate(year = as.factor(year)) %>%
  mutate(salinity = as.factor(salinity)) %>%
  mutate(temp = as.factor(temp)) %>%
  mutate(nitrate = as.factor(nitrate)) %>%
  mutate(run_combo = as.factor(run_combo)) %>%
  mutate(plant_id = as.factor(plant_id))

glimpse(combined_growth_all)
  
#make a list to be used for plots that orders chronologically first then by nitrate value low to high
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "a"] <- "a) 0.5 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "b"] <- "c) 14 μmol" 
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "c"] <- "d) 27 μmol" 
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "d"] <- "f) 53 μmol" 
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "e"] <- "g) 53 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "f"] <- "i) 80 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "g"] <- "b) 0.5 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "h"] <- "e) 37 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "i"] <- "h) 53 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "j"] <- "j) 86 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "k"] <- "k) 80 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "l"] <- "l) 245 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "m"] <- "m) 748 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "n"] <- "n) 2287 μmol"

combined_growth_all$salinity_graph[combined_growth_all$salinity == "35"] <- "a) 35 ppt"
combined_growth_all$salinity_graph[combined_growth_all$salinity == "28"] <- "b) 28 ppt"
combined_growth_all$salinity_graph[combined_growth_all$salinity == "27"] <- "c) 27 ppt"
combined_growth_all$salinity_graph[combined_growth_all$salinity == "22"] <- "d) 22 ppt" 
combined_growth_all$salinity_graph[combined_growth_all$salinity == "21"] <- "e) 21 ppt" 
combined_growth_all$salinity_graph[combined_growth_all$salinity == "18"] <- "f) 18 ppt"
combined_growth_all$salinity_graph[combined_growth_all$salinity == "15"] <- "g) 15 ppt"
combined_growth_all$salinity_graph[combined_growth_all$salinity == "11"] <- "h) 11 ppt"

#Subset by species
combined_growth_u <- combined_growth_all %>%
  filter(species == "u")

combined_growth_h <- combined_growth_all %>%
  filter(species == "h")

xtabs(~ salinity + temp, data = combined_growth_u)
#Run lmm model for Ulva_________________________________________________________
growth_model_u <- lmer(formula = growth_d9 ~ nitrate + temp + (1 | plant_id) + 
                         (1 | illumination) + (1 | run_combo), 
                       data = combined_growth_all, REML = TRUE)

hist(resid(growth_model_u))
plot(resid(growth_model_u) ~ fitted(growth_model_u))
qqnorm(resid(growth_model_u))
qqline(resid(growth_model_u))

#check the performance of the model for dataset: u
performance::check_model(growth_model_u)
check_collinearity(growth_model_u)
rsquared(growth_model_u)
summary(growth_model_u)

#view random effects levels
tab_model(growth_model_u, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, 
          show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_u))

#construct null model to perform likelihood ratio test REML must be FALSE
u_growth_nitrate_null <- lmer(formula = growth_d9 ~ temp + 
                              (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                              data = combined_growth_u, REML = FALSE)
u_growth_model2 <- lmer(formula = growth_d9 ~ nitrate + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_nitrate_null, u_growth_model2)


u_growth_temp_null <- lmer(formula = growth_d9 ~ nitrate + 
                             (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                           data = combined_growth_u, REML = FALSE)
u_growth_model4 <- lmer(formula = growth_d9 ~ nitrate + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_temp_null, u_growth_model4)


#Run lmm model for Hypnea_________________________________________________________
growth_model_h <- lmer(formula = growth_d9 ~ nitrate + temp +
                         (1 | plant_id) + (1 | run_combo) + (1 | illumination), 
                       data = combined_growth_h, REML = TRUE)

hist(resid(growth_model_h))
plot(resid(growth_model_h) ~ fitted(growth_model_h))
qqnorm(resid(growth_model_h))
qqline(resid(growth_model_h))

#check the performance of the model for dataset: u
performance::check_model(growth_model_h)
rsquared(growth_model_h)
summary(growth_model_h)

#view random effects levels
tab_model(growth_model_h, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, 
          show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_h))

#construct null model to perform likelihood ratio test REML must be FALSE
u_growth_nitrate_null <- lmer(formula = growth_d9 ~ salinity + temp + 
                                (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                              data = combined_growth_u, REML = FALSE)
u_growth_model2 <- lmer(formula = growth_d9 ~ nitrate + salinity + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_nitrate_null, u_growth_model2)

u_growth_salinity_null <- lmer(formula = growth_d9 ~ nitrate + temp + 
                                 (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                               data = combined_growth_u, REML = FALSE)
u_growth_model3 <- lmer(formula = growth_d9 ~ nitrate + salinity + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_salinity_null, u_growth_model3)

u_growth_temp_null <- lmer(formula = growth_d9 ~ nitrate + salinity + 
                             (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                           data = combined_growth_u, REML = FALSE)
u_growth_model4 <- lmer(formula = growth_d9 ~ nitrate + salinity + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_temp_null, u_growth_model4)

#MAKE PLOTS_______________________________________________________________________

library(paletteer)
# All data included in this plot
combo_growth_species_plot <- combined_growth_all %>% 
  ggplot(aes(treatment_graph, growth_d9, color = salinity, shape = temp)) + 
  geom_boxplot(linewidth=0.25) + 
  geom_point(alpha = 0.75,
             size = 2,
             position = position_jitter(width = 0.3), show.legend = TRUE) + 
  labs(x="nitrate (μmols)", y= "9-Day Growth (%)", title = "By Species") + 
  facet_wrap(~species, labeller = labeller(species = c("h" = "Hypnea", "u" = "Ulva"))) +
  scale_x_discrete(labels = c("0.5", "0.5", "14", "27", "37", "53", "53", "53", "80", "86", "80", "245", "748", "2287")) + 
  ylim(-90, 380) + stat_mean() + 
  geom_hline(yintercept=0, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  scale_color_paletteer_d("MetBrewer::Cross") +
  theme(plot.title = element_text(face = "bold", vjust = -25, hjust = 0.05), # By Species title
        legend.position = c(0.40,0.82), #move it to center of plot
        legend.box = "horizontal", #make the two legend boxes stay side-by-side
        legend.text = element_text(size = 12), #increase size for legibility
        strip.text = element_text(face = "italic", size = 12), #increase size of Hypnea and Ulva in facet headers
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
combo_growth_species_plot
ggsave("combo_all_growth_byspecies_plot.png", path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots")

# Facet to see the different salinities
combined_growth_all <- combined_growth_all %>%
  mutate(
    sal_group = case_when(
      salinity %in% c(35, 28)        ~ "High (35-28)",
      salinity %in% c(22, 21, 18)    ~ "Medium (22–18)",
      salinity %in% c(15, 11)        ~ "Low (15–11)",
      TRUE ~ NA_character_
    )
  )


growth_species_salinity_group_plot <- combined_growth_all %>% 
  ggplot(aes(salinity_graph, growth_d9, color = nitrate, shape = temp)) + 
  geom_boxplot(linewidth=0.5, size = 0.5) + 
  geom_point(alpha = 1,
             size = 2,
             position = position_jitter(width = 0.3), show.legend = TRUE) + 
  labs(x="salinity (‰)", y= "9-Day Growth (%)", title = "By Species") + 
  facet_wrap(~species, labeller = labeller(species = c("h" = "Hypnea", "u" = "Ulva"))) +
  scale_x_discrete(labels = c("35 ppt", "28 ppt", "27 ppt", "22 ppt", "21 ppt", 
                              "18 ppt", "15 ppt", "11 ppt")) + 
  ylim(-90, 380) + stat_mean() + 
  geom_hline(yintercept=0, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("#F4C40FFF", "#FE9B00FF", "#D8443CFF", "#9B3441FF", "#DE597CFF", "#AA7AA1FF", 
                       "#633372FF", "#1F6E9CFF", "#2B9B81FF", "#92C051FF")) +
  theme(plot.title = element_text(face = "bold", vjust = -25, hjust = 0.05), # By Species title
        legend.position = c(0.40,0.74), #move it to center of plot
        legend.box = "horizontal", #make the two legend boxes stay side-by-side
        legend.text = element_text(size = 12), #increase size for legibility
        strip.text = element_text(face = "italic", size = 12), #increase size of Hypnea and Ulva in facet headers
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
growth_species_salinity_group_plot

combined_growth_all %>% group_by(species, nitrate) %>% summarise_at(vars(growth_d9), list(mean = mean))
combined_growth_all %>% group_by(species, salinity) %>% summarise_at(vars(growth_d9), list(mean = mean))


#ULVA
#subset data by species
ulva <- subset(combo_growth, species == "u" & temp != 30 & treatment == )

#create subsets for the plots
ulva$treatment_graph[ulva$treatment == 0] <- "a) 0.5 μmol"
ulva$treatment_graph[ulva$treatment == 1] <- "b) 14 μmol" 
ulva$treatment_graph[ulva$treatment == 2] <- "c) 27 μmol" 
ulva$treatment_graph[ulva$treatment == 2.5] <- "d) 53 μmol" 
ulva$treatment_graph[ulva$treatment == 3] <- "e) 53 μmol"
ulva$treatment_graph[ulva$treatment == 4] <- "f) 80 μmol"
ulva$treatment_graph[ulva$treatment == 11] <- "g) 80 μmol"
ulva$treatment_graph[ulva$treatment == 13] <- "h) 245 μmol"
ulva$treatment_graph[ulva$treatment == 15] <- "i) 748 μmol"
ulva$treatment_graph[ulva$treatment == 17] <- "j) 2287 μmol"

#make a histogram of the data for ulva
hist(ulva$growth_percent, main = paste("Ulva lactuca Growth Rate (%)"), col = "olivedrab3", labels = TRUE)
#or
ulva %>% ggplot(aes(growth_percent)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction
growth_model_ulva <- lmer(formula = growth_percent ~ nitrate + temp +
                            (1 | plant_id) + (1 | run), data = ulva, REML = TRUE)

hist(resid(growth_model_ulva))
plot(resid(growth_model_ulva) ~ fitted(growth_model_ulva))
qqnorm(resid(growth_model_ulva))
qqline(resid(growth_model_ulva))

#check the performance of the model for dataset: ulva
performance::check_model(growth_model_ulva)
rsquared(growth_model_ulva)
summary(growth_model_ulva)
#view random effects levels
ranef(growth_model_ulva)
tab_model(growth_model_ulva, show.intercept = FALSE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_ulva))

#construct null model to perform likelihood ratio test REML must be FALSE
ulva_growth_nitrate_null <- lmer(formula = growth_percent ~ temp + salinity + (1 | run) + (1 | plant_id), data = ulva, REML = FALSE)
ulva_growth_model2 <- lmer(formula = growth_percent ~ nitrate + temp + salinity + (1 | run) + (1 | plant_id), data = ulva, REML = FALSE)
anova(ulva_growth_nitrate_null, ulva_growth_model2)

ulva_growth_temperature_null <- lmer(formula = growth_percent ~ nitrate + salinity +  (1 | run) + (1 | plant_id), data = ulva, REML = FALSE)
ulva_growth_model3 <- lmer(formula = growth_percent ~ nitrate + salinity + temp + (1 | run) + (1 | plant_id), data = ulva, REML = FALSE)
anova(ulva_growth_temperature_null, ulva_growth_model3)

ulva_growth_salinity_null <- lmer(formula = growth_percent ~ nitrate + temp + (1 | run) + (1 | plant_id), data = ulva, REML = FALSE)
ulva_growth_model4 <- lmer(formula = growth_percent ~ nitrate + salinity + temp + (1 | run) + (1 | plant_id), data = ulva, REML = FALSE)
anova(ulva_growth_salinity_null, ulva_growth_model4)

#plots
ulva %>% ggplot(aes(treatment_graph, growth_percent, color = salinity, shape = temp)) + 
  geom_boxplot(size=0.25) + 
  geom_point(alpha = 0.75,
             size = 3,
             position = position_jitter(width = 0.3),
             show.legend = TRUE) + 
  labs(x = "nitrate",
       y = "8-Day Growth (%)",
       title = "A",
       subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("0.5μmols", "14μmols", "27μmols", "53μmols", "53μmols",
                              "80μmols", "80μmols", "245μmols", "748μmols", "2287μmols")) + 
  ylim(-90, 380) + stat_mean() + 
  geom_hline(yintercept=0, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  scale_color_scico_d(palette = "batlow", begin = 0.2) +
  theme(legend.position = c(0.88,0.88), 
        legend.box = "horizontal",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 16, vjust = -20, hjust = 0.05),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
ggsave("ulva_growth.png", path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots")

#summarize the means
ulva %>% group_by(treatment) %>% summarise_at(vars(growth_percent), list(mean = mean))
ulva %>% group_by(temp) %>% summarise_at(vars(growth_percent), list(mean = mean))
ulva %>% group_by(salinity) %>% summarise_at(vars(growth_percent), list(mean = mean))

#HYPNEA______________________________________________________________________

hypnea <- subset(combo_growth, species == "h" & temp != 30)

#create subsets for the plots
hypnea$treatment_graph[hypnea$treatment == 0] <- "a) 0.5 μmol"
hypnea$treatment_graph[hypnea$treatment == 1] <- "b) 14 μmol" 
hypnea$treatment_graph[hypnea$treatment == 2] <- "c) 27 μmol" 
hypnea$treatment_graph[hypnea$treatment == 2.5] <- "d) 53 μmol" 
hypnea$treatment_graph[hypnea$treatment == 3] <- "e) 53 μmol"
hypnea$treatment_graph[hypnea$treatment == 4] <- "f) 80 μmol"
hypnea$treatment_graph[hypnea$treatment == 11] <- "g) 80 μmol"
hypnea$treatment_graph[hypnea$treatment == 13] <- "h) 245 μmol"
hypnea$treatment_graph[hypnea$treatment == 15] <- "i) 748 μmol"
hypnea$treatment_graph[hypnea$treatment == 17] <- "j) 2287 μmol"

hypnea %>% ggplot(aes(treatment_graph, growth_percent, color = salinity, shape = temp)) + 
  geom_boxplot(size=0.25) + 
  geom_point(alpha = 0.75, 
             size = 3, 
             position = position_jitter(width = 0.3),
             show.legend = TRUE) + 
  labs(x="nitrate",
       y= "8-Day Growth Rate (%)",
       title= "B", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("0.5μmols", "14μmols", "27μmols", "53μmols", "53μmols", "80μmols",
                              "80μmols", "245μmols", "748μmols", "2287μmols")) + 
  ylim(-90, 380) + stat_mean() + 
  #scale_color_manual(values = c("#9C0627", "#BB589F", "#F4B4E2")) +
  geom_hline(yintercept=0, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  scale_color_scico_d(palette = "batlow") +
  theme(legend.position = c(0.88,0.88), 
        legend.box = "horizontal",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 16, vjust = -20, hjust = 0.05),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
ggsave("hypnea_growth.png", path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots")

hypnea %>% group_by(treatment) %>% summarise_at(vars(growth_percent), list(mean = mean))
hypnea %>% group_by(temp) %>% summarise_at(vars(growth_percent), list(mean = mean))
hypnea %>% group_by(salinity) %>% summarise_at(vars(growth_percent), list(mean = mean))

#BELOW HAS NOT YET BEEN MODIFIED
#make a histogram of the data for hypnea
hist(hypnea$growth_rate_percent, main = paste("Hypnea musciformis Growth Rate (%)"), col = "maroon", labels = TRUE)

hypnea %>% ggplot(aes(growth_rate_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction
growth_model_hypnea <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant.ID) +
                              (1 | run) + (1 | RLC.order), data = hypnea, REML = TRUE)


#plot residuals for hypnea growth model
plot(resid(growth_model_hypnea) ~ fitted(growth_model_hypnea))
qqnorm(resid(growth_model_hypnea))
qqline(resid(growth_model_hypnea))

#check the performance of the model
performance::check_model(growth_model_hypnea)
r.squaredGLMM(growth_model_hypnea)
summary(growth_model_hypnea)
ranef(growth_model_hypnea)
tab_model(growth_model_hypnea, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_hypnea))

#construct null model to perform likelihood ratio test REML must be FALSE
hypnea_growth_treatment_null <- lmer(formula = growth_rate_percent ~ temperature + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
hypnea_growth_model2 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
anova(hypnea_growth_treatment_null, hypnea_growth_model2)
hypnea_growth_temperature_null <- lmer(formula = growth_rate_percent ~ treatment + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
hypnea_growth_model3 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
anova(hypnea_growth_temperature_null, hypnea_growth_model3)

hypnea %>% ggplot(aes(treatment_graph, growth_percent, color = salinity, shape = temp)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 1, size = 3, position = "jitter", show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "8-Day Growth Rate (%)", title= "B", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("0.5μmolN", "14μmolN", "27μmolN", "53μmolN", "53μmolN", "80μmolN", "80μmolN", "245μmolN", "748μmolN", "2287μmolN")) + 
  ylim(-90, 380) + stat_mean() + 
  #scale_color_manual(values = c("#9C0627", "#BB589F", "#F4B4E2")) +
  geom_hline(yintercept=0, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  scale_color_scico_d(palette = "batlow") +
  theme(legend.position = c(0.88,0.88),
        legend.box = "horizontal", plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
ggsave("hypnea_growth.png", path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots")


#summarize the means
hypnea %>% group_by(treatment) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))
hypnea %>% group_by(temperature) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))

