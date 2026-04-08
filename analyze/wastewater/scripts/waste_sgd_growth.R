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
# Filter data before proceeding
combined_growth_all <- combined_growth_all %>%
  filter(
    # keep:
    (
      year %in% c(2021, 2022) & 
        nitrate == "80" & 
        temp != 30
    ) |
      year == 2025
  )  
#make a list to be used for plots that orders chronologically first then by nitrate value low to high
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "f"] <- "i) 80 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "k"] <- "k) 80 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "l"] <- "l) 245 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "m"] <- "m) 748 μmol"
combined_growth_all$treatment_graph[combined_growth_all$treat_letter == "n"] <- "n) 2287 μmol"

combined_growth_all$salinity_graph[combined_growth_all$salinity == "27"] <- "c) 27 ppt"
combined_growth_all$salinity_graph[combined_growth_all$salinity == "22"] <- "d) 22 ppt" 
combined_growth_all$salinity_graph[combined_growth_all$salinity == "18"] <- "f) 18 ppt"
combined_growth_all$salinity_graph[combined_growth_all$salinity == "11"] <- "h) 11 ppt"

#Subset by species
combined_growth_u <- combined_growth_all %>%
  filter(species == "u")

combined_growth_h <- combined_growth_all %>%
  filter(species == "h")

xtabs(~ salinity + temp + year, data = combined_growth_u)
xtabs(~ salinity + temp + year, data = combined_growth_h)

#Run lmm model for Ulva_________________________________________________________
growth_model_u <- lmer(formula = growth_d9 ~ nitrate + salinity + temp + (1 | run_combo),
                       data = combined_growth_u, REML = TRUE)

hist(resid(growth_model_u))
plot(resid(growth_model_u) ~ fitted(growth_model_u))
qqnorm(resid(growth_model_u))
qqline(resid(growth_model_u))

#check the performance of the model for dataset: u
performance::check_model(growth_model_u)
check_collinearity(growth_model_u)
r.squaredGLMM(growth_model_u)
summary(growth_model_u)

#view random effects levels
tab_model(growth_model_u, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, 
          show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_u))

#construct null model to perform likelihood ratio test REML must be FALSE
# Does nitrate affect Ulva growth?
u_growth_nitrate_null <- lmer(formula = growth_d9 ~ temp + salinity +
                              (1 | run_combo), 
                              data = combined_growth_u, REML = FALSE)
u_growth_model2 <- lmer(formula = growth_d9 ~ nitrate + temp + salinity +
                          (1 | run_combo), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_nitrate_null, u_growth_model2)

# Does temperature affect Ulva growth?
u_growth_temp_null <- lmer(formula = growth_d9 ~ nitrate + salinity +
                             (1 | run_combo), 
                           data = combined_growth_u, REML = FALSE)
u_growth_model4 <- lmer(formula = growth_d9 ~ nitrate + temp + salinity +
                          (1 | run_combo), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_temp_null, u_growth_model4)

# Does salinity affect Ulva growth?
u_growth_salinity_null <- lmer(formula = growth_d9 ~ nitrate + temp +
                             (1 | run_combo), 
                           data = combined_growth_u, REML = FALSE)
u_growth_model4 <- lmer(formula = growth_d9 ~ nitrate + temp + salinity +
                          (1 | run_combo), 
                        data = combined_growth_u, REML = FALSE)
anova(u_growth_salinity_null, u_growth_model4)


#Run lmm model for Hypnea_________________________________________________________
growth_model_h <- lmer(formula = growth_d9 ~ nitrate + temp + salinity +
                         (1 | plant_id) + (1 | run_combo) + (1 | illumination), 
                       data = combined_growth_h, REML = TRUE)

hist(resid(growth_model_h))
plot(resid(growth_model_h) ~ fitted(growth_model_h))
qqnorm(resid(growth_model_h))
qqline(resid(growth_model_h))

#check the performance of the model for dataset: u
performance::check_model(growth_model_h)
r.squaredGLMM(growth_model_h)
summary(growth_model_h)

#view random effects levels
tab_model(growth_model_h, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, 
          show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_h))

#construct null model to perform likelihood ratio test REML must be FALSE
# Does nitrate affect Hypnea growth?
h_growth_nitrate_null <- lmer(formula = growth_d9 ~ salinity + temp + 
                                (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                              data = combined_growth_h, REML = FALSE)
h_growth_model2 <- lmer(formula = growth_d9 ~ nitrate + salinity + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_h, REML = FALSE)
anova(h_growth_nitrate_null, h_growth_model2)

#Does salinity affect Hypnea growth?
h_growth_salinity_null <- lmer(formula = growth_d9 ~ nitrate + temp + 
                                 (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                               data = combined_growth_h, REML = FALSE)
h_growth_model3 <- lmer(formula = growth_d9 ~ nitrate + salinity + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_h, REML = FALSE)
anova(h_growth_salinity_null, h_growth_model3)

h_growth_temp_null <- lmer(formula = growth_d9 ~ nitrate + salinity + 
                             (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                           data = combined_growth_h, REML = FALSE)
h_growth_model4 <- lmer(formula = growth_d9 ~ nitrate + salinity + temp + 
                          (1 | run_combo) + (1 | plant_id) + (1 | illumination), 
                        data = combined_growth_h, REML = FALSE)
anova(h_growth_temp_null, h_growth_model4)

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
  scale_x_discrete(labels = c("80", "80", "245", "748", "2287")) + 
  ylim(-50, 380) + stat_mean() + 
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

#_________________________________________________________________________________
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


