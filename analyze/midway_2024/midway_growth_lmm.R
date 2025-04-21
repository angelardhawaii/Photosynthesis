#This script brings in growth data for Midway Chondria tumulosa analysis
#By Angela Richards Donà
#Date: 9/13/24

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
mid_growth <- read.csv("../data/midway_2024/input/midway_growth.csv")

#make a new column for weight change (difference final from initial)
mid_growth$d9_growth_percent <- round(((mid_growth$final_weight - mid_growth$initial_weight) / mid_growth$initial_weight * 100), digits = 2)
mid_growth$d5_growth_percent <- round(((mid_growth$d5_weight - mid_growth$initial_weight) / mid_growth$initial_weight * 100), digits = 2)

#add a new column that gets rid of characters in ID column
mid_growth$plant_ID <- as.factor(substr(mid_growth$ID, 4, 5))

#make a new column for run using the third value in the ID, then change the 1 to 2 and 2 to 3
mid_growth$run <- substr(mid_growth$ID, 3, 3)
mid_growth$run <- as.factor(ifelse(mid_growth$run == 1, 2, 3)) #actual run numbers are 2 and 3

#assigns treatment as characters from integers then to factors
mid_growth$nitrate <- as.factor(as.character(mid_growth$treatment))

#get rid of the ppt in salinity
mid_growth$salinity <- substr(mid_growth$salinity, 1, 2)

#combine nitrate and salinity for a treatment number
mid_growth <- mid_growth %>%
  mutate(treatment = case_when(
    nitrate == 1 & salinity == 35 ~ "1",
    nitrate == 1 & salinity == 28 ~ "2",
    nitrate == 2 & salinity == 35 ~ "3",
    nitrate == 2 & salinity == 28 ~ "4",
    nitrate == 3 & salinity == 35 ~ "5",
    nitrate == 3 & salinity == 28 ~ "6",
    nitrate == 4 & salinity == 35 ~ "7",
    nitrate == 4 & salinity == 28 ~ "8",
  ))

#create subsets for the plots
#toggle between the plant_part for output. Use Day 9 for final analysis
canopy_g <- subset(mid_growth, plant_part == "canopy" & d9_growth_percent!= -42.45)
canopy_g$treatment_graph[canopy_g$treatment == 1] <- "1) 0.5umol/35 ppt"
canopy_g$treatment_graph[canopy_g$treatment == 2] <- "2) 0.5umol/28 ppt"
canopy_g$treatment_graph[canopy_g$treatment == 3] <- "3) 2umol/35 ppt" 
canopy_g$treatment_graph[canopy_g$treatment == 4] <- "4) 2umol/28 ppt"
canopy_g$treatment_graph[canopy_g$treatment == 5] <- "5) 4umol/35 ppt"
canopy_g$treatment_graph[canopy_g$treatment == 6] <- "6) 4umol/28 ppt"
canopy_g$treatment_graph[canopy_g$treatment == 7] <- "7) 8umol/35 ppt"
canopy_g$treatment_graph[canopy_g$treatment == 8] <- "8) 8umol/28 ppt"


under_g <- subset(mid_growth, plant_part == "under" & d9_growth_percent != -62.82& 
                    d9_growth_percent != -32.50)
under_g$treatment_graph[under_g$treatment == 1] <- "1) 0.5umol/35 ppt"
under_g$treatment_graph[under_g$treatment == 2] <- "2) 0.5umol/28 ppt"
under_g$treatment_graph[under_g$treatment == 3] <- "3) 2umol/35 ppt" 
under_g$treatment_graph[under_g$treatment == 4] <- "4) 2umol/28 ppt"
under_g$treatment_graph[under_g$treatment == 5] <- "5) 4umol/35 ppt"
under_g$treatment_graph[under_g$treatment == 6] <- "6) 4umol/28 ppt"
under_g$treatment_graph[under_g$treatment == 7] <- "7) 8umol/35 ppt"
under_g$treatment_graph[under_g$treatment == 8] <- "8) 8umol/28 ppt"

#add new column to subsets from time over 28 summary dataset
time_over_28 <- read_csv("../data/midway_2024/transformed/mins_28_plant_part.csv") #load dataset

mean_mins28 <- time_over_28 %>%
  select(plant_part, run, mean_mins28_plant_part) %>% #keep only relevant columns of data
  mutate(run = as.factor(run), mean_mins28_plant_part = as.factor(mean_mins28_plant_part))

canopy_g <- canopy_g %>%
  left_join(mean_mins28, by = c("run", "plant_part")) #join the datasets

canopy_g <- canopy_g %>%
  rename(mean_mins28 = mean_mins28_plant_part) #name is too long
glimpse(canopy_g)

under_g <- under_g %>%
  left_join(mean_mins28, by = c("run", "plant_part")) #join the datasets

under_g <- under_g %>%
  rename(mean_mins28 = mean_mins28_plant_part) #name is too long
glimpse(under_g)

#Is there a sig difference between canopy and understory growth?
#CANOPY

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
canopy_growth_plot <- canopy_g %>% 
  ggplot(aes(treatment_graph, d9_growth_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = salinity), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "A", subtitle = "Chondria tumulosa -- Canopy") + 
  scale_x_discrete(labels = c("35 ppt/0.5 μmol N", "28 ppt/0.5 μmol N", "35 ppt/2 μmol N", "28 ppt/2 μmol N", 
                              "35 ppt/4 μmol N", "28 ppt/4 μmol N", "35 ppt/8 μmol N", "28 ppt/8 μmol N")) + 
  ylim(-45, 50) + stat_mean() + 
  scale_color_manual(values = c("goldenrod1", "darkgoldenrod3")) +
  geom_hline(yintercept=0, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position.inside = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
canopy_growth_plot
ggsave("canopy_growth_plot.png", path = "midway_2024/plots/")

#summarize the means for canopy
mid_growth %>% group_by(nitrate, plant_part) %>% summarise_at(vars(d9_growth_percent), list(mean = mean))
mid_growth %>% group_by(salinity, plant_part) %>% summarise_at(vars(d9_growth_percent), list(mean = mean))



#UNDERSTORY

#-----------make a histogram of the growth rate data
under_g %>% ggplot(aes(d5_growth_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()
under_g %>% ggplot(aes(d9_growth_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", linewidth = 0.25, alpha = 0.85) +
  theme_bw()

#run model without RLC_order as this has little effect
#mean_min28 used as random effect in lieu of run
growth_model_under <- lmer(formula = d9_growth_percent ~ salinity + nitrate +
                              (1 | plant_ID) + (1 | run), data = under_g, REML = TRUE)

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
                              (1 | plant_ID) + (1 | run), data = under_g, REML = FALSE)
under_model2 <- lmer(formula = d9_growth_percent ~ nitrate + salinity +
                        (1 | plant_ID) + (1 | run), data = under_g, REML = FALSE)
anova(under_nitrate_null, under_model2)
under_salinity_null <- lmer(formula = d9_growth_percent ~ nitrate + (1 | plant_ID) + (1 | run), data = under_g, REML = FALSE)
under_model3 <- lmer(formula = d9_growth_percent ~ nitrate + salinity + (1 | plant_ID) + (1 | run), data = under_g, REML = FALSE)
anova(under_salinity_null, under_model3)

#plot under growth
under_growth_plot <- under_g %>% 
  ggplot(aes(treatment_graph, d9_growth_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = salinity), position = "jitter", show.legend = TRUE) + 
  labs(x="treatment", y= "9-Day Growth (%)", title= "B", subtitle = "Chondria tumulosa -- Understory") + 
  scale_x_discrete(labels = c("35 ppt/0.5 μmol N", "28 ppt/0.5 μmol N", "35 ppt/2 μmol N", "28 ppt/2 μmol N", 
                              "35 ppt/4 μmol N", "28 ppt/4 μmol N", "35 ppt/8 μmol N", "28 ppt/8 μmol N")) + 
  ylim(-45, 50) + stat_mean() + 
  scale_color_manual(values = c("maroon", "maroon2")) +
  geom_hline(yintercept=0, color = "red", linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
under_growth_plot
ggsave("under_growth_plot.png", path = "midway_2024/plots/")

#summarize the means for under
mid_growth %>% group_by(nitrate, plant_part) %>% summarise_at(vars(d9_growth_percent), list(mean = mean))
mid_growth %>% group_by(salinity, plant_part) %>% summarise_at(vars(d9_growth_percent), list(mean = mean))







