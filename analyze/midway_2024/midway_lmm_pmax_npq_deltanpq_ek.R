#Script to run model for photosynthesis data for Midway data on Chondria tumulosa
#By Angela Richards Donà
#Date: 9/13/24
#Modified with revamped functions to run all basic variables 11/25/24

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
library(hms)
library(lubridate)
library(ggeffects)


#load this file for normalized to quantum efficiency of photosynthesis per Silsbe and Kromkamp
mid_ps <- read.csv("../data/midway_2024/transformed/midway_ek_alpha_normalized_2024.csv")
glimpse(mid_ps)

mid_ps <- mid_ps %>% 
  rename(ID = specimenID)

mid_ps$date <- ymd(mid_ps$Date)

# make sure time is hms
mid_ps$rlc_end_time <- as_hms(mid_ps$rlc_end_time)

# assign run as a factor
mid_ps$run <- as.factor(mid_ps$run)

#assign temperature as a factor
mid_ps$salinity <- as.factor(mid_ps$salinity)

#assigns treatment as characters from integers then to factors
mid_ps$nitrate <- as.factor(as.character(mid_ps$treatment))

# assign deltaNPQ as a numeric
mid_ps$deltaNPQ <- as.numeric(mid_ps$deltaNPQ)

#assign new column for chronological ranking of individuals in each day_group
mid_ps <- mid_ps %>%
  group_by(day_group) %>%
  mutate(rank = rank(rlc_end_time)) %>%
  ungroup()

#use ranking to make smaller groups of ~15 minutes for rlc_order
mid_ps <- mid_ps %>%
  group_by(day_group) %>%
  mutate(rlc_order = floor((rank -1)/3) +1)

#combine nitrate and salinity for a treatment number
mid_ps <- mid_ps %>%
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

#subset for days during the experiments
mid_ps_d9 <- subset(mid_ps, rlc_day == 9) #do not separate by plant_part
mid_ps_d1 <- subset(mid_ps, rlc_day == 1) #do not separate by plant_part

#subset the plant_part for output. Use Day 9 for analysis
canopy <- subset(mid_ps_d9, plant_part == "canopy")
#OR use day 1
canopy <- subset(mid_ps_d1, plant_part == "canopy")
glimpse(canopy)

canopy$treatment_graph[canopy$treatment == 1] <- "1) 0.5μmol/35 ppt"
canopy$treatment_graph[canopy$treatment == 2] <- "2) 0.5μmol/28 ppt"
canopy$treatment_graph[canopy$treatment == 3] <- "3) 2μmol/35 ppt" 
canopy$treatment_graph[canopy$treatment == 4] <- "4) 2μmol/28 ppt"
canopy$treatment_graph[canopy$treatment == 5] <- "5) 4μmol/35 ppt"
canopy$treatment_graph[canopy$treatment == 6] <- "6) 4μmol/28 ppt"
canopy$treatment_graph[canopy$treatment == 7] <- "7) 8μmol/35 ppt"
canopy$treatment_graph[canopy$treatment == 8] <- "8) 8μmol/28 ppt"
glimpse(canopy)


under <- subset(mid_ps, rlc_day == 9 & plant_part == "under")
#OR
under <- subset(mid_ps, rlc_day == 1 & plant_part == "under")

under$treatment_graph[under$treatment == 1] <- "1) 0.5μmol/35 ppt"
under$treatment_graph[under$treatment == 2] <- "2) 0.5μmol/28 ppt"
under$treatment_graph[under$treatment == 3] <- "3) 2μmol/35 ppt" 
under$treatment_graph[under$treatment == 4] <- "4) 2μmol/28 ppt"
under$treatment_graph[under$treatment == 5] <- "5) 4μmol/35 ppt"
under$treatment_graph[under$treatment == 6] <- "6) 4μmol/28 ppt"
under$treatment_graph[under$treatment == 7] <- "7) 8μmol/35 ppt"
under$treatment_graph[under$treatment == 8] <- "8) 8μmol/28 ppt"


#add new column to subsets from time over 28 summary dataset
time_over_28 <- read_csv("../data/midway_2024/transformed/mins_28_plant_part.csv") #load dataset

mean_mins28 <- time_over_28 %>%
  select(plant_part, run, mean_mins28_plant_part) %>% #keep only relevant columns of data
  mutate(run = as.factor(run), mean_mins28_plant_part = as.factor(mean_mins28_plant_part))
  
canopy <- canopy %>%
  left_join(mean_mins28, by = c("run", "plant_part")) #join the datasets

canopy <- canopy %>%
  rename(mean_mins28 = mean_mins28_plant_part) #name is too long

under <- under %>%
  left_join(mean_mins28, by = c("run", "plant_part")) #join the datasets

under <- under %>%
  rename(mean_mins28 = mean_mins28_plant_part) #name is too long

#add growth rate from growth dataset
growth_rate <- read.csv("../data/midway_2024/input/midway_growth.csv")
growth_rate$treatment <- as.factor(growth_rate$treatment)

#make a new column for weight change (difference final from initial)
growth_rate$d9_growth_percent <- (growth_rate$final_weight - growth_rate$initial_weight) / growth_rate$initial_weight * 100
d9 <- growth_rate %>%
  select(ID, d9_growth_percent)#keep only relevant columns of data
glimpse(d9)  

canopy <- canopy %>%
  left_join(d9, by = c("ID")) #join the datasets   

under <- under %>%
  left_join(d9, by = c("ID"))


#Datasets are ready, now let's run models

check_model_fit <- function(model, terms) {
  hist(resid(model))
  plot(resid(model) ~ fitted(model))
  qqnorm(resid(model))
  qqline(resid(model))
  print(performance::check_model(model))
  print(r.squaredGLMM(model))
  #print(tab_model(model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, 
                  #show.df = TRUE, show.zeroinf = TRUE))
  print(model.frame(model))
  plot(ggpredict(model, terms = terms))
  
}

#get predictor means
get_data_means <- function(data, predictors_list, response){
  data %>% 
    group_by(across(all_of(c("nitrate", setdiff(predictors_list, "nitrate"))))) %>%
    summarise_at(all_of(response), list(mean = mean), .groups = "drop")
}


#Make a function to run the the models for likelihood ratio tests where REML must be FALSE
#construct null model to perform likelihood ratio test REML must be FALSE
compare_lmer_models <- function(data, response, predictors_list, random_effects) {
  # Ensure response is a string
  response <- as.character(response)
  
  # Create a list to store the models
  models <- list()
  
  #initialize anova
  anova_result <- list()
  

  #random effects used for all
  random_effects_part <- paste0("(1 | ", random_effects, ")", collapse = " + ")
  
  #make a formula for all predictors
  fixed_effects_all <- paste(predictors_list, collapse = " + ")
  formula <- as.formula(paste(response, "~", fixed_effects_all, "+", random_effects_part))
  
  model_all <- lmer(formula = formula, data = data, REML = FALSE)
  
 
  # Loop through predictors_list to create and fit models
  for (i in seq_along(predictors_list)) {
    fixed_effects <- paste(predictors_list[[i]], collapse = " + ")
    formula <- as.formula(paste(response, "~", fixed_effects, "+", random_effects_part))
    
    # Fit the lmer model and store in the list
    models[[i]] <- lmer(formula = formula, data = data, REML = FALSE)
    
    anova_result[[i]] <- anova(models[[i]], model_all)
  }

  data_means <- get_data_means(data, predictors_list, response)
  
  # Return results as a list
  return(list(
    models = models,
    model_all = model_all,
    anova_result = anova_result,
    data_means = data_means
  ))
}

predictors_list <- c("salinity", "nitrate")


#Pmax_____________________________________
#Inputs for Canopy/Pmax
canopy_pmax_results <- compare_lmer_models(
  data = canopy,
  response = "pmax",
  predictors_list,
  random_effects = c("plantID", "rlc_order")
)

# Access results
summary(canopy_pmax_results$models[[1]]) # Model 1 summary (salinity)
summary(canopy_pmax_results$models[[2]]) # Model 2 summary (nitrate)
summary(canopy_pmax_results$model_all) # Model 3 summary (salinity and nitrate)
canopy_pmax_results$anova_result[[1]]        # ANOVA did nitrate have an effect on Pmax?
canopy_pmax_results$anova_result[[2]]         # ANOVA did salinity have an effect on Pmax?
check_model_fit(canopy_pmax_results$model_all, terms = predictors_list)
canopy_pmax_results$data_means


#inputs for Understory/Pmax
under_pmax_results <- compare_lmer_models(
  data = under,
  response = "pmax",
  predictors_list,
  random_effects = c("mean_mins28", "rlc_order")
)

# Access results
summary(under_pmax_results$models[[1]]) # Model 1 summary
summary(under_pmax_results$models[[2]]) # Model 2 summary
summary(under_pmax_results$model_all) # Model 3 summary
under_pmax_results$anova_result[[1]]        # ANOVA did nitrate have an effect on Pmax?
under_pmax_results$anova_result[[2]]         # ANOVA did salinity have an effect on Pmax?
check_model_fit(under_pmax_results$model_all, terms = predictors_list)
under_pmax_results$data_means

#NPQmax______________________________
#inputs for Canopy/NPQmax
canopy_npqmax_results <- compare_lmer_models(
  data = canopy,
  response = "maxNPQ_Ypoint1",
  predictors_list,
  random_effects = c("plantID", "mean_mins28", "rlc_order")
)

# Access results
summary(canopy_npqmax_results$models[[1]]) # Model 1 summary
summary(canopy_npqmax_results$models[[2]]) # Model 2 summary
summary(canopy_npqmax_results$model_all) # Model 3 summary
canopy_npqmax_results$anova[[1]]        # ANOVA did nitrate have an effect on NPQmax?
canopy_npqmax_results$anova[[2]]        # ANOVA did salinity have an effect on NPQmax?
check_model_fit(canopy_npqmax_results$model_all, terms = predictors_list)
canopy_npqmax_results$data_means

#inputs for Understory/NPQmax
under_npqmax_results <- compare_lmer_models(
  data = under,
  response = "maxNPQ_Ypoint1",
  predictors_list,
  random_effects = c("plantID")
)

# Access results
summary(under_npqmax_results$models[[1]]) # Model 1 summary
summary(under_npqmax_results$models[[2]]) # Model 2 summary
summary(under_npqmax_results$model_all) # Model 3 summary
under_npqmax_results$anova[[1]]         # ANOVA did nitrate have an effect on NPQmax?
under_npqmax_results$anova[[2]]         # ANOVA did salinity have an effect on NPQmax?
check_model_fit(under_npqmax_results$model_all, terms = predictors_list)
under_npqmax_results$data_means

#deltaNPQ____________________________

#Inputs for Canopy/deltaNPQ
canopy_delta_npq_results <- compare_lmer_models(
  data = canopy,
  response = "deltaNPQ",
  predictors_list,
  random_effects = c("plantID")
)

# Access results
summary(canopy_delta_npq_results$models[[1]]) # Model 1 summary
summary(canopy_delta_npq_results$models[[2]]) # Model 2 summary
summary(canopy_delta_npq_results$model_all) # Model 3 summary
canopy_delta_npq_results$anova[[1]]         # ANOVA did nitrate have an effect on deltaNPQ?
canopy_delta_npq_results$anova[[2]]         # ANOVA did salinity have an effect on deltaNPQ?
check_model_fit(canopy_delta_npq_results$model_all, terms = predictors_list)
canopy_delta_npq_results$data_means

#Inputs for Understory/deltaNPQ
under_delta_npq_results <- compare_lmer_models(
  data = under,
  response = "deltaNPQ",
  predictors_list,
  random_effects = c("plantID")
)

# Access results
summary(under_delta_npq_results$models[[1]]) # Model 1 summary
summary(under_delta_npq_results$models[[2]]) # Model 2 summary
summary(under_delta_npq_results$model_all) # Model 3 summary
under_delta_npq_results$anova[[1]]         # ANOVA did nitrate have an effect on deltaNPQ?
under_delta_npq_results$anova[[2]]         # ANOVA did salinity have an effect on deltaNPQ?
check_model_fit(under_delta_npq_results$model_all, terms = predictors_list)
under_delta_npq_results$data_means

#Ek____________________________________________
#Inputs for Canopy/Ek
canopy_ek_results <- compare_lmer_models(
  data = canopy,
  response = "ek.est",
  predictors_list,
  random_effects = c("plantID", "mean_mins28", "rlc_order")
)

# Access results
summary(canopy_ek_results$models[[1]]) # Model 1 summary
summary(canopy_ek_results$models[[2]]) # Model 2 summary
summary(canopy_ek_results$model_all) # Model 3 summary
canopy_ek_results$anova[[1]]         # ANOVA did nitrate have an effect on Ek?
canopy_ek_results$anova[[2]]         # ANOVA did salinity have an effect on Ek?
check_model_fit(canopy_ek_results$model_all, terms = predictors_list)
canopy_ek_results$data_means

#Inputs for Understory/Ek
under_ek_results <- compare_lmer_models(
  data = under,
  response = "ek.est",
  predictors_list,
  random_effects = c("mean_mins28", "rlc_order")
)

# Access results
summary(under_ek_results$models[[1]]) # Model 1 summary
summary(under_ek_results$models[[2]]) # Model 2 summary
summary(under_ek_results$model_all) # Model 3 summary
under_ek_results$anova[[1]]         # ANOVA did nitrate have an effect on Ek?
under_ek_results$anova[[2]]         # ANOVA did salinity have an effect on Ek?
check_model_fit(under_ek_results$model_all, terms = predictors_list)
under_ek_results$data_means

#ALL OF THE ABOVE RANDOM EFFECTS HAVE BEEN CHECKED AND MAXIMIZED FOR FIT


#HISTOGRAMS and PLOTS_____________________________________
#function for raw data plots
raw_plots <- function(data, response, response2, label, pretty_color, aescolor, x, x2, y, y2, title, title2, subtitle, ylim1, ylim2, 
                      color1, color2, yint, vjust_t, hjust_t, vjust_s, hjust_s, labels1) {
  histo <- data %>%
    ggplot(aes(x = {{response}})) +
    geom_histogram(fill = pretty_color, color = "black", size = 0.25, alpha = 0.85) +
    labs(title = label) +
    theme_bw()

  plot <- data %>% 
    ggplot(aes(x = treatment_graph, y = {{response}})) + 
    geom_boxplot(size=0.5) + 
    geom_point(alpha = 0.75, size = 3, aes(color = {{aescolor}}), position = "jitter", show.legend = TRUE) + 
    labs(x = x, y = y, title = title, subtitle = subtitle) + 
    #scale_x_discrete(labels = c("0.5μmolN", "2μmolN", "4μmolN", "8μmolN")) + 
    ylim(ylim1, ylim2) + stat_mean() + 
    scale_color_manual(values = c(color1, color2)) +
    scale_x_discrete(labels = labels1) +
    geom_hline(yintercept=yint, color = "red", size = 0.5, alpha = 0.5) +
    theme_bw() +
    theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = vjust_t, hjust = hjust_t), 
          plot.subtitle = element_text(face = "italic", size = 14, vjust = vjust_s, hjust = hjust_s))
  
  lin_regr <- data %>%
    ggplot(aes(x = {{response}}, 
               y = {{response2}})) + 
    geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
    geom_smooth(method = "lm", col = "black") + 
    theme_bw() + 
    labs(title = title2, 
         x = x2, 
         y = y2) + 
    stat_regline_equation(label.x = 1, label.y = max(data$response2) * 1.1) + 
    stat_cor(label.x = 1, label.y = max(data$response2) * 1.05)
  
  return(list(
    histo = histo,
    plot = plot,
    lin_regr = lin_regr)) 
}

#INPUTS for histo and plots-----------------

#Pmax ------------------------

#Inputs for canopy/pmax
pmax_canopy <- raw_plots(
  data = canopy,
  response = pmax,
  #response2 = d1_growth_percent,
  label = "Canopy",
  pretty_color = "goldenrod1",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "Pmax (μmols electrons m-2 s-1)",
  y = "Day 1 Pmax (μmols electrons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "A - Canopy",
  title2 = "Chondria tumulosa Canopy --- Pmax vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 0, 
  ylim2 = 175,
  color1 = "goldenrod3",
  color2 = "darkgoldenrod1",
  yint = 100,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(pmax_canopy$histo)
plot(pmax_canopy$plot)
plot(pmax_canopy$lin_regr)
ggsave("pmax_canopy_d1.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#inputs for under/pmax
pmax_under <- raw_plots(
  data = under,
  response = pmax,
  #response2 = d9_growth_percent,
  label = "Understory",
  pretty_color = "maroon4",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "Pmax (μmols electrons m-2 s-1)",
  y = "Day 1 Pmax (μmols electrons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "B - Understory",
  title2 = "Chondria tumulosa Understory --- Pmax vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 0, 
  ylim2 = 175,
  color1 = "deeppink4",
  color2 = "deeppink3",
  yint = 100,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(pmax_under$histo)
plot(pmax_under$plot)
plot(pmax_under$lin_regr)
ggsave("pmax_under_d1.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#NPQmax----------------

#inputs for canopy/NPQmax
npqmax_canopy <- raw_plots(
  data = canopy,
  response = maxNPQ_Ypoint1,
  response2 = d9_growth_percent,
  label = "Canopy",
  pretty_color = "goldenrod2",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "NPQmax",
  y = "Day 9 NPQmax",
  y2 = "9-Day Growth (%)",
  title = "C - Canopy",
  title2 = "Chondria tumulosa Canopy --- NPQmax vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 0, 
  ylim2 = 2,
  color1 = "gold1",
  color2 = "gold3",
  yint = 1,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(npqmax_canopy$histo)
plot(npqmax_canopy$plot)
plot(npqmax_canopy$lin_regr)
ggsave("npqmax_canopy.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)



#inputs for under/NPQmax
npqmax_under <- raw_plots(
  data = under,
  response = maxNPQ_Ypoint1,
  response2 = d9_growth_percent,
  label = "Understory",
  pretty_color = "maroon1",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "NPQmax",
  y = "Day 9 NPQmax",
  y2 = "9-Day Growth (%)",
  title = "D - Understory",
  title2 = "Chondria tumulosa Understory --- NPQmax vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 0, 
  ylim2 = 2,
  color1 = "hotpink2",
  color2 = "hotpink4",
  yint = 1,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(npqmax_under$histo)
plot(npqmax_under$plot)
plot(npqmax_under$lin_regr)
ggsave("npqmax_under.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#deltaNPQ----------------

#inputs for canopy/deltaNPQ
delta_npq_canopy <- raw_plots(
  data = canopy,
  response = deltaNPQ,
  response2 = d9_growth_percent,
  label = "Canopy",
  pretty_color = "goldenrod3",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "ΔNPQ",
  y = "Day 9 ΔNPQ",
  y2 = "9-Day Growth (%)",
  title = "E - Canopy",
  title2 = "Chondria tumulosa Canopy --- ΔNPQ vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 0, 
  ylim2 = 2,
  color1 = "lightgoldenrod2",
  color2 = "lightgoldenrod3",
  yint = 0.5,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(delta_npq_canopy$histo)
plot(delta_npq_canopy$plot)
plot(delta_npq_canopy$lin_regr)
ggsave("delta_npq_canopy.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#inputs for under/deltaNPQ
delta_npq_under <- raw_plots(
  data = under,
  response = deltaNPQ,
  response2 = d9_growth_percent,
  label = "Understory",
  pretty_color = "maroon3",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "ΔNPQ",
  y = "Day 9 ΔNPQ",
  y2 = "9-Day Growth (%)",
  title = "F - Understory",
  title2 = "Chondria tumulosa Understory --- ΔNPQ vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 0, 
  ylim2 = 2,
  color1 = "hotpink1",
  color2 = "hotpink3",
  yint = 0.5,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(delta_npq_under$histo)
plot(delta_npq_under$plot)
plot(delta_npq_under$lin_regr)
ggsave("delta_npq_under.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#Ek----------------

#inputs for canopy/Ek
ek_canopy <- raw_plots(
  data = canopy,
  response = ek.est,
  response2 = d9_growth_percent,
  label = "Canopy",
  pretty_color = "goldenrod4",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "Ek (μmols electrons m-2 s-1)",
  y = "Day 9 Ek (μmols electrons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "G - Canopy",
  title2 = "Chondria tumulosa Canopy --- Ek vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 50, 
  ylim2 = 300,
  color1 = "khaki4",
  color2 = "khaki3",
  yint = 150,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(ek_canopy$histo)
plot(ek_canopy$plot)
plot(ek_canopy$lin_regr)
ggsave("ek_canopy.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#inputs for under/Ek
ek_under <- raw_plots(
  data = under,
  response = ek.est,
  response2 = d9_growth_percent,
  label = "Understory",
  pretty_color = "maroon3",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "Pmax (μmols electrons m-2 s-1)",
  y = "Day 9 Ek (μmols electrons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "H - Understory",
  title2 = "Chondria tumulosa Understory --- Ek vs 9-Day Growth (%)",
  subtitle = "Chondria tumulosa",
  ylim1 = 50, 
  ylim2 = 300,
  color1 = "palevioletred3",
  color2 = "palevioletred4",
  yint = 150,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("0.5 μmol", "0.5 μmol", 
              "2 μmol", "2 μmol",
              "4 μmol", "4 μmol",
              "8 μmol", "8 μmol")
)
plot(ek_under$histo)
plot(ek_under$plot)
plot(ek_under$lin_regr)
ggsave("ek_under.png", path = "midway_2024/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)



#Linear regression Pmax vs Growth--------------------------------------------------------------------------------

#plot a regression between the photosynthetic independent variables of interest and growth rate


under_growth_pmax_plot <- under %>%
  ggplot(aes(x=pmax, 
             y=d9_growth_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + 
  theme_bw() + 
  labs(title = "Chondria tumulosa Understory --- Pmax vs 9-Day Growth (%)", 
       x = "Pmax (μmols electrons m-2 s-1)", 
       y = "9-Day Growth (%)") + 
  stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
under_growth_pmax_plot


#--------------------NPQmax--------------------------






#summarize the means for pmax
canopy %>% group_by(nitrate, salinity) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))
under %>% group_by(nitrate, salinity) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))


#Linear regression NPQ vs Growth and Apices
#plot a regression between the photosynthetic independent variables of interest and growth rate
#canopy_growth_NPQmax_graph <- ggplot(canopy_sub, aes(x=maxNPQ_Ypoint1, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Chondria tumulosa NPQmax vs 9-Day Growth (%)", x = " NPQmax (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
#canopy_growth_NPQmax_graph



#-----------------delta NPQ------------------


#summarize the means for deltaNPQ
canopy %>% group_by(treatment) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
canopy %>% group_by(temp) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
#canopy %>% group_by(treatment, rlc_day) %>% summarise_at(vars(pmax), list(mean = mean))

#Linear regression deltaNPQ vs Growth and Apices

#plot a regression between the photosynthetic independent variables of interest and growth rate
canopy_growth_deltaNPQ_graph <- ggplot(canopy_sub, aes(x=deltaNPQ, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Chondria tumulosa Delta NPQ vs 9-Day Growth (%)", x = " Delta NPQ (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
canopy_growth_deltaNPQ_graph


