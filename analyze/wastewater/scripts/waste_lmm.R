#Script to run model for photosynthesis data for Hypnea and Ulva wastewater data
#By Angela Richards Donà
#Date: 1/16/25
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
waste_p1_ps <- read.csv("../data/wastewater/transformed/waste_p1_ek_alpha_normalized_2025.csv")
glimpse(waste_p1_ps)

waste_p1_ps$Date <- ymd(waste_p1_ps$Date)

# make sure time is hms
waste_p1_ps$rlc_end_time <- as_hms(waste_p1_ps$rlc_end_time)

# assign run as a factor
#waste_p1_ps$run <- as.factor(waste_p1_ps$run)

#assign temperature as a factor
#waste_p1_ps$salinity <- as.factor(waste_p1_ps$salinity)

#assigns treatment as characters from integers then to factors
waste_p1_ps$treatment <- as.factor(as.character(waste_p1_ps$treatment))

# assign deltaNPQ as a numeric
waste_p1_ps$deltaNPQ <- as.numeric(waste_p1_ps$deltaNPQ)

#assign new column for chronological ranking of individuals in each day_group
waste_p1_ps <- waste_p1_ps %>%
  group_by(Date) %>%
  mutate(rank = rank(rlc_end_time))


#use ranking to make smaller groups of ~15 minutes for rlc_order
waste_p1_ps <- waste_p1_ps %>%
  group_by(Date) %>%
  mutate(rlc_order = floor((rank -1)/2) +1)

#waste_p1_ps_d9 <- subset(waste_p1_ps, rlc_day == 9) #only d9 is included in dataset

#toggle between the species for output.
ulva_ps <- subset(waste_p1_ps, species == "u")
ulva_ps$treatment_graph[ulva_ps$treatment == 1] <- "1) 80 μmol"
ulva_ps$treatment_graph[ulva_ps$treatment == 2] <- "2) 140 μmol"
ulva_ps$treatment_graph[ulva_ps$treatment == 3] <- "3) 245 μmol" 
ulva_ps$treatment_graph[ulva_ps$treatment == 4] <- "4) 428 μmol"
ulva_ps$treatment_graph[ulva_ps$treatment == 5] <- "5) 748 μmol"
ulva_ps$treatment_graph[ulva_ps$treatment == 6] <- "6) 1308 μmol"
ulva_ps$treatment_graph[ulva_ps$treatment == 7] <- "7) 2287 μmol"
ulva_ps$treatment_graph[ulva_ps$treatment == 8] <- "8) 4000 μmol"
glimpse(ulva_ps)

hypnea_ps <- subset(waste_p1_ps, species == "h")
hypnea_ps$treatment_graph[hypnea_ps$treatment == 1] <- "1) 80 μmol"
hypnea_ps$treatment_graph[hypnea_ps$treatment == 2] <- "2) 140 μmol"
hypnea_ps$treatment_graph[hypnea_ps$treatment == 3] <- "3) 245 μmol" 
hypnea_ps$treatment_graph[hypnea_ps$treatment == 4] <- "4) 428 μmol"
hypnea_ps$treatment_graph[hypnea_ps$treatment == 5] <- "5) 748 μmol"
hypnea_ps$treatment_graph[hypnea_ps$treatment == 6] <- "6) 1308 μmol"
hypnea_ps$treatment_graph[hypnea_ps$treatment == 7] <- "7) 2287 μmol"
hypnea_ps$treatment_graph[hypnea_ps$treatment == 8] <- "8) 4000 μmol"

#SKIP models for phase I
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
hypnea_ps_pmax_results <- compare_lmer_models(
  data = hypnea_ps,
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

#ALL OF THE ABOVE RANDOM EFFECTS MuST BE CHECKED AND MAXIMIZED FOR FIT


#HISTOGRAMS and PLOTS_____________________________________
#function for raw data plots
raw_plots <- function(data, response, response2, label, pretty_color, aescolor, x, x2, y, y2, title, title2, subtitle, ylim1, ylim2, 
                      color1, color2, yint, vjust_t, hjust_t, vjust_s, hjust_s, labels1) {
  histo <- data %>%
    ggplot(aes(x = {{response}})) +
    geom_histogram(fill = pretty_color, color = "black", linewidth = 0.25, alpha = 0.85) +
    labs(title = label) +
    theme_bw()
  
  plot <- data %>% 
    ggplot(aes(x = treatment_graph, y = {{response}})) + 
    geom_boxplot(size=0.5) + 
    geom_point(alpha = 0.75, size = 3, aes(color = {{aescolor}}), position = "jitter", show.legend = TRUE) + 
    labs(x = x, y = y, title = title, subtitle = subtitle) + 
    #scale_x_discrete(labels = c("0.5μmolN", "2μmolN", "4μmolN", "8μmolN")) + 
    ylim(ylim1, ylim2) + stat_mean() + 
    #scale_color_manual(values = c(color1, color2)) +
    scale_x_discrete(labels = labels1) +
    geom_hline(yintercept=yint, color = "red", linewidth = 0.5, alpha = 0.5) +
    theme_bw(axis.text.x = element_text(size = 14)) +
    theme(legend.position.inside = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = vjust_t, hjust = hjust_t), 
          plot.subtitle = element_text(face = "italic", size = 14, vjust = vjust_s, hjust = hjust_s))
  
 # lin_regr <- data %>%
  #  ggplot(aes(x = {{response}}, 
   #            y = {{response2}})) + 
  #  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  #  geom_smooth(method = "lm", col = "black") + 
  #  theme_bw() + 
  #  labs(title = title2, 
   #      x = x2, 
    #     y = y2) + 
  #  stat_regline_equation(label.x = 1, label.y = max(data$response2) * 1.1) + 
  #  stat_cor(label.x = 1, label.y = max(data$response2) * 1.05)
  
  return(list(
    #histo = histo,
    plot = plot))
    #lin_regr = lin_regr)) 
}

#INPUTS for histo and plots-----------------

#Pmax ------------------------

#Inputs for Ulva/pmax
pmax_ulva <- raw_plots(
  data = ulva_ps,
  response = pmax,
  response2 = growth,
  label = "Ulva lactuca",
  pretty_color = "goldenrod1",
  aescolor = treatment,
  x = "Nitrate",
  x2 = "Pmax (μmols electrons m-2 s-1)",
  y = "Day 9 Pmax (μmols electrons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "A - Ulva lactuca",
  title2 = "Ulva lactuca --- Pmax vs 9-Day Growth (%)",
  subtitle = "Pmax",
  ylim1 = 0, 
  ylim2 = 30,
  #color1 = "goldenrod3",
  #color2 = "darkgoldenrod1",
  yint = 20,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("80 μmol", "140 μmol", 
              "245 μmol", "428 μmol",
              "748 μmol", "1308 μmol",
              "2287 μmol", "4000 μmol")
)
plot(pmax_ulva$histo)
plot(pmax_ulva$plot)
plot(pmax_ulva$lin_regr)
ggsave("pmax_ulva_n6.png", path = "wastewater/plots/")#, 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#inputs for hypnea/pmax
pmax_hypnea <- raw_plots(
  data = hypnea_ps,
  response = pmax,
  response2 = growth,
  label = "Hypnea musciformis",
  pretty_color = "maroon4",
  aescolor = treatment,
  x = "Nitrate",
  x2 = "Pmax (μmols electrons m-2 s-1)",
  y = "Day 9 Pmax (μmols electrons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "B - Hypnea musciformis",
  title2 = "Hypnea musciformis --- Pmax vs 9-Day Growth (%)",
  subtitle = "Pmax",
  ylim1 = 0, 
  ylim2 = 30,
  color1 = "deeppink4",
  color2 = "deeppink3",
  yint = 20,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("80 μmol", "140 μmol", 
              "245 μmol", "428 μmol",
              "748 μmol", "1308 μmol",
              "2287 μmol", "4000 μmol")
)
plot(pmax_under$histo)
plot(pmax_hypnea$plot)
plot(pmax_under$lin_regr)
ggsave("pmax_hypnea_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#NPQmax----------------

#inputs for canopy/NPQmax
npqmax_ulva <- raw_plots(
  data = ulva_ps,
  response = maxNPQ_Ypoint1,
  response2 = growth,
  label = "Ulva lactuca",
  pretty_color = "goldenrod2",
  aescolor = treatment,
  x = "Nitrate",
  x2 = "NPQmax",
  y = "Day 9 NPQmax",
  y2 = "9-Day Growth (%)",
  title = "C - Ulva lactuca",
  title2 = "Ulva lactuca --- NPQmax vs 9-Day Growth (%)",
  subtitle = "NPQmax",
  ylim1 = 0, 
  ylim2 = 1.4,
  color1 = "gold1",
  color2 = "gold3",
  yint = 0.3,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("80 μmol", "140 μmol", 
              "245 μmol", "428 μmol",
              "748 μmol", "1308 μmol",
              "2287 μmol", "4000 μmol")
)
plot(npqmax_canopy$histo)
plot(npqmax_ulva$plot)
plot(npqmax_canopy$lin_regr)
ggsave("npqmax_ulva_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)



#inputs for under/NPQmax
npqmax_hypnea <- raw_plots(
  data = hypnea_ps,
  response = maxNPQ_Ypoint1,
  response2 = growth,
  label = "Hypnea musciformis",
  pretty_color = "maroon1",
  aescolor = treatment,
  x = "Nitrate",
  x2 = "NPQmax",
  y = "Day 9 NPQmax",
  y2 = "9-Day Growth (%)",
  title = "D - Hypnea musciformis",
  title2 = "Hypnea musciformis --- NPQmax vs 9-Day Growth (%)",
  subtitle = "NPQmax",
  ylim1 = 0, 
  ylim2 = 1.4,
  color1 = "hotpink2",
  color2 = "hotpink4",
  yint = 0.84,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("80 μmol", "140 μmol", 
              "245 μmol", "428 μmol",
              "748 μmol", "1308 μmol",
              "2287 μmol", "4000 μmol")
)
plot(npqmax_under$histo)
plot(npqmax_hypnea$plot)
plot(npqmax_under$lin_regr)
ggsave("npqmax_hypnea_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#deltaNPQ----------------

#inputs for canopy/deltaNPQ
delta_npq_ulva <- raw_plots(
  data = ulva_ps,
  response = deltaNPQ,
  response2 = growth,
  label = "Ulva lactuca",
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
ggsave("delta_npq_canopy.png", path = "wastewater/plots/", 
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
ggsave("delta_npq_under.png", path = "wastewater/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#Ek----------------

#inputs for canopy/Ek
ek_ulva <- raw_plots(
  data = ulva_ps,
  response = ek.est,
  response2 = growth,
  label = "Ulva lactuca",
  pretty_color = "goldenrod4",
  aescolor = treatment,
  x = "Nitrate",
  x2 = "Ek (μmols photons m-2 s-1)",
  y = "Day 9 Ek (μmols photons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "G - Ulva lactuca",
  title2 = "Ulva lactuca --- Ek vs 9-Day Growth (%)",
  subtitle = "Ek",
  ylim1 = 0, 
  ylim2 = 50,
  color1 = "khaki4",
  color2 = "khaki3",
  yint = 25,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("80 μmol", "140 μmol", 
              "245 μmol", "428 μmol",
              "748 μmol", "1308 μmol",
              "2287 μmol", "4000 μmol")
)
plot(ek_canopy$histo)
plot(ek_ulva$plot)
plot(ek_canopy$lin_regr)
ggsave("ek_ulva_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#inputs for under/Ek
ek_hypnea <- raw_plots(
  data = hypnea_ps,
  response = ek.est,
  response2 = growth,
  label = "Hypnea musciformis",
  pretty_color = "maroon3",
  aescolor = treatment,
  x = "Nitrate",
  x2 = "Ek (μmols photons m-2 s-1)",
  y = "Day 9 Ek (μmols photons m-2 s-1)",
  y2 = "9-Day Growth (%)",
  title = "H - Hypnea musciformis",
  title2 = "Hypnea musciformis --- Ek vs 9-Day Growth (%)",
  subtitle = "Ek",
  ylim1 = 0, 
  ylim2 = 50,
  color1 = "palevioletred3",
  color2 = "palevioletred4",
  yint = 25,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -20,
  hjust_s = 0.05,
  labels1 = c("80 μmol", "140 μmol", 
              "245 μmol", "428 μmol",
              "748 μmol", "1308 μmol",
              "2287 μmol", "4000 μmol")
)
plot(ek_under$histo)
plot(ek_hypnea$plot)
plot(ek_under$lin_regr)
ggsave("ek_hypnea_n6.png", path = "wastewater/plots/") 
      # width = 7, height = 6, units = "in", dpi = 300, scale = 1)



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






#summarize the means for NPQmax
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


