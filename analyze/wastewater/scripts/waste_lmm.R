#Script to run model for photosynthesis data for Hypnea and Ulva wastewater data
# This script pulls in the combined dataset for both the 2023 SGD and wastewater data
#By Angela Richards Donà
#Date: 1/16/25
#Modified with revamped functions to run all basic variables 11/25/24
#Modified for combined data 08/26/25

#load libraries for analysis and wrangling
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
library(gt)
library(purrr)
library(stringr)
library(tidyr)
library(hms)
library(lubridate)
library(ggeffects)


#load this file for normalized to quantum efficiency of photosynthesis per Silsbe and Kromkamp
ps_combo <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/2021_2025_combo_ek_alpha_normalized.csv")


ps_combo <- ps_combo %>%
  mutate(date = ymd(date))

#unique_id looking funky, redo
ps_combo <- ps_combo %>%
  mutate(unique_id = paste(date, tolower(id), sep = "-"))

# make sure time is hms
ps_combo$rlc_end_time <- as_hms(ps_combo$rlc_end_time)
ps_combo$rlc_time <- as_hms(ps_combo$rlc_time)

# assign list of variables as factor
factor_vars <- c("treat_letter", "treatment", "salinity", "temp", "rlc_day", "nitrate")

ps_combo <- ps_combo %>%
  mutate(across(all_of(factor_vars), as.factor))

#assign new column for chronological ranking of individuals in each day_group
#use ranking to make smaller groups of ~15 minutes for rlc_order
ps_combo <- ps_combo %>%
  group_by(date) %>%
  arrange(rlc_end_time, .by_group = TRUE) %>%
  mutate(
    rank = dplyr::row_number(),                 # stable rank by time
    rlc_order1 = floor((rank - 1) / 3) + 1      # 3 steps ≈ 15 minutes
  ) %>%
  ungroup()

#Treatment Graph is set up to have predictors in order of my choosing
ps_combo <- ps_combo %>%
  mutate(treatment_graph = case_when(
    treat_letter == "a" ~ "1) 0.5 μmol",
    treat_letter == "b" ~ "2) 14 μmol",
    treat_letter == "c" ~ "3) 27 μmol",
    treat_letter == "d" ~ "4) 53 μmol",
    treat_letter == "e" ~ "5) 53 μmol",
    treat_letter == "f" ~ "6) 80 μmol",
    treat_letter == "k" ~ "7) 80 μmol",
    treat_letter == "l" ~ "8) 245 μmol",
    treat_letter == "m" ~ "9) 748 μmol",
    treat_letter == "n" ~ "10) 2287 μmol"
  ))


#subset the nitrate levels needed removing those used for 2023 paper irrelevant to wastewater
# and keep only day 9 for analysis
# now that i have spent an eternity combining the two datasets, I will filter out the 2021 
#data because it is giving major problems of singularity and rank deficiency
ps_combo <- ps_combo %>%
  mutate(year = lubridate::year(date))
ps_combo_d9 <- ps_combo %>%
  filter(!year %in% c(2021, 2022)) %>%
  filter(rlc_day == 9) #only d9 is included in dataset
  
# subset the species
ulva_ps <- subset(ps_combo_d9, species == "u")

hypnea_ps <- subset(ps_combo_d9, species == "h")



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
  group_vars <- c("nitrate", setdiff(predictors_list, "nitrate"))
  data %>% 
    group_by(across(all_of(group_vars))) %>%
    summarise(!!paste0(response, "_mean") := mean(.data[[response]], na.rm = TRUE),
              .groups = "drop")
}

#Make a function to run the the models for likelihood ratio tests where REML must be FALSE
#construct null model to perform likelihood ratio test REML must be FALSE

compare_lmer_models <- function(data, response, predictors_list, random_effects) {
    # Keep only rows usable by all models
    needed <- c(response, predictors_list, random_effects)
    data_clean <- tidyr::drop_na(data, all_of(unique(c(response, predictors_list))))
    data_clean <- droplevels(data_clean)  # drop unused factor levels 
  # Build random-effects string
    re_str <- paste0("(1 | ", random_effects, ")", collapse = " + ")
  
  # Full model
    fixed_all <- paste(predictors_list, collapse = " + ")
    f_full <- as.formula(paste(response, "~", fixed_all, "+", re_str))
    model_all <- lme4::lmer(f_full, data = data_clean, REML = FALSE)
  
  # Drop-one models + LRTs
    models <- list()
    anova_result <- list()
    
    for (dropped in predictors_list) {
      fixed_reduced <- paste(setdiff(predictors_list, dropped), collapse = " + ")
      f_reduced <- as.formula(paste(response, "~", fixed_reduced, "+", re_str))
      models[[dropped]] <- lme4::lmer(f_reduced, data = data_clean, REML = FALSE)
      anova_result[[dropped]] <- anova(models[[dropped]], model_all, test = "Chisq")
    }
  
    data_means <- get_data_means(data_clean, predictors_list, response)
  
  list(
    models = models,          # named by the predictor that was dropped
    model_all = model_all,
    anova_result = anova_result,
    data_means = data_means,
    used_data = data_clean
  )
}

predictors_list <- c("salinity", "nitrate")

#Pmax_____________________________________
#Inputs for ulva_ps/Pmax
ulva_ps_pmax_results <- compare_lmer_models(
  data = ulva_ps,
  response = "pmax",
  predictors_list,
  random_effects = c("plant_id", "rlc_order1")
)

# Access results
# Full model
summary(ulva_ps_pmax_results$model_all)

# Drop-one summaries (if you want to inspect them)
summary(ulva_ps_pmax_results$models[["salinity"]])
summary(ulva_ps_pmax_results$models[["nitrate"]])

# Does salinity effect Ulva Pmax?
ulva_ps_pmax_results$anova_result[["salinity"]]
# Does nitrate effect Ulva Pmax?
ulva_ps_pmax_results$anova_result[["nitrate"]]
check_model_fit(ulva_ps_pmax_results$model_all, terms = predictors_list)
ulva_ps_pmax_results$data_means

#inputs for hypnea_ps/Pmax
hypnea_ps_pmax_results <- compare_lmer_models(
  data = hypnea_ps,
  response = "pmax",
  predictors_list,
  random_effects = c("plant_id", "rlc_order1")
)

# Access results
summary(hypnea_ps_pmax_results$models[[1]]) # Model 1 summary
summary(hypnea_ps_pmax_results$models[[2]]) # Model 2 summary
summary(hypnea_ps_pmax_results$model_all) # Model 3 summary
hypnea_ps_pmax_results$anova_result[[1]]        # ANOVA did nitrate have an effect on Pmax?
hypnea_ps_pmax_results$anova_result[[2]]         # ANOVA did salinity have an effect on Pmax?
check_model_fit(hypnea_ps_pmax_results$model_all, terms = predictors_list)
hypnea_ps_pmax_results$data_means

#NPQmax______________________________
#inputs for ulva_ps/NPQmax
ulva_ps_npqmax_results <- compare_lmer_models(
  data = ulva_ps,
  response = "npq_max",
  predictors_list,
  random_effects = c("plant_id", "rlc_order", "illumination", "run_combo")
)

# Access results
summary(ulva_ps_npqmax_results$models[[1]]) # Model 1 summary
summary(ulva_ps_npqmax_results$models[[2]]) # Model 2 summary
summary(ulva_ps_npqmax_results$model_all) # Model 3 summary
ulva_ps_npqmax_results$anova[[1]]        # ANOVA did nitrate have an effect on NPQmax?
ulva_ps_npqmax_results$anova[[2]]        # ANOVA did salinity have an effect on NPQmax?
check_model_fit(ulva_ps_npqmax_results$model_all, terms = predictors_list)
ulva_ps_npqmax_results$data_means

#inputs for hypnea_ps/NPQmax
hypnea_ps_npqmax_results <- compare_lmer_models(
  data = hypnea_ps,
  response = "npq_max",
  predictors_list,
  random_effects = c("plant_id", "rlc_order", "illumination", "run_combo")
)

# Access results
summary(hypnea_ps_npqmax_results$models[[1]]) # Model 1 summary
summary(hypnea_ps_npqmax_results$models[[2]]) # Model 2 summary
summary(hypnea_ps_npqmax_results$model_all) # Model 3 summary
hypnea_ps_npqmax_results$anova[[1]]         # ANOVA did nitrate have an effect on NPQmax?
hypnea_ps_npqmax_results$anova[[2]]         # ANOVA did salinity have an effect on NPQmax?
check_model_fit(hypnea_ps_npqmax_results$model_all, terms = predictors_list)
hypnea_ps_npqmax_results$data_means

#deltaNPQ____________________________

#Inputs for ulva_ps/deltaNPQ
ulva_ps_delta_npq_results <- compare_lmer_models(
  data = ulva_ps,
  response = "deltaNPQ",
  predictors_list,
  random_effects = c("plantID")
)

# Access results
summary(ulva_ps_delta_npq_results$models[[1]]) # Model 1 summary
summary(ulva_ps_delta_npq_results$models[[2]]) # Model 2 summary
summary(ulva_ps_delta_npq_results$model_all) # Model 3 summary
ulva_ps_delta_npq_results$anova[[1]]         # ANOVA did nitrate have an effect on deltaNPQ?
ulva_ps_delta_npq_results$anova[[2]]         # ANOVA did salinity have an effect on deltaNPQ?
check_model_fit(ulva_ps_delta_npq_results$model_all, terms = predictors_list)
ulva_ps_delta_npq_results$data_means

#Inputs for hypnea_ps/deltaNPQ
hypnea_ps_delta_npq_results <- compare_lmer_models(
  data = hypnea_ps,
  response = "deltaNPQ",
  predictors_list,
  random_effects = c("plantID")
)

# Access results
summary(hypnea_ps_delta_npq_results$models[[1]]) # Model 1 summary
summary(hypnea_ps_delta_npq_results$models[[2]]) # Model 2 summary
summary(hypnea_ps_delta_npq_results$model_all) # Model 3 summary
hypnea_ps_delta_npq_results$anova[[1]]         # ANOVA did nitrate have an effect on deltaNPQ?
hypnea_ps_delta_npq_results$anova[[2]]         # ANOVA did salinity have an effect on deltaNPQ?
check_model_fit(hypnea_ps_delta_npq_results$model_all, terms = predictors_list)
hypnea_ps_delta_npq_results$data_means

#Ek____________________________________________
#Inputs for ulva_ps/Ek
ulva_ps_ek_results <- compare_lmer_models(
  data = ulva_ps,
  response = "ek.est",
  predictors_list,
  random_effects = c("plantID", "mean_mins28", "rlc_order")
)

# Access results
summary(ulva_ps_ek_results$models[[1]]) # Model 1 summary
summary(ulva_ps_ek_results$models[[2]]) # Model 2 summary
summary(ulva_ps_ek_results$model_all) # Model 3 summary
ulva_ps_ek_results$anova[[1]]         # ANOVA did nitrate have an effect on Ek?
ulva_ps_ek_results$anova[[2]]         # ANOVA did salinity have an effect on Ek?
check_model_fit(ulva_ps_ek_results$model_all, terms = predictors_list)
ulva_ps_ek_results$data_means

#Inputs for hypnea_ps/Ek
hypnea_ps_ek_results <- compare_lmer_models(
  data = hypnea_ps,
  response = "ek.est",
  predictors_list,
  random_effects = c("mean_mins28", "rlc_order")
)

# Access results
summary(hypnea_ps_ek_results$models[[1]]) # Model 1 summary
summary(hypnea_ps_ek_results$models[[2]]) # Model 2 summary
summary(hypnea_ps_ek_results$model_all) # Model 3 summary
hypnea_ps_ek_results$anova[[1]]         # ANOVA did nitrate have an effect on Ek?
hypnea_ps_ek_results$anova[[2]]         # ANOVA did salinity have an effect on Ek?
check_model_fit(hypnea_ps_ek_results$model_all, terms = predictors_list)
hypnea_ps_ek_results$data_means

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
plot(pmax_hypnea_ps$histo)
plot(pmax_hypnea_ps$plot)
plot(pmax_hypnea_ps$lin_regr)
ggsave("pmax_hypnea_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#NPQmax----------------

#inputs for ulva_ps/NPQmax
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
plot(npqmax_ulva_ps$histo)
plot(npqmax_ulva$plot)
plot(npqmax_ulva_ps$lin_regr)
ggsave("npqmax_ulva_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)



#inputs for hypnea_ps/NPQmax
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
plot(npqmax_hypnea_ps$histo)
plot(npqmax_hypnea_ps$plot)
plot(npqmax_hypnea_ps$lin_regr)
ggsave("npqmax_hypnea_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#deltaNPQ----------------

#inputs for ulva_ps/deltaNPQ
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
  title = "E - ulva_ps",
  title2 = "Chondria tumulosa ulva_ps --- ΔNPQ vs 9-Day Growth (%)",
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
plot(delta_npq_ulva_ps$histo)
plot(delta_npq_ulva_ps$plot)
plot(delta_npq_ulva_ps$lin_regr)
ggsave("delta_npq_ulva_ps.png", path = "wastewater/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#inputs for hypnea_ps/deltaNPQ
delta_npq_hypnea_ps <- raw_plots(
  data = hypnea_ps,
  response = deltaNPQ,
  response2 = d9_growth_percent,
  label = "hypnea_psstory",
  pretty_color = "maroon3",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "ΔNPQ",
  y = "Day 9 ΔNPQ",
  y2 = "9-Day Growth (%)",
  title = "F - hypnea_psstory",
  title2 = "Chondria tumulosa hypnea_psstory --- ΔNPQ vs 9-Day Growth (%)",
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
plot(delta_npq_hypnea_ps$histo)
plot(delta_npq_hypnea_ps$plot)
plot(delta_npq_hypnea_ps$lin_regr)
ggsave("delta_npq_hypnea_ps.png", path = "wastewater/plots/", 
       width = 7, height = 6, units = "in", dpi = 300, scale = 1)


#Ek----------------

#inputs for ulva_ps/Ek
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
plot(ek_ulva_ps$histo)
plot(ek_ulva$plot)
plot(ek_ulva_ps$lin_regr)
ggsave("ek_ulva_n6.png", path = "wastewater/plots/") 
       #width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#inputs for hypnea_ps/Ek
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
plot(ek_hypnea_ps$histo)
plot(ek_hypnea$plot)
plot(ek_hypnea_ps$lin_regr)
ggsave("ek_hypnea_n6.png", path = "wastewater/plots/") 
      # width = 7, height = 6, units = "in", dpi = 300, scale = 1)



#Linear regression Pmax vs Growth--------------------------------------------------------------------------------

#plot a regression between the photosynthetic independent variables of interest and growth rate


hypnea_ps_growth_pmax_plot <- hypnea_ps %>%
  ggplot(aes(x=pmax, 
             y=d9_growth_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + 
  theme_bw() + 
  labs(title = "Chondria tumulosa hypnea_psstory --- Pmax vs 9-Day Growth (%)", 
       x = "Pmax (μmols electrons m-2 s-1)", 
       y = "9-Day Growth (%)") + 
  stat_regline_equation(label.x = 10, label.y = 65) + stat_cor(label.x = 10, label.y = 60)
hypnea_ps_growth_pmax_plot


#--------------------NPQmax--------------------------






#summarize the means for NPQmax
ulva_ps %>% group_by(nitrate, salinity) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))
hypnea_ps %>% group_by(nitrate, salinity) %>% summarise_at(vars(maxNPQ_Ypoint1), list(mean = mean))


#Linear regression NPQ vs Growth and Apices
#plot a regression between the photosynthetic independent variables of interest and growth rate
#ulva_ps_growth_NPQmax_graph <- ggplot(ulva_ps_sub, aes(x=maxNPQ_Ypoint1, y=growth_rate_percent)) + 
geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Chondria tumulosa NPQmax vs 9-Day Growth (%)", x = " NPQmax (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
#ulva_ps_growth_NPQmax_graph



#-----------------delta NPQ------------------


#summarize the means for deltaNPQ
ulva_ps %>% group_by(treatment) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
ulva_ps %>% group_by(temp) %>% summarise_at(vars(deltaNPQ), list(mean = mean))
#ulva_ps %>% group_by(treatment, rlc_day) %>% summarise_at(vars(pmax), list(mean = mean))

#Linear regression deltaNPQ vs Growth and Apices

#plot a regression between the photosynthetic independent variables of interest and growth rate
ulva_ps_growth_deltaNPQ_graph <- ggplot(ulva_ps_sub, aes(x=deltaNPQ, y=growth_rate_percent)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Chondria tumulosa Delta NPQ vs 9-Day Growth (%)", x = " Delta NPQ (rel. units)", 
       y = "9-Day Growth (%)") + stat_regline_equation(label.x = 0.25, label.y = 38) + stat_cor(label.x = 0.25, label.y = 40)
ulva_ps_growth_deltaNPQ_graph


