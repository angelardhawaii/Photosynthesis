# Script to run model for Hsat and DSPI analysis for Hypnea and Ulva wastewater data
# This script pulls in the combined dataset for the 2021, 2023 SGD and 2025 wastewater data
# sal_bin predictor (35 ppt = 1, 28, 27 ppt = 2, 18-22 ppt = 3, 11-15 ppt = 4)
# By Angela Richards Donà
# Date: 3/4/2026

#load libraries
# data wrangling
library(tidyverse) #includes ggplot2, readr, forcats, tibble, purrr, tidyr, dplyr, stringr
library(hms)
library(lubridate)
# mixed models
library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(multcomp)
library(ggeffects)
library(performance)
library(ggpubr)
library(sjPlot)
# optional for plots
library(RColorBrewer)

#load this file for normalized to quantum efficiency of photosynthesis per Silsbe and Kromkamp
ww_hsat_lmm <- read_csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/ww_combo_hsat_dspi.csv")

# remove individual that died during experiment
ww_hsat_lmm <- ww_hsat_lmm %>%
  filter(id != "ul09")

glimpse(ww_hsat_lmm)

TZ <- "Pacific/Honolulu"

ww_hsat_lmm <- ww_hsat_lmm %>%
  mutate(
    start_date = as.Date(start_date),
    start_time = hms::as_hms(start_time),
    end_time   = hms::as_hms(end_time)
  )

# Make unique_id column for troubleshooting individuals
ww_hsat_lmm <- ww_hsat_lmm %>%
  mutate(unique_id = paste(start_date, id, sep = "_"))

# assign list of variables as factor
factor_vars <- c("sal_bin", "temp", "nitrate", "run_combo", "salinity")

ww_hsat_lmm <- ww_hsat_lmm %>%
  mutate(across(all_of(factor_vars), as.factor))

#assign new column for chronological ranking of individuals in each day_group
#use ranking to make smaller groups of ~15 minutes for rlc_order
ww_hsat_lmm <- ww_hsat_lmm %>%
  group_by(start_date, rlc_day) %>%
  arrange(end_time, .by_group = TRUE) %>%
  mutate(
    rank = dplyr::row_number(),                 # stable rank by time
    rlc_order1 = floor((rank - 1) / 3) + 1      # 3 steps ≈ 15 minutes
  ) %>%
  ungroup()
# add a year column for subsetting later
ww_hsat_lmm <- ww_hsat_lmm %>%
  mutate(year = lubridate::year(start_date))

# Convert DSPI from micromols to mols m-2 (dspi_day_mean computed below after run_days is known)
ww_hsat_lmm <- ww_hsat_lmm %>%
  mutate(
    dspi_mol = round(dspi_total / 1e6, 2)
  )
ww_hsat_lmm <- ww_hsat_lmm %>%
  mutate(
    start_dt = lubridate::as_datetime(start_date, tz = TZ) + start_time,
    end_dt   = lubridate::as_datetime(start_date + days(8), tz = TZ) + end_time,  # all experiments ran day 1–9, so 8 elapsed days
    run_hours = as.numeric(difftime(end_dt, start_dt, units = "hours")),
    run_days  = run_hours / 24
  ) %>%
  mutate(
    # mean per day + round to 2 decimals
    hsat_day_mean = round(hsat_total / run_days, 2),
    relhsat_day_mean = round(relhsat_total / run_days, 2),
    dspi_day_mean = round(dspi_mol / run_days, 2)
  )

# Check data before moving on
ww_hsat_lmm %>%
  summarise(
    min_days = min(run_days, na.rm = TRUE),
    max_days = max(run_days, na.rm = TRUE),
    n_bad = sum(is.na(run_days) | run_days <= 0, na.rm = TRUE)
  )

# Join growth rates (growth_d9) from combined_data_all.
# run_combo numbering differs between the two files (different schemes), so
# it cannot be used as a join key. id alone is not unique either (same IDs
# reused across years and across multiple runs within the same year).
# The only collision-free key is id + end_date, where end_date = start_date + 8 days
# in ww_combo_hsat_dspi, which matches the 'date' column in combined_data_all.
growth_key <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/combined_data_all.csv") %>%
  dplyr::select(id, date, growth_d9) %>%
  mutate(id      = tolower(trimws(id)),
         end_date = as.Date(date))

ww_hsat_lmm <- ww_hsat_lmm %>%
  mutate(end_date = as.Date(start_date) + days(8)) %>%
  left_join(growth_key, by = c("id", "end_date")) %>%
  mutate(run_combo = as.factor(run_combo))

# new dataset with overlaps from 2021, 2022, 2023 and all of 2025
# this keeps all of 2025, and only the 53 and 80 nitrate treatments from 2021-2023, which are the only ones that overlap with 2025.
# It also removes the 30 degree temp treatment from 2021-2023, which is not present in 2025.
sgd_ww_hsat <- ww_hsat_lmm %>%
  filter(
    # keep:
    (
      year %in% c(2021, 2022) & 
        nitrate == "80" & 
        temp != 30
    ) |
      year == 2025
  )

# subset the species — combined 2021/2022 + 2025 dataset
ulva_ww_hsat <- subset(sgd_ww_hsat, species == "u")

hypnea_ww_hsat <- subset(sgd_ww_hsat, species == "h")

# 2025-only subsets for comparison
# Note: all 2025 data is at 22 degrees C, so temp cannot be used as a predictor;
# models use nitrate + salinity only.
ulva_2025   <- subset(ww_hsat_lmm, species == "u" & year == 2025)
hypnea_2025 <- subset(ww_hsat_lmm, species == "h" & year == 2025)

#Check the data
xtabs(~ salinity + temp + year, data = ulva_ww_hsat)
xtabs(~ salinity + temp + year, data = hypnea_ww_hsat)
#MODEL FUNCTION__________________________________________________________________

check_model_fit <- function(model) {
  hist(resid(model))
  plot(resid(model) ~ fitted(model))
  qqnorm(resid(model))
  qqline(resid(model))
  print(performance::check_model(model))
  print(MuMIn::r.squaredGLMM(model))
}

#get predictor means
get_data_means <- function(data, predictors_list, response){
  data %>%
    group_by(across(all_of(predictors_list))) %>%
    summarise(!!paste0(response, "_mean") := mean(.data[[response]], na.rm = TRUE),
              .groups = "drop")
}

#Make a function to run the models for likelihood ratio tests where REML must be FALSE
#construct null model to perform likelihood ratio test REML must be FALSE

compare_lmer_models <- function(data, response, predictors_list, random_effects) {
  # remove NA/NaN/Inf ONLY for the response variable
  data_clean <- data %>%
    filter(is.finite(.data[[response]])) %>%              # drop NA / NaN / Inf in response
    tidyr::drop_na(all_of(predictors_list)) %>%           # drop rows missing predictor values
    droplevels()

  # Build random-effects string
  re_str <- paste(paste0("(1 | ", random_effects, ")"), collapse = " + ")

  # Optimizer control: bobyqa with more iterations suppresses non-convergence
  # warnings common with small numbers of random effect levels
  lmer_ctrl <- lme4::lmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 2e5))

  # Full model (REML = FALSE required for LRTs)
  fixed_all <- paste(predictors_list, collapse = " + ")
  f_full <- as.formula(paste(response, "~", fixed_all, "+", re_str))
  model_all <- lme4::lmer(f_full, data = data_clean, REML = FALSE, control = lmer_ctrl)

  # Refit full model with REML = TRUE for reporting coefficients and diagnostics
  model_all_reml <- lme4::lmer(f_full, data = data_clean, REML = TRUE, control = lmer_ctrl)

  # Drop-one models + LRTs
  models <- list()
  anova_result <- list()

  for (dropped in predictors_list) {
    fixed_reduced <- paste(setdiff(predictors_list, dropped), collapse = " + ")
    f_reduced <- as.formula(paste(response, "~", fixed_reduced, "+", re_str))
    models[[dropped]] <- lme4::lmer(f_reduced, data = data_clean, REML = FALSE,
                                    control = lmer_ctrl)
    anova_result[[dropped]] <- anova(models[[dropped]], model_all, test = "Chisq")
  }

  data_means <- get_data_means(data_clean, predictors_list, response)

  list(
    models = models,          # named by the predictor that was dropped
    model_all = model_all,    # ML fit — use only for LRTs
    model_all_reml = model_all_reml,  # REML fit — use for reporting and diagnostics
    anova_result = anova_result,
    data_means = data_means,
    used_data = data_clean
  )
}

predictors_list <- c("nitrate", "temp")
predictors_list_species <- c("nitrate", "temp", "species")

# MODEL SELECTION NOTE — combined 2021+2022 and 2025 dataset:
# Full fixed effects tested: salinity + nitrate + temp
# Salinity was dropped from all combined-dataset models: including it caused rank
# deficiency in the fixed-effect model matrix across every response variable and
# both species (confirmed by "fixed-effect model matrix is rank deficient" warning).
# Rank deficiency arises because salinity treatment structure is collinear with
# temperature across the two datasets — certain salinity-temperature combinations
# exist only in 2021 or only in 2025, making their effects inseparable.
# For DSPI models, rlc_order1 was also dropped from random effects to resolve
# singular fit (variance component collapsing to zero).
# Final fixed effects retained: nitrate + temp
# Salinity IS included as a predictor in the 2025-only models below, where it
# varies independently of temperature.

#Hsat_____________________________________
#Inputs for ulva_ps/Hsat
ulva_ww_hsat_results <- compare_lmer_models(
  data = ulva_ww_hsat,
  response = "hsat_day_mean",
  predictors_list,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)

# Access results
# Full model
summary(ulva_ww_hsat_results$model_all_reml)

# Does nitrate affect Ulva hsat?
ulva_ww_hsat_results$anova_result[["nitrate"]]
#Does temp affect Ulva hsat?
ulva_ww_hsat_results$anova_result[["temp"]]
check_model_fit(ulva_ww_hsat_results$model_all_reml)
ulva_ww_hsat_results$data_means


#inputs for hypnea_ww/hsat
hypnea_ww_hsat_results <- compare_lmer_models(
  data = hypnea_ww_hsat,
  response = "hsat_day_mean",
  predictors_list,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)

# Access results ALL 3 PREDICTORS significant
summary(hypnea_ww_hsat_results$model_all_reml)

# Does nitrate affect Hypnea hsat?
hypnea_ww_hsat_results$anova_result[["nitrate"]]
#Does temp affect Hypnea hsat?
hypnea_ww_hsat_results$anova_result[["temp"]]

check_model_fit(hypnea_ww_hsat_results$model_all_reml)
hypnea_ww_hsat_results$data_means

#Inputs for Hsat species comparison
ww_hsat_species_results <- compare_lmer_models(
  data = sgd_ww_hsat,
  response = "hsat_day_mean",
  predictors_list_species,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)

# Access results
# Full model
summary(ww_hsat_species_results$model_all_reml)

# Does species affect hsat?
ww_hsat_species_results$anova_result[["species"]]

check_model_fit(ww_hsat_species_results$model_all_reml)
ww_hsat_species_results$data_means

#relHsat____________________________________

# Inputs for Ulva relHsat
ulva_ww_relhsat_results <- compare_lmer_models(
  data = ulva_ww_hsat,
  response = "relhsat_day_mean",
  predictors_list,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)

# Access results
# Full model
summary(ulva_ww_relhsat_results$model_all_reml)

# Does nitrate affect Ulva relhsat?
ulva_ww_relhsat_results$anova_result[["nitrate"]]
#Does temp affect Ulva relhsat?
ulva_ww_relhsat_results$anova_result[["temp"]]
check_model_fit(ulva_ww_relhsat_results$model_all_reml)
ulva_ww_relhsat_results$data_means

# Inputs for Hypnea_relhsat
hypnea_ww_relhsat_results <- compare_lmer_models(
  data = hypnea_ww_hsat,
  response = "relhsat_day_mean",
  predictors_list,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)

# Access results
summary(hypnea_ww_relhsat_results$model_all_reml)

# Does nitrate affect Hypnea relhsat?
hypnea_ww_relhsat_results$anova_result[["nitrate"]]
#Does temp affect Hypnea relhsat?
hypnea_ww_relhsat_results$anova_result[["temp"]]

check_model_fit(hypnea_ww_relhsat_results$model_all_reml)
hypnea_ww_relhsat_results$data_means

# Inputs for relHsat species comparison
ww_relhsat_species_results <- compare_lmer_models(
  data = sgd_ww_hsat,
  response = "relhsat_day_mean",
  predictors_list_species,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)

# Access results
# Full model
summary(ww_relhsat_species_results$model_all_reml)

# Drop-one summaries (if you want to inspect them)
summary(ww_relhsat_species_results$models[["species"]])

# Does species affect relhsat?
ww_relhsat_species_results$anova_result[["species"]]

check_model_fit(ww_relhsat_species_results$model_all_reml)
ww_relhsat_species_results$data_means

#DSPI_____________________________________
#Inputs for Ulva DSPI
ulva_ww_dspi_results <- compare_lmer_models(
  data = ulva_ww_hsat,
  response = "dspi_day_mean",
  predictors_list,
  random_effects = c("plant_id", "run_combo")
)

# Access results
summary(ulva_ww_dspi_results$model_all_reml)

# Does nitrate affect Ulva dspi?
ulva_ww_dspi_results$anova_result[["nitrate"]]
#Does temp affect Ulva dspi?
ulva_ww_dspi_results$anova_result[["temp"]]
check_model_fit(ulva_ww_dspi_results$model_all_reml)
ulva_ww_dspi_results$data_means

#inputs for hypnea_ww/dspi
hypnea_ww_dspi_results <- compare_lmer_models(
  data = hypnea_ww_hsat,
  response = "dspi_day_mean",
  predictors_list,
  random_effects = c("plant_id", "run_combo")
)

# Access results
summary(hypnea_ww_dspi_results$model_all_reml)

# Does nitrate affect Hypnea dspi?
hypnea_ww_dspi_results$anova_result[["nitrate"]]
#Does temp affect Hypnea dspi?
hypnea_ww_dspi_results$anova_result[["temp"]]

check_model_fit(hypnea_ww_dspi_results$model_all_reml)
hypnea_ww_dspi_results$data_means

#Inputs for dspi species comparison
ww_dspi_species_results <- compare_lmer_models(
  data = sgd_ww_hsat,
  response = "dspi_day_mean",
  predictors_list_species,
  random_effects = c("plant_id", "run_combo")
)

# Access results
# Full model
summary(ww_dspi_species_results$model_all_reml)

# Drop-one summaries (if you want to inspect them)
summary(ww_dspi_species_results$models[["species"]])

# Does species affect dspi?
ww_dspi_species_results$anova_result[["species"]]

check_model_fit(ww_dspi_species_results$model_all_reml)
ww_dspi_species_results$data_means

# Model Summary Tables_____________________________________________________

# This one is like the old printouts
print_model <- function(model, species, response) {
  sjPlot::tab_model(
    model,
    title = paste(species, "-", response),
    show.intercept = TRUE,
    show.se = TRUE,
    show.stat = TRUE,
    show.df = TRUE
  )
}
print_model(ulva_ww_hsat_results$model_all_reml, "Ulva", "Hsat day mean")
print_model(hypnea_ww_hsat_results$model_all_reml, "Hypnea", "Hsat day mean")

print_model(ulva_ww_relhsat_results$model_all_reml, "Ulva", "Relative Hsat")
print_model(hypnea_ww_relhsat_results$model_all_reml, "Hypnea", "Relative Hsat")

print_model(ulva_ww_dspi_results$model_all_reml, "Ulva", "DSPI day mean")
print_model(hypnea_ww_dspi_results$model_all_reml, "Hypnea", "DSPI day mean")

#HISTOGRAMS and PLOTS_____________________________________
#function for raw data plots
raw_plots <- function(data, response, response2, label, pretty_color, aescolor, x, x2, y, y2, title, title2, subtitle, ylim1, ylim2, 
                      color1, color2, color3, yint, vjust_t, hjust_t, vjust_s, hjust_s, labels1) {
  histo <- data %>%
    ggplot(aes(x = {{response}})) +
    geom_histogram(fill = pretty_color, color = "black", linewidth = 0.25, alpha = 0.85) +
    labs(title = label) +
    theme_bw()
  
  plot <- data %>% 
    ggplot(aes(x = nitrate, y = {{response}})) + 
    geom_boxplot(size=0.5) + 
    geom_point(alpha = 0.75, size = 3, aes(color = {{aescolor}}), position = "jitter", show.legend = TRUE) + 
    labs(x = x, y = y, title = title, subtitle = subtitle) + 
    #scale_x_discrete(labels = c("0.5μmolN", "2μmolN", "4μmolN", "8μmolN")) + 
    ylim(ylim1, ylim2) + stat_mean() + 
    scale_color_manual(values = c(color1, color2, color3)) +
    scale_x_discrete(labels = labels1) +
    geom_hline(yintercept=yint, color = "red", linewidth = 0.5, alpha = 0.5) +
    theme_bw() +
    theme(legend.position = c(0.90,0.80), plot.title = element_text(face = "italic", vjust = vjust_t, hjust = hjust_t), 
          plot.subtitle = element_text(face = "bold", size = 14, vjust = vjust_s, hjust = hjust_s))
  
  return(list(
    histo = histo,
    plot = plot)) 
}

#INPUTS for histo and plots-----------------

#Hsat ------------------------

#Inputs for Ulva/Hsat
hsat_ulva <- raw_plots(
  data = ulva_ww_hsat,
  response = hsat_day_mean,
  response2 = growth,
  label = "Ulva lactuca",
  pretty_color = "goldenrod1",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "Hsat (min/day)",
  y = "Day 9 Hsat (min/day)",
  y2 = "9-Day Growth (%)",
  title = "A - Ulva lactuca",
  title2 = "Ulva lactuca --- Hsat (min/day) vs 9-Day Growth (%)",
  subtitle = "Hsat",
  ylim1 = 100, 
  ylim2 = 1000,
  color1 = "forestgreen",
  color2 = "darkolivegreen3",
  color3 = "olivedrab4",
  yint = 400,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -10,
  hjust_s = 0.95,
  labels1 = c("80 μmol", 
              "245 μmol",
              "748 μmol",
              "2287 μmol")
)
plot(hsat_ulva$histo)
plot(hsat_ulva$plot)
#ggsave("hsat_ulva_sgdww.png", path = "wastewater/plots/", 
#      width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#inputs for hypnea/hsat
hsat_hypnea <- raw_plots(
  data = hypnea_ww_hsat,
  response = hsat_day_mean,
  response2 = growth,
  label = "Hypnea musciformis",
  pretty_color = "maroon4",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "Hsat (min/day)",
  y = "Day 9 Hsat (min/day)",
  y2 = "9-Day Growth (%)",
  title = "B - Hypnea musciformis",
  title2 = "Hypnea musciformis --- Hsat vs 9-Day Growth (%)",
  subtitle = "Hsat",
  ylim1 = 100, 
  ylim2 = 1000,
  color1 = "deeppink4",
  color2 = "deeppink3",
  color3 = "hotpink",
  yint = 20,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -10,
  hjust_s = 0.95,
  labels1 = c("80 μmol", 
              "245 μmol",
              "748 μmol",
              "2287 μmol")
)
plot(hsat_hypnea$histo)
plot(hsat_hypnea$plot)
#ggsave("hsat_hypnea_sgdww.png", path = "wastewater/plots/",
#      width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#Hsat ------------------------

#Inputs for Ulva/Hsat
dspi_ulva <- raw_plots(
  data = ulva_ww_hsat,
  response = dspi_day_mean,
  response2 = growth,
  label = "Ulva lactuca",
  pretty_color = "goldenrod1",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "DSPI (mol/day)",
  y = "Day 9 DSPI (mol/day)",
  y2 = "9-Day Growth (%)",
  title = "A - Ulva lactuca",
  title2 = "Ulva lactuca --- DSPI (mol/day) vs 9-Day Growth (%)",
  subtitle = "DSPI",
  ylim1 = 0, 
  ylim2 = 5,
  color1 = "forestgreen",
  color2 = "darkolivegreen3",
  color3 = "olivedrab4",
  yint = 2,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -10,
  hjust_s = 0.95,
  labels1 = c("80 μmol", 
              "245 μmol",
              "748 μmol",
              "2287 μmol")
)
plot(dspi_ulva$histo)
plot(dspi_ulva$plot)
#ggsave("hsat_ulva_sgdww.png", path = "wastewater/plots/", 
#      width = 7, height = 6, units = "in", dpi = 300, scale = 1)

#inputs for hypnea/hsat
dspi_hypnea <- raw_plots(
  data = hypnea_ww_hsat,
  response = dspi_day_mean,
  response2 = growth,
  label = "Hypnea musciformis",
  pretty_color = "maroon4",
  aescolor = salinity,
  x = "Nitrate",
  x2 = "DSPI (mol/day)",
  y = "Day 9 DSPI (mol/day)",
  y2 = "9-Day Growth (%)",
  title = "B - Hypnea musciformis",
  title2 = "Hypnea musciformis --- DSPI vs 9-Day Growth (%)",
  subtitle = "DSPI",
  ylim1 = 0, 
  ylim2 = 5,
  color1 = "deeppink4",
  color2 = "deeppink3",
  color3 = "hotpink",
  yint = 2,
  vjust_t = -15,
  hjust_t = 0.05,
  vjust_s = -10,
  hjust_s = 0.95,
  labels1 = c("80 μmol", 
              "245 μmol",
              "748 μmol",
              "2287 μmol")
)
plot(dspi_hypnea$histo)
plot(dspi_hypnea$plot)
#ggsave("hsat_hypnea_sgdww.png", path = "wastewater/plots/",
#      width = 7, height = 6, units = "in", dpi = 300, scale = 1)

# 2025-ONLY COMPARISON MODELS ________________________________________________
# Predictors: nitrate + salinity only (temp is constant at 22C in 2025)
predictors_list_2025        <- c("nitrate", "salinity")
predictors_list_2025_species <- c("nitrate", "salinity", "species")

# Hsat
ulva_2025_hsat_results <- compare_lmer_models(
  data = ulva_2025,
  response = "hsat_day_mean",
  predictors_list_2025,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)
summary(ulva_2025_hsat_results$model_all_reml)
ulva_2025_hsat_results$anova_result[["nitrate"]]
ulva_2025_hsat_results$anova_result[["salinity"]]
check_model_fit(ulva_2025_hsat_results$model_all_reml)
ulva_2025_hsat_results$data_means

hypnea_2025_hsat_results <- compare_lmer_models(
  data = hypnea_2025,
  response = "hsat_day_mean",
  predictors_list_2025,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)
summary(hypnea_2025_hsat_results$model_all_reml)
hypnea_2025_hsat_results$anova_result[["nitrate"]]
hypnea_2025_hsat_results$anova_result[["salinity"]]
check_model_fit(hypnea_2025_hsat_results$model_all_reml)
hypnea_2025_hsat_results$data_means

# relHsat
ulva_2025_relhsat_results <- compare_lmer_models(
  data = ulva_2025,
  response = "relhsat_day_mean",
  predictors_list_2025,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)
summary(ulva_2025_relhsat_results$model_all_reml)
ulva_2025_relhsat_results$anova_result[["nitrate"]]
ulva_2025_relhsat_results$anova_result[["salinity"]]
check_model_fit(ulva_2025_relhsat_results$model_all_reml)
ulva_2025_relhsat_results$data_means

hypnea_2025_relhsat_results <- compare_lmer_models(
  data = hypnea_2025,
  response = "relhsat_day_mean",
  predictors_list_2025,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)
summary(hypnea_2025_relhsat_results$model_all_reml)
hypnea_2025_relhsat_results$anova_result[["nitrate"]]
hypnea_2025_relhsat_results$anova_result[["salinity"]]
check_model_fit(hypnea_2025_relhsat_results$model_all_reml)
hypnea_2025_relhsat_results$data_means

# DSPI
ulva_2025_dspi_results <- compare_lmer_models(
  data = ulva_2025,
  response = "dspi_day_mean",
  predictors_list_2025,
  random_effects = c("run_combo")
)
summary(ulva_2025_dspi_results$model_all_reml)
ulva_2025_dspi_results$anova_result[["nitrate"]]
ulva_2025_dspi_results$anova_result[["salinity"]]
check_model_fit(ulva_2025_dspi_results$model_all_reml)
ulva_2025_dspi_results$data_means

hypnea_2025_dspi_results <- compare_lmer_models(
  data = hypnea_2025,
  response = "dspi_day_mean",
  predictors_list_2025,
  random_effects = c("run_combo")
)
summary(hypnea_2025_dspi_results$model_all_reml)
hypnea_2025_dspi_results$anova_result[["nitrate"]]
hypnea_2025_dspi_results$anova_result[["salinity"]]
check_model_fit(hypnea_2025_dspi_results$model_all_reml)
hypnea_2025_dspi_results$data_means

# Species comparison — 2025 only
ww_2025_hsat_species_results <- compare_lmer_models(
  data = subset(ww_hsat_lmm, year == 2025),
  response = "hsat_day_mean",
  predictors_list_2025_species,
  random_effects = c("plant_id", "run_combo", "rlc_order1")
)
ww_2025_hsat_species_results$anova_result[["species"]]

ww_2025_dspi_species_results <- compare_lmer_models(
  data = subset(ww_hsat_lmm, year == 2025),
  response = "dspi_day_mean",
  predictors_list_2025_species,
  random_effects = c("run_combo")
)
ww_2025_dspi_species_results$anova_result[["species"]]

# Model summary tables — 2025 only
print_model(ulva_2025_hsat_results$model_all_reml,   "Ulva 2025",   "Hsat day mean")
print_model(hypnea_2025_hsat_results$model_all_reml, "Hypnea 2025", "Hsat day mean")
print_model(ulva_2025_dspi_results$model_all_reml,   "Ulva 2025",   "DSPI day mean")
print_model(hypnea_2025_dspi_results$model_all_reml, "Hypnea 2025", "DSPI day mean")

# DSPI vs Growth linear regression ____________________________________________
# DSPI (mean daily, mol O2 m-2 day-1) as predictor of 9-day growth (%)
# Uses the full dataset (ww_hsat_lmm) — all years and all nitrate levels —
# so the regression is not limited to the sgd_ww_hsat LMM subset.
# Points coloured by nitrate treatment; single regression line across all treatments.
# Annotation shows equation with slope, R2, and p-value.

dspi_growth_ulva <- ulva_ww_hsat %>%
  filter(!is.na(dspi_day_mean), !is.na(growth_d9)) %>%
  ggplot(aes(x = dspi_day_mean, y = growth_d9)) +
  geom_point(aes(color = nitrate), alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    label.x.npc = "left", label.y.npc = 0.97, size = 3.8
  ) +
  stat_cor(
    aes(label = ..p.label..),
    label.x.npc = "left", label.y.npc = 0.90, size = 3.8
  ) +
  scale_color_viridis_d(option = "viridis", direction = -1) +
  labs(
    x     = "Mean daily DSPI (mol O2 m-2 day-1)",
    y     = "9-day growth (%)",
    title = "Ulva lactuca",
    color = "Nitrate (μmol)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic", size = 13),
        axis.title = element_text(size = 12))
plot(dspi_growth_ulva)

dspi_growth_hypnea <- hypnea_ww_hsat %>%
  filter(!is.na(dspi_day_mean), !is.na(growth_d9)) %>%
  ggplot(aes(x = dspi_day_mean, y = growth_d9)) +
  geom_point(aes(color = nitrate), alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    label.x.npc = "left", label.y.npc = 0.97, size = 3.8
  ) +
  stat_cor(
    aes(label = ..p.label..),
    label.x.npc = "left", label.y.npc = 0.90, size = 3.8
  ) +
  scale_color_viridis_d(option = "magma", direction = -1) +
  labs(
    x     = "Mean daily DSPI (mol O2 m-2 day-1)",
    y     = "9-day growth (%)",
    title = "Hypnea musciformis",
    color = "Nitrate (\u00b5mol)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic", size = 13),
        axis.title = element_text(size = 12))
plot(dspi_growth_hypnea)

# Both species side by side
dspi_growth_both <- ggarrange(
  dspi_growth_ulva, dspi_growth_hypnea,
  ncol = 2, nrow = 1,
  labels = c("A", "B")
)
plot(dspi_growth_both)

# ggsave("dspi_vs_growth_ulva.png",   plot = dspi_growth_ulva,   path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots", width = 6, height = 5, dpi = 300)
# ggsave("dspi_vs_growth_hypnea.png", plot = dspi_growth_hypnea, path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots", width = 6, height = 5, dpi = 300)
# ggsave("dspi_vs_growth_both.png",   plot = dspi_growth_both,   path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots", width = 11, height = 5, dpi = 300)
