# Experimenting with predictive modeling
# Boosted regression trees
# By Angela Richards Donà
# July 09, 2025

#load libraries
library(randomForest)
library(caret)
library(e1071)
library(dplyr)
library(lubridate)
library(gbm)
library(dismo)
library(ggplot2)
library(gridExtra)
library(janitor)
library(stringr)

#Use combined dataset transformed in growth_combine_waste_sgd.R script
combo_waste_sgd <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/combined_data.csv")

combo_waste_sgd$species <- as.factor(combo_waste_sgd$species)
combo_waste_sgd$salinity <- as.factor(combo_waste_sgd$salinity)
combo_waste_sgd$nitrate <- round(as.numeric(combo_waste_sgd$n_g_m3_2d), digits = 3)
combo_waste_sgd$phosphate <- round(as.numeric(combo_waste_sgd$p_g_m3_2d), digits = 3)
combo_waste_sgd$temp <- as.factor(combo_waste_sgd$temp)

combined_growth <- combined_growth %>%
  mutate(run_combined = paste(year, run, sep = "_")) %>%  # Make a unique label per year and run
  mutate(run_id = dense_rank(run_combined)) %>%  # Assign unique numbers from 1 to 10
  select(-run_combined)

#MOVED TO COMBINE_________________________________

#Assign new run values to avoid issues of repeats and to get true chronology
combo_waste_sgd <- combo_waste_sgd %>%
  mutate(run_combo = case_when(
    date == as.Date("2021-09-28") ~ 1,
    date == as.Date("2021-09-29") ~ 1,
    date == as.Date("2021-10-08") ~ 2,
    date == as.Date("2021-10-09") ~ 2,
    date == as.Date("2021-10-19") ~ 3,
    date == as.Date("2021-10-29") ~ 4,
    date == as.Date("2021-11-09") ~ 5,
    date == as.Date("2021-11-12") ~ 5,
    date == as.Date("2022-02-19") ~ 6,
    date == as.Date("2022-03-01") ~ 7,
    date == as.Date("2022-04-29") ~ 8,
    date == as.Date("2025-01-31") ~ 9,
    date == as.Date("2025-02-25") ~ 10,
    date == as.Date("2025-03-14") ~ 11,
    date == as.Date("2025-03-18") ~ 11,
    date == as.Date("2025-03-25") ~ 12,
    date == as.Date("2025-03-28") ~ 12
  ))
combo_waste_sgd$run_combo <- as.factor(combo_waste_sgd$run_combo)

combo_waste_sgd <- combo_waste_sgd %>%
  mutate(t_combo = case_when(
    year == 2021 & treatment == 1 ~ 2,
    year == 2021 & treatment == 2 ~ 3,
    year == 2021 & treatment == 3 ~ 5,
    year == 2021 & treatment == 4 ~ 6,
    year == 2022 & treatment == 2.5 ~ 4,
    year == 2022 & treatment == 0 ~ 1,
    year == 2025 & treatment == 1 ~ 7,
    year == 2025 & treatment == 3 ~ 8,
    year == 2025 & treatment == 5 ~ 9,
    year == 2025 & treatment == 7 ~ 10
  ))

combo_waste_sgd$t_combo <- as.factor(combo_waste_sgd$t_combo)

combo_waste_sgd <- combo_waste_sgd %>%
  mutate(
    species = case_when(
      species %in% c("Hm", "h") ~ "h",
      species %in% c("Ul", "u") ~ "u"
    )
  )

#subset to include only certain treatments - those with some overlap
combo_waste_sgd <- combo_waste_sgd %>%
  filter(t_combo %in% 5:10)

glimpse(combo_waste_sgd)

#subset by species
combo_waste_sgd_u <- combo_waste_sgd %>% 
  filter(species == "u")

combo_waste_sgd_h <- combo_waste_sgd %>% 
  filter(species == "h")



#Run brt
#Combined data Ulva
#plot histogram
histo_u <- combo_waste_sgd_u %>%
  ggplot(aes(growth_d9)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_u 

#transform for better normality
combo_waste_sgd_u <- combo_waste_sgd_u %>%
  mutate(
    sqrt_growth = sqrt(growth_d9 + 57 + 0.01)  # shift + small constant
  )

histo_u_sqrt <- combo_waste_sgd_u %>%
  ggplot(aes(sqrt_growth)) +
  geom_histogram(binwidth=1, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_u_sqrt 

set.seed(42)

brt_combo_ulva <- gbm.step(
  data = combo_waste_sgd_u,
  gbm.x = c("salinity", "nitrate", "phosphate", "temp"),  # predictors
  gbm.y = "sqrt_growth", # response
  family = "gaussian",
  tree.complexity = 3,
  learning.rate = 0.001,
  bag.fraction = 0.5,
  n.folds = 5
)

#Ask for the contributions of each predictor
brt_combo_ulva$contributions

#Plotting these contributions. Adapt the number of plots to the number of predictors.
gbm.plot(brt_combo_ulva, n.plots=4, write.title = F,
         common.scale = TRUE, y.label="Fitted function",
         show.contrib = TRUE)

#You can also obtain the values of each observation on the fitted function by using
gbm.plot.fits(brt_combo_ulva)



#Combined data Hypnea
#Run brt
#plot histogram
histo_h <- combo_waste_sgd_h %>%
  ggplot(aes(growth_d9)) +
  geom_histogram(binwidth=5, fill = "deeppink4", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_h 

#transform for better normality
combo_waste_sgd_h <- combo_waste_sgd_h %>%
  mutate(
    sqrt_growth_h = sqrt(growth_d9 + 88 + 0.01)  # shift + small constant
  )

histo_h_sqrt <- combo_waste_sgd_h %>%
  ggplot(aes(sqrt_growth_h)) +
  geom_histogram(binwidth=1, fill = "deeppink4", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_h_sqrt 

set.seed(42)

brt_combo_hypnea <- gbm.step(
  data = combo_waste_sgd_h,
  gbm.x = c("salinity", "nitrate", "phosphate", "temp"),  # predictors
  gbm.y = "sqrt_growth_h", # response
  family = "gaussian",
  tree.complexity = 3,
  learning.rate = 0.011,
  bag.fraction = 0.5,
  n.folds = 5
)

#Ask for the contributions of each predictor
brt_combo_hypnea$contributions

#Plotting these contributions. Adapt the number of plots to the number of predictors.
gbm.plot(brt_combo_hypnea, n.plots=4, write.title = F,
         common.scale = TRUE, y.label="Fitted function",
         show.contrib = TRUE)

#You can also obtain the values of each observation on the fitted function by using
gbm.plot.fits(brt_combo_hypnea)



#______________________________________________________________________________
#OLD DATA from 2022 WRR Paper
#Try to reproduce the results Jade got from the old data I provided her
wg_old <- read.csv("/Users/angela/src/Photosynthesis/data/hyp_ulv_old/ulva_hyp_growth_all_051122.csv")

wg_old <- clean_names(wg_old) #clean up ugly names
wg_old$treatment <- as.factor(wg_old$treatment)
wg_old$species <- as.factor(wg_old$species)
wg_old$salinity <- as.factor(wg_old$salinity_ppt)
wg_old$run <- as.factor(wg_old$run)
wg_old$nitrate <- round(as.numeric(wg_old$n_g_m3_2d), digits = 3)
wg_old$phosphate <- round(as.numeric(wg_old$p_g_m3_2d), digits = 3)
wg_old$temperature_c <- as.factor(wg_old$temperature_c)

wg_old <- wg_old %>%
  mutate(date = mdy(date))

wg_old <- wg_old %>%
  mutate(unique_id = paste(date, id, sep = "_")) 

wg_old <- wg_old %>%
  dplyr::select(unique_id, species, run, treatment, temperature_c, growth_rate, nitrate, phosphate, salinity)

glimpse(wg_old)

#subset by species
wg_old_u <- wg_old %>% 
  filter(species == "Ul")

wg_old_h <- wg_old %>% 
  filter(species == "Hm")

 #Run brt
#Old data Ulva
histo_old_u <- wg_old_u %>%
  ggplot(aes(growth_rate)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_old_u 

#transform for better normality
wg_old_u <- wg_old_u %>%
  mutate(
    sqrt_growth = sqrt(growth_rate + 57 + 0.01)  # shift + small constant
  )

histo_u_sqrt <- wg_old_u %>%
  ggplot(aes(sqrt_growth)) +
  geom_histogram(binwidth=1, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_u_sqrt 

#run model
set.seed(42)

brt_old_ulva <- gbm.step(
  data = wg_old_u,
  gbm.x = c("salinity", "nitrate", "phosphate", "temperature_c"),  # predictors
  gbm.y = "sqrt_growth", # response
  family = "gaussian",
  tree.complexity = 2,
  learning.rate = 0.001,
  bag.fraction = 0.5,
  n.folds = 5
)

#Ask for the contributions of each predictor
brt_old_ulva$contributions

#Plotting these contributions. Adapt the number of plots to the number of predictors.
gbm.plot(brt_old_ulva, n.plots=4, write.title = F,
         common.scale = TRUE, y.label="Fitted function",
         show.contrib = TRUE)




#Old data Hypnea
histo_old_h <- wg_old_h %>%
  ggplot(aes(growth_rate)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_old_h

#transform for better normality
wg_old_h <- wg_old_h %>%
  mutate(
    sqrt_growth = sqrt(growth_rate + 88 + 0.01)  # shift + small constant
  )

histo_h_sqrt <- wg_old_h %>%
  ggplot(aes(sqrt_growth)) +
  geom_histogram(binwidth=1, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_h_sqrt

set.seed(42)

brt_old_hypnea <- gbm.step(
  data = wg_old_h,
  gbm.x = c("salinity", "nitrate", "phosphate", "temperature_c"),  # predictors
  gbm.y = "sqrt_growth", # response
  family = "gaussian",
  tree.complexity = 2,
  learning.rate = 0.001,
  bag.fraction = 0.5,
  n.folds = 5
)

#Ask for the contributions of each predictor
brt_old_hypnea$contributions

#Plotting these contributions. Adapt the number of plots to the number of predictors.
gbm.plot(brt_old_hypnea, n.plots=4, write.title = F,
         common.scale = TRUE, y.label="Fitted function",
         show.contrib = TRUE)





#New Data **
#**there was no difference in salinity or temperature in this dataset, only N matters, so no need to use brt as N = 100%
waste_growth_2025 <- read.csv("../data/wastewater/input/growth_2025/p2_growth_all.csv")
waste_growth_2025$treatment <- as.factor(waste_growth_2025$treatment)
waste_growth_2025$species <- as.factor(waste_growth_2025$species)
waste_growth_2025$plant_id <- as.factor(waste_growth_2025$plant_id)
waste_growth_2025$salinity <- as.factor(waste_growth_2025$salinity)
waste_growth_2025$run <- as.factor(waste_growth_2025$run)
waste_growth_2025$nitrate <- round(as.numeric(waste_growth_2025$N_g_m3_2d), digits = 1)
waste_growth_2025$phosphate <- round(as.numeric(waste_growth_2025$P_g_m3_2d), digits = 1)


waste_growth_2025 <- waste_growth_2025 %>%
  mutate(date = mdy(date))

glimpse(waste_growth_2025)

wg_u <- waste_growth_2025 %>% 
  filter(species == "u")

wg_h <- waste_growth_2025 %>% 
  filter(species == "h")



# ULVA
histo_u <- wg_u %>%
  ggplot(aes(growth_d9)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_u

wg_u <- wg_u %>%
  mutate(
    sqrt_growth = sqrt(growth_d9 + 24.11 + 0.01)  # shift + small constant
  )

histo_u_sqrt <- wg_u %>%
  ggplot(aes(sqrt_growth)) +
  geom_histogram(binwidth=.5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
histo_u_sqrt

set.seed(42)

brt_ulva <- gbm.step(
  data = wg_u,
  gbm.x = c("nitrate", "phosphate"),  # predictors
  gbm.y = "sqrt_growth", # response
  family = "gaussian",
  tree.complexity = 2,
  learning.rate = 0.001,
  bag.fraction = 0.5,
  n.folds = 5
)

# Fit a BRT model for Ulva (wg_u)
set.seed(42)
#Find out the best iteration
brt_waste_ulva <- gbm.step(
  data = wg_u,
  gbm.x = c("nitrate", "phosphate"),  # predictors
  gbm.y = "growth_d9", # response
  family = "gaussian",
  tree.complexity = 3,
  learning.rate = 0.001,
  bag.fraction = 0.5,
  n.folds = 5
)


#Ask for the contributions of each predictor
brt_waste_ulva$contributions

#Plotting these contributions. Adapt the number of plots to the number of predictors.
gbm.plot(brt_waste_ulva, n.plots=2, write.title = F,
         common.scale = TRUE, y.label="Fitted function",
         show.contrib = TRUE)


#Hypnea

#Ask for the contributions of each predictor
brt_waste_ulva$contributions

#Plotting these contributions. Adapt the number of plots to the number of predictors.
gbm.plot(brt_waste_ulva, n.plots=4, write.title = F,
         common.scale = TRUE, y.label="Fitted function",
         show.contrib = TRUE)






#ALL BELOW FOR MORE DETAILED PLOTTING WITH HORIZONTAL LINE AT 0 Per Jade's output

# Extract partial dependence data
pd_sal <- plot.gbm(brt_ulva, i.var = "salinity", return.grid = TRUE)
pd_N <- plot.gbm(brt_ulva, i.var = "nitrate", return.grid = TRUE)
pd_P <- plot.gbm(brt_ulva, i.var = "phosphate", return.grid = TRUE)

# Center y-values around 0
pd_sal$y <- pd_sal$y - mean(pd_sal$y)
pd_N$y <- pd_N$y - mean(pd_N$y)
pd_P$y <- pd_P$y - mean(pd_P$y)

sal_contrib <- brt_ulva$contributions %>%
  filter(var == "salinity") %>%
  pull(rel.inf)

N_contrib <- brt_ulva$contributions %>%
  filter(var == "nitrate") %>%
  pull(rel.inf)

P_contrib <- brt_ulva$contributions %>%
  filter(var == "phosphate") %>%
  pull(rel.inf)

# Then use:
labs(title = paste0("Salinity (", round(sal_contrib, 1), "%)"))
labs(title = paste0("Nitrate (", round(N_contrib, 1), "%)"))
labs(title = paste0("Phosphate (", round(P_contrib, 1), "%)"))

# PDP: Salinity
p_sal <- ggplot(pd_sal, aes(x = salinity, y = y, group = 1)) +
  geom_step(color = "black", size = 1, direction = "vh") +  # ← step instead of line
  geom_point(size = 2, color = "black") +  # optional: show points
  geom_rug(data = wg_u, aes(x = salinity), inherit.aes = FALSE, sides = "b", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkorchid3") +
  ylim(min(c(pd_sal$y, pd_treat$y)), max(c(pd_sal$y, pd_treat$y))) +
  labs(title = paste0("Salinity (", round(sal_contrib, 1), "%)"),
       x = "Salinity", y = "Partial effect on Growth")
p_sal

# PDP: Treatment
p_nit <- ggplot(pd_N, aes(x = nitrate, y = y, group = 1)) +
  geom_step(color = "black", size = 1, direction = "hv") +  # ← step instead of line
  geom_hline(yintercept = 0, linetype = "dashed", color = "deeppink3") +
  geom_point(color = "black", size = 2) +
  ylim(min(c(pd_nit$y, pd_sal$y)), max(c(pd_nit$y, pd_sal$y))) +
  theme_bw() +
  labs(title = paste0("Nitrate (", round(nit_contrib, 1), "%)"),
       x = "Nitrate", y = "Partial effect on Growth")
p_nit

