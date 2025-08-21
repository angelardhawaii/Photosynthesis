# Growth analysis script modified for use with wastewater dataset in May 2025
# Modifying both to be able to combine
# all_growth_waste is 2025 dataset, all_growth is 2021-2022 dataset, all_2023 is 2023 dataset
# Workflow may not be perfectly in order
# By Angela Richards Dona
# 5/22/25

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
library(piecewiseSEM)
library(lubridate)
library(hms)
library(janitor)
library(lunar)
library(suncalc)
library(scico)

#open weight dataset and make columns for growth rate from initial and final weights
#all_growth <- read.csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/input/all_runs_growth_011723.csv")
all_growth_waste <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/input/growth_2025/p2_growth_all.csv")
all_growth <- read.csv("/Users/angela/src/Photosynthesis/data/hyp_ulv_old/all_runs_growth_011723.csv")
all_2023 <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/input/growth_2025/hyp_ulv_2023_growth.csv")

#Scrubby scrubby all datasets
all_growth <- clean_names(all_growth) #clean up ugly names in 2021-2022 dataset
all_2023 <- clean_names(all_2023)
all_growth_waste <- clean_names(all_growth_waste)

# 2021-2022 DATA
#add a column to the all_growth dataset for growth_d9
all_growth <- all_growth %>%
  mutate(growth_d9 = round(((final_weight - initial_weight)/initial_weight) * 100, 4)
         )

#remove the C in the values for temperature in old dataset
all_growth <- all_growth %>%
  mutate(temp = str_sub(temperature, 1,2)) %>%
  select(-temperature)

#Work on oldest dataset
all_growth$temp <- as.numeric(all_growth$temp)
all_growth$treatment <- as.character(all_growth$treatment)

#remove .0s from treatment column in old dataset (integers, using modulus does NOT work)
all_growth <- all_growth %>%
  mutate(treat_letter = case_when(
    treatment == 0   ~ "a",
    treatment == 1   ~ "b",
    treatment == 2   ~ "c",
    treatment == 2.5 ~ "d",
    treatment == 3   ~ "e",
    treatment == 4   ~ "f",
    TRUE ~ treatment  # keep all other values as-is
  ))


#add three columns for salinity, Nitrate (n_g_m3_2d) and Phosphate (p_g_m3_2d)
all_growth <- all_growth %>%
  mutate(n_g_m3_2d = case_when(
    treatment == 0 ~ 0.0607,
    treatment == 1 ~ 1.7363,
    treatment == 2 ~ 3.2905,
    treatment == 2.5 ~ 6.4232,
    treatment == 3 ~ 6.4232,
    treatment == 4 ~ 9.7137
  ))

all_growth <- all_growth %>%
  mutate(p_g_m3_2d = case_when(
    treatment == 0 ~ 0.0012,
    treatment == 1 ~ 0.0012,
    treatment == 2 ~ 0.0373,
    treatment == 2.5 ~ 0.4081,
    treatment == 3 ~ 0.4081,
    treatment == 4 ~ 0.943
  ))

all_growth <- all_growth %>%
  mutate(salinity = case_when(
    treatment == 0 ~ 35,
    treatment == 1 ~ 35,
    treatment == 2 ~ 28,
    treatment == 2.5 ~ 28,
    treatment == 3 ~ 18,
    treatment == 4 ~ 11
  ))

all_growth <- all_growth %>%
  mutate(nitrate = case_when(
    treatment == 0  ~ 0.5,
    treatment == 1  ~ 14,
    treatment == 2  ~ 27,
    treatment == 2.5  ~ 53,
    treatment == 3  ~ 53,
    treatment == 4  ~ 80,
    TRUE ~ NA_real_ #fallback if no match
  ))

#Change names of weight variables
all_growth <- all_growth %>%
  rename(i_weight = initial_weight, f_weight = final_weight)

#add prefix to plant_id and run values to distinguish between the datasets since overlap occurs
all_growth <- all_growth %>%
  mutate(plant_id = paste0(plant_id, "_sgd")) %>%
  mutate(run = paste0(run, "_sgd"))




# 2023 DATA_________________________________________________________________
#add species column
all_2023 <- all_2023 %>%
  mutate(species = str_sub(id, 1, 1))

#add temperature column for all_2023 dataset, the first run was all hyp then all ulv
#so data has to be treated in blocks of 8 (per bin) run 1 and 4 per bin runs 2, 3 
#first add id_num column
all_2023 <- all_2023 %>%
  mutate(id_num = str_sub(id, 3,4))

all_2023$id_num <- as.integer(all_2023$id_num) 
all_2023$run <- as.integer(all_2023$run) 

# Define the repeating temp pattern
temps <- c(25, 27, 30)

pick_from <- function(block, pattern) pattern[((block - 1) %% length(pattern)) + 1]

#For run 1, groups are blocks of 8 IDs; for runs 2 & 3, blocks of 4 IDs
all_2023 <- all_2023 %>%
  mutate(
    block8 = ceiling(id_num / 16),   # used for run 1
    block4 = ceiling(id_num / 8),   # used for runs 2 & 3
    temp = case_when(
      run == 1 ~ pick_from(block8, temps),
      run %in% c(2, 3) ~ pick_from(block4, temps),
      TRUE ~ NA_real_
    )
  ) %>%
  select(-block8, -block4)

#add column for salinity
all_2023 <- all_2023 %>%
  mutate(salinity = case_when(
    treatment == 1 ~ 35,
    treatment == 2 ~ 27,
    treatment == 3 ~ 21,
    treatment == 4 ~ 15
  ))

#add column for N
all_2023 <- all_2023 %>%
  mutate(n_g_m3_2d = case_when(
    treatment == 1 ~ 0.0607,
    treatment == 2 ~ 10.4422,
    treatment == 3 ~ 6.4353,
    treatment == 4 ~ 4.4926
  ))
#add column for P
all_2023 <- all_2023 %>%
  mutate(p_g_m3_2d = case_when(
    treatment == 1 ~ 0.0012,
    treatment == 2 ~ 2.0652,
    treatment == 3 ~ 1.7169,
    treatment == 4 ~ 0.423
  ))

# Add column for nitrate in micromols
all_2023 <- all_2023 %>%
  mutate(nitrate = case_when(
    treatment == 1 ~ 0.5,
    treatment == 2 ~ 86,
    treatment == 3 ~ 53,
    treatment == 4 ~ 37
  ))

#Add column for lettered treatments - Order these from lowest to highest, 1, 3, 4, 2
all_2023 <- all_2023 %>%
  mutate(treat_letter = case_when(
    treatment == 1   ~ "g",
    treatment == 2   ~ "j",
    treatment == 3   ~ "i",
    treatment == 4   ~ "h"
  ))

#add prefix to plant_id and run values to distinguish between the datasets since overlap occurs
all_2023 <- all_2023 %>%
  mutate(plant_id = paste0(plant_id, "_23")) %>%
  mutate(run = paste0(run, "_23"))




# 2025 DATA___________________________________________________________________
#add column for lettered treatments
all_growth_waste <- all_growth_waste %>%
  mutate(treat_letter = case_when(
    treatment == 1   ~ "k",
    treatment == 3   ~ "l",
    treatment == 5   ~ "m",
    treatment == 7   ~ "n"
  ))

# add column for nitrate to wastewater dataset
all_growth_waste <- all_growth_waste %>%
  mutate(nitrate = case_when(
    treatment == 1  ~ 80,
    treatment == 3  ~ 245,
    treatment == 5  ~ 748,
    treatment == 7  ~ 2287,
    TRUE ~ NA_real_ #fallback if no match
  ))
#add prefix to plant_id and run values to distinguish between the datasets since overlap occurs
all_growth_waste <- all_growth_waste %>%
  mutate(plant_id = paste0(plant_id, "_waste")) %>%
  mutate(run = paste0(run, "_waste"))

#JOIN ALL DATASETS
# Show all unique column names across the three datasets
all_names <- unique(c(names(all_growth), names(all_2023), names(all_growth_waste)))
all_names

# See which columns are missing from each dataset
setdiff(all_names, names(all_growth))  # columns missing from all_growth
setdiff(all_names, names(all_2023))  # columns missing from all_2023
setdiff(all_names, names(all_growth_waste))  # columns missing from all_growth_waste




#COMBINE DATASETS-----------------------------------------------
# 1) Columns common to all three
common_cols <- Reduce(intersect, list(names(all_growth), names(all_2023), names(all_growth_waste)))
common_cols

# force character for any type mismatched columns
to_char <- c("treatment", "treat_letter")
all_growth <- all_growth %>% mutate(across(all_of(to_char), as.character))
all_2023 <- all_2023 %>% mutate(across(all_of(to_char), as.character))
all_growth_waste <- all_growth_waste %>% mutate(across(all_of(to_char), as.character))

# 2) Select those columns in the same order, then stack
combined_growth <- bind_rows(
  all_growth %>% select(all_of(common_cols)),
  all_2023 %>% select(all_of(common_cols)),
  all_growth_waste %>% select(all_of(common_cols))
)

# Modify combined dataset
#Add columns for species and change to "h" and "u" for Hypnea and Ulva
combined_growth <- combined_growth %>%
  mutate(
    species = case_when(
      species %in% c("Hm", "h") ~ "h",
      species %in% c("Ul", "u") ~ "u",
      TRUE ~ species  # Keep other values unchanged
    )
  )

#make date a POSIX date and add a year column
combined_growth <- combined_growth %>%
  mutate(date = mdy(date)) %>%
  mutate(year = year(date))

combined_growth <- combined_growth %>%
  mutate(unique_id = paste(date, id, sep = "_")) 


#Assign new run values to avoid issues of repeats and to get true chronology
combined_growth <- combined_growth %>%
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
    date == as.Date("2022-10-12") ~ 9,
    date == as.Date("2022-10-21") ~ 10,
    date == as.Date("2022-10-22") ~ 10,
    date == as.Date("2023-03-17") ~ 11,
    date == as.Date("2023-03-18") ~ 11,
    date == as.Date("2023-03-28") ~ 12,
    date == as.Date("2023-03-31") ~ 12,
    date == as.Date("2023-04-11") ~ 13,
    date == as.Date("2023-04-14") ~ 13,
    date == as.Date("2025-01-31") ~ 14,
    date == as.Date("2025-02-25") ~ 15,
    date == as.Date("2025-03-14") ~ 16,
    date == as.Date("2025-03-18") ~ 16,
    date == as.Date("2025-03-25") ~ 17,
    date == as.Date("2025-03-28") ~ 17
  ))


#add lunar phase to the combined dataset
library(lunar)
library(suncalc)
combined_growth <- combined_growth %>%
  mutate(lunar_phase = lunar.phase(date, name = TRUE))

#for more precision on lunar illumination
moon_data <- getMoonIllumination(combined_growth$date)

combined_growth <- combined_growth %>%
  mutate(illumination = moon_data$fraction,
         lunar_illumination = case_when(
           illumination < 0.1 ~ "new",
           illumination >= 0.1 & illumination < 0.4 ~ "crescent",
           illumination >= 0.4 & illumination < 0.6 ~ "first_q / last_q",
           illumination >= 0.6 & illumination < 0.9 ~ "gibbous",
           illumination >= 0.9 ~ "full"
         ))

#Save this new combined dataset
write.csv(combined_growth, file = "/Users/angela/src/Photosynthesis/data/wastewater/transformed/combined_data_all.csv")


