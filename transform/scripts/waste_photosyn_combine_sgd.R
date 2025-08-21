# Photosynthesis script for combining 2023 SGD dataset with 2025 wastewater dataset
# Both will be cleaned by pam_data_cleaner_tidy.R function
# Will add some of the higher N values from 2021 dataset
# By Angela Richards Dona
# 8/20/25



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

#open cleaned ps datasets and make columns for the various meta data
ps_2025 <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/pam_waste_p2_clean_tidy.csv")
ps_2023 <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/pam_2023_clean_tidy.csv")
#ps_2021 <- read.csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/transformed/hyp_ulva_all_runs_clean2.csv")

glimpse(ps_2023)
glimpse(ps_2025)

#Scrubby scrubby all datasets
ps_2025 <- clean_names(ps_2025) #clean up ugly names
ps_2023 <- clean_names(ps_2023)



# 2023 DATA_________________________________________________________________
#add species column
ps_2023 <- ps_2023 %>%
  mutate(species = str_sub(id, 1, 1))

# column for run
ps_2023 <- ps_2023 %>%
  mutate(run = case_when(
    date %in% c("2023-03-09", "2023-03-13", "2023-03-17") ~ 1,
    date %in% c("2023-03-10", "2023-03-14", "2023-03-18") ~ 1,
    date %in% c("2023-03-20", "2023-03-24", "2023-03-28") ~ 2,
    date %in% c("2023-03-23", "2023-03-27", "2023-03-31") ~ 2,
    date %in% c("2023-04-03", "2023-04-07", "2023-04-11") ~ 3,
    date %in% c("2023-04-06", "2023-04-10", "2023-04-14") ~ 3,
    TRUE ~ NA_real_   # fallback for all other dates
    ))
#add temperature column for ps_2023 dataset, the first run was all hyp then all ulv
#so data has to be treated in blocks of 8 (per bin) run 1 and 4 per bin runs 2, 3 
#first add id_num column
ps_2023 <- ps_2023 %>%
  mutate(id_num = as.integer(str_sub(id, -2, -1)))

ps_2023$id_num <- as.integer(ps_2023$id_num) 
ps_2023$run <- as.integer(ps_2023$run) 

# Define the repeating temp and treatment patterns
temps <- c(25, 27, 30)

pick_from <- function(block, pattern) pattern[((block - 1) %% length(pattern)) + 1]

#For run 1, groups are blocks of 16 IDs; for runs 2 & 3, blocks of 8 IDs
ps_2023 <- ps_2023 %>%
  mutate(
    block16 = ceiling(id_num / 16),   # used for run 1
    block8 = ceiling(id_num / 8),   # used for runs 2 & 3
    temp = case_when(
      run == 1 ~ pick_from(block16, temps),
      run %in% c(2, 3) ~ pick_from(block8, temps),
      TRUE ~ NA_real_
    )) %>% 
  select(-block16, -block8)

# Add a column for treatment from id_num
  ps_2023 <- ps_2023 %>%
    mutate(
      id_num = as.integer(id_num),
      run    = as.integer(run),
      treatment = case_when(
        # Run 1: pairs of IDs (1–2→1, 3–4→2, 5–6→3, 7–8→4, then repeat)
        run == 1        ~ ((ceiling(id_num / 2) - 1) %% 4) + 1L,
        
        # Runs 2 & 3: every 4th ID by position (1,5,9,…→1; 2,6,10,…→2; etc.)
        run %in% c(2,3) ~ ((id_num - 1) %% 4) + 1L,
        
        TRUE ~ NA_integer_
      )
    )

#add column for salinity
ps_2023 <- ps_2023 %>%
  mutate(salinity = case_when(
    treatment == 1 ~ 35,
    treatment == 2 ~ 27,
    treatment == 3 ~ 21,
    treatment == 4 ~ 15
  ))

#add column for N
ps_2023 <- ps_2023 %>%
  mutate(n_g_m3_2d = case_when(
    treatment == 1 ~ 0.0607,
    treatment == 2 ~ 10.4422,
    treatment == 3 ~ 6.4353,
    treatment == 4 ~ 4.4926
  ))
#add column for P
ps_2023 <- ps_2023 %>%
  mutate(p_g_m3_2d = case_when(
    treatment == 1 ~ 0.0012,
    treatment == 2 ~ 2.0652,
    treatment == 3 ~ 1.7169,
    treatment == 4 ~ 0.423
  ))

# Add column for nitrate in micromols
ps_2023 <- ps_2023 %>%
  mutate(nitrate = case_when(
    treatment == 1 ~ 0.5,
    treatment == 2 ~ 86,
    treatment == 3 ~ 53,
    treatment == 4 ~ 37
  ))

#Add column for lettered treatments - Order these from lowest to highest, 1, 3, 4, 2
ps_2023 <- ps_2023 %>%
  mutate(treat_letter = case_when(
    treatment == 1   ~ "g",
    treatment == 2   ~ "j",
    treatment == 3   ~ "i",
    treatment == 4   ~ "h"
  ))


# 2025 DATA_________________________________________________________________
ps_2025 <- ps_2025 %>%
  mutate(species = str_sub(id, 1, 1))

# column for run
ps_2025 <- ps_2025 %>%
  mutate(run = case_when(
    date %in% c("2025-01-23", "2025-01-27", "2025-01-31") ~ 1,
    date %in% c("2025-02-17", "2025-02-21", "2025-02-25") ~ 2,
    date %in% c("2025-03-06", "2025-03-10", "2025-03-14", "2025-03-18") ~ 3,
    date %in% c("2025-03-17", "2025-03-21", "2025-03-25") ~ 4,
    date %in% c("2025-03-20", "2025-03-24", "2025-03-28") ~ 4,
    TRUE ~ NA_real_   # fallback for all other dates
  ))

# add id_num column
ps_2025 <- ps_2025 %>%
  mutate(id_num = as.integer(str_sub(id, -2, -1)))

ps_2025$run <- as.integer(ps_2025$run) 

# Add a column for temp
ps_2025 <- ps_2025 %>%
  mutate(temp = 22) # All runs in 2025 were at 22°C

# Add a column for treatment from id_num
ps_2025 <- ps_2025 %>%
  mutate(treatment =
     # every 4th ID by position (1,5,9,…→1; 2,6,10,…→2; etc.)
      ((id_num - 1) %% 4) + 1L
    )

#add column for salinity
ps_2025 <- ps_2025 %>%
  mutate(
    salinity = if_else(((ceiling(id_num / 4) - 1) %% 2) == 0, 18, 22)
  )
#add column for lettered treatments
ps_2025 <- ps_2025 %>%
  mutate(treat_letter = case_when(
    treatment == 1   ~ "k",
    treatment == 2   ~ "l",
    treatment == 3   ~ "m",
    treatment == 4   ~ "n"
  ))

# add column for nitrate to wastewater dataset
ps_2025 <- ps_2025 %>%
  mutate(nitrate = case_when(
    treatment == 1  ~ 80,
    treatment == 2  ~ 245,
    treatment == 3  ~ 748,
    treatment == 4  ~ 2287,
    TRUE ~ NA_real_ #fallback if no match
  ))

#add column for N
ps_2025 <- ps_2025 %>%
  mutate(n_g_m3_2d = case_when(
    treatment == 1 ~ 9.714,
    treatment == 2 ~ 29.748,
    treatment == 3 ~ 90.823,
    treatment == 4 ~ 277.690
  ))
#add column for P
ps_2025 <- ps_2025 %>%
  mutate(p_g_m3_2d = case_when(
    treatment == 1 ~ 0.995,
    treatment == 2 ~ 3.981,
    treatment == 3 ~ 15.676,
    treatment == 4 ~ 62.455
  ))

#COMBINE DATASETS-----------------------------------------------
# Show all unique column names across the three datasets
all_names <- unique(c(names(ps_2025), names(ps_2023)))
all_names

# See which columns are missing from each dataset
setdiff(all_names, names(ps_2023))  # columns missing from all_growth
setdiff(all_names, names(ps_2025))  # columns missing from all_2023


# Columns common to both datasets
common_cols <- Reduce(intersect, list(names(ps_2023), names(ps_2023)))
common_cols

# force character for any type mismatched columns
to_char <- c("treatment", "treat_letter")
ps_2023 <- ps_2023 %>% mutate(across(all_of(to_char), as.character))
ps_2025 <- ps_2025 %>% mutate(across(all_of(to_char), as.character))

# Select columns in the same order, then stack
ps_combo <- bind_rows(
  ps_2023 %>% select(all_of(common_cols)),
  ps_2025 %>% select(all_of(common_cols))
)

# Modify combined dataset

#make date a POSIX date and add a year column
ps_combo <- ps_combo %>%
  mutate(date = ymd(date))

ps_combo <- ps_combo %>%
  mutate(unique_id = paste(date, id, sep = "_")) 


#Assign new run values to avoid issues of repeats and to get true chronology
ps_combo <- ps_combo %>%
  mutate(run_combo = case_when(
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
ps_combo <- ps_combo %>%
  mutate(lunar_phase = lunar.phase(date, name = TRUE))

#for more precision on lunar illumination
moon_data <- getMoonIllumination(ps_combo$date)

ps_combo <- ps_combo %>%
  mutate(illumination = moon_data$fraction,
         lunar_illumination = case_when(
           illumination < 0.1 ~ "new",
           illumination >= 0.1 & illumination < 0.4 ~ "crescent",
           illumination >= 0.4 & illumination < 0.6 ~ "first_q / last_q",
           illumination >= 0.6 & illumination < 0.9 ~ "gibbous",
           illumination >= 0.9 ~ "full"
         ))

#Save this new combined dataset
write.csv(ps_combo, file = "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ps_combo_2023_2025.csv")
