# Photosynthesis script for combining 2023 SGD dataset with 2025 wastewater dataset
# Both were cleaned by pam_data_cleaner_tidy.R function THUS the _tidy in the names
#Important to remember this
# Script may need to add some variables
# By Angela Richards Dona
# 8/20/25

# load libraries for plots and tables
library(ggplot2)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(sjPlot)
library(sjmisc)
# for tidying and joining and wrangling data in general
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
library(dplyr)
library(readr)

#open cleaned ps datasets and make columns for the various meta data
ps_2025 <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/pam_waste_p2_clean_tidy.csv")
#ps_2023 <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/pam_2023_clean_tidy.csv")
#ps_2021 <- read.csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/transformed/pam_2021data_clean_tidy.csv")
ps_2021 <- read.csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/transformed/hyp_ulva_all_runs_clean.csv")

glimpse(ps_2025)
glimpse(ps_2021)

#Scrubby scrubby all datasets
ps_2025 <- clean_names(ps_2025) #clean up ugly names
ps_2021 <- clean_names(ps_2021)

#2021 DATA_________________________________________________________________
ps_2021 <- ps_2021 %>%
  mutate(id = tolower(id)) #make all ids lowercase to match other datasets

ps_2021 <- ps_2021 %>%
  rename(retr = "r_etr", yield = "y_ii", fm_prime = "fm_2") 

# add unique_id column)
ps_2021 <- ps_2021 %>%
  mutate(unique_id = paste(date, tolower(id), sep = "_"))  

ps_2021 <- ps_2021 %>%
  mutate(date = mdy(date))

ps_2021 <- ps_2021 %>%
  mutate(time = as_hms(time))

# Fill in F0 from f
# since f0 is character, coerce once
ps_2021 <- ps_2021 %>%
  mutate(
    f      = as.numeric(f),
    epar   = as.numeric(epar),
    f0     = suppressWarnings(as.numeric(f)) * NA_real_ + suppressWarnings(as.numeric(f0)) # keep numeric type
  ) %>%
  group_by(unique_id) %>%
  mutate(
    f0 = coalesce(
      f0,
      as.numeric(round(first(f[epar == 0], na_rm = TRUE)))
    )
  ) %>%
  ungroup()

# Fill in Fm from fm_prime
ps_2021 <- ps_2021 %>%
  mutate(
    fm_prime      = as.numeric(fm_prime),
    epar   = as.numeric(epar),
    fm     = suppressWarnings(as.numeric(fm_prime)) * NA_real_ + suppressWarnings(as.numeric(fm)) # keep numeric type
  ) %>%
  group_by(unique_id) %>%
  mutate(
    fm = coalesce(
      fm,
      as.numeric(round(first(fm_prime[epar == 0], na_rm = TRUE)))
    )
  ) %>%
  ungroup()

# calculate fv_fm
ps_2021 <- ps_2021 %>%
  mutate(fv_fm = round((fm - f0) / fm, 3))

# column for run
ps_2021 <- ps_2021 %>%
  mutate(run = case_when(
    date %in% c("2021-09-20", "2021-09-24", "2021-09-28") ~ 1,
    date %in% c("2021-09-21", "2021-09-25", "2021-09-29") ~ 1,
    date %in% c("2021-09-30", "2021-10-04", "2021-10-08") ~ 2,
    date %in% c("2021-10-01", "2021-10-05", "2021-10-09") ~ 2,
    date %in% c("2021-10-11", "2021-10-15", "2021-10-19") ~ 3,
    date %in% c("2021-10-21", "2021-10-25", "2021-10-29") ~ 4,
    date %in% c("2021-11-01", "2021-11-05", "2021-11-09") ~ 5,
    date %in% c("2021-11-04", "2021-11-08", "2021-11-12") ~ 5,
    date %in% c("2022-02-11", "2022-02-15", "2022-02-19") ~ 6,
    date %in% c("2022-02-21", "2022-02-25", "2022-03-01") ~ 7,
    date %in% c("2022-04-21", "2022-04-25", "2022-04-29") ~ 8,
    date %in% c("2022-10-04", "2022-10-08", "2022-10-12") ~ 9,
    date %in% c("2022-10-13", "2022-10-17", "2022-10-21") ~ 10,
    date %in% c("2022-10-14", "2022-10-18", "2022-10-22") ~ 10,
    TRUE ~ NA_real_   # fallback for all other dates
  ))

#add species column
ps_2021 <- ps_2021 %>%
  mutate(species = str_sub(id, 1, 1))

ps_2021 <- ps_2021 %>%
  mutate(treatment = as.factor(str_sub(id2, 5, 7)))

#add column for lettered treatments
ps_2021 <- ps_2021 %>%
  mutate(treat_letter = case_when(
    treatment == 0   ~ "a",
    treatment == 1   ~ "b",
    treatment == 2   ~ "c",
    treatment == 2.5   ~ "d",
    treatment == 3   ~ "e",
    treatment == 4   ~ "f"
  ))

# add column for nitrate to dataset
ps_2021 <- ps_2021 %>%
  mutate(nitrate = case_when(
    treatment == 0   ~ "0.5",
    treatment == 1   ~ "14",
    treatment == 2   ~ "27",
    treatment == 2.5   ~ "53",
    treatment == 3   ~ "53",
    treatment == 4   ~ "80"
  ))

#add column for N
ps_2021 <- ps_2021 %>%
  mutate(n_g_m3_2d = case_when(
    treatment == 0   ~ "0.06",
    treatment == 1   ~ "1.7",
    treatment == 2   ~ "3.3",
    treatment == 2.5   ~ "6.4",
    treatment == 3   ~ "6.4",
    treatment == 4   ~ "9.7"
  ))

#add column for P
ps_2021 <- ps_2021 %>%
  mutate(p_g_m3_2d = case_when(
    treatment == 0   ~ "0.001",
    treatment == 1   ~ "0.001",
    treatment == 2   ~ "0.037",
    treatment == 2.5   ~ "0.41",
    treatment == 3   ~ "0.41",
    treatment == 4   ~ "0.94"
  ))

#add column for salinity
ps_2021 <- ps_2021 %>%
  mutate(salinity = case_when(
    treatment == 0   ~ "35",
    treatment == 1   ~ "35",
    treatment == 2   ~ "28",
    treatment == 2.5   ~ "28",
    treatment == 3   ~ "18",
    treatment == 4   ~ "11"
  ))

# 2025 DATA_________________________________________________________________
ps_2025 <- ps_2025 %>%
  mutate(date = as.Date(date))
ps_2025 <- ps_2025 %>%
  mutate(time = as_hms(time))
 
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

ps_2025 <- ps_2025 %>%
  mutate(unique_id = paste(date, tolower(id), sep = "_")) 

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
# add column for plant_id based on key
ps_2025_plantid_key <- read_csv("/Users/angela/src/Photosynthesis/data/wastewater/input/waste_id_plant_key.csv")

ps_2025 <- ps_2025 %>%
  left_join(ps_2025_plantid_key, by = "id")



#COMBINE DATASETS-----------------------------------------------
# Show all unique column names across the three datasets
all_names <- unique(c(names(ps_2025), names(ps_2021)))
all_names

# See which columns are missing from each dataset
setdiff(all_names, names(ps_2021))  # columns missing from ps_2021
setdiff(all_names, names(ps_2025))  # columns missing from ps_2023


# Columns common to both datasets
common_cols <- Reduce(intersect, list(names(ps_2025), names(ps_2021)))
common_cols

glimpse(ps_2025)
glimpse(ps_2021)

# force character for any type mismatched columns
to_factor <- c("treatment", "nitrate", "n_g_m3_2d", "p_g_m3_2d", "salinity", "run", "plant_id", "treat_letter")
ps_2021 <- ps_2021 %>% mutate(across(all_of(to_factor), as.factor))
ps_2025 <- ps_2025 %>% mutate(across(all_of(to_factor), as.factor))
ps_2025 <- ps_2025 %>% mutate(delta_npq = as.character(delta_npq))
# Select columns in the same order, then stack
ps_combo <- bind_rows(ps_2021, ps_2025) 


# Modify combined dataset

#make date a POSIX date and add a year column
ps_combo <- ps_combo %>%
  mutate(date = ymd(date))


#Assign new run values to avoid issues of repeats and to get true chronology
ps_combo <- ps_combo %>%
  mutate(run_combo = as.factor(case_when(
    date %in% c("2021-09-20", "2021-09-24", "2021-09-28") ~ 1,
    date %in% c("2021-09-21", "2021-09-25", "2021-09-29") ~ 1,
    date %in% c("2021-09-30", "2021-10-04", "2021-10-08") ~ 2,
    date %in% c("2021-10-01", "2021-10-05", "2021-10-09") ~ 2,
    date %in% c("2021-10-11", "2021-10-15", "2021-10-19") ~ 3,
    date %in% c("2021-10-21", "2021-10-25", "2021-10-29") ~ 4,
    date %in% c("2021-11-01", "2021-11-05", "2021-11-09") ~ 5,
    date %in% c("2021-11-04", "2021-11-08", "2021-11-12") ~ 5,
    date %in% c("2022-02-11", "2022-02-15", "2022-02-19") ~ 6,
    date %in% c("2022-02-21", "2022-02-25", "2022-03-01") ~ 7,
    date %in% c("2022-04-21", "2022-04-25", "2022-04-29") ~ 8,
    date %in% c("2022-10-04", "2022-10-08", "2022-10-12") ~ 9,
    date %in% c("2022-10-13", "2022-10-17", "2022-10-21") ~ 10,
    date %in% c("2022-10-14", "2022-10-18", "2022-10-22") ~ 10,
    date >= as.Date("2025-01-23") & date <= as.Date("2025-01-31") ~ 11,
    date >= as.Date("2025-02-17") & date <= as.Date("2025-02-25") ~ 12,
    date >= as.Date("2025-03-06") & date <= as.Date("2025-03-18") ~ 13,
    date >= as.Date("2025-03-17") & date <= as.Date("2025-03-28") ~ 14,
    TRUE ~ NA_real_  # Optional: keep NA for unmatched dates
  )))

ps_combo <- ps_combo %>%
  mutate(rlc_day = case_when(
    date == as.Date("2025-03-10") & time < hms::as_hms("12:00:00") ~ 5,
    date == as.Date("2025-03-10") & time >= hms::as_hms("12:00:00") ~ 1,
    date == as.Date("2025-03-14") & time < hms::as_hms("12:00:00") ~ 9,
    date == as.Date("2025-03-14") & time >= hms::as_hms("12:00:00") ~ 5,
    date %in% as.Date(c("2021-09-20", "2021-09-21", "2013-09-30", "2021-10-01", "2021-10-11", 
                        "2021-10-21", "2021-11-01", "2021-11-04", "2022-02-11", "2022-02-21",
                        "2022-04-21", "2022-10-04", "2022-10-13", "2022-10-14", "2025-01-23", 
                        "2025-02-17", "2025-03-06", "2025-03-17", "2025-03-20")) ~ 1,
    date %in% as.Date(c("2021-09-24", "2021-09-25", "2021-10-04", "2021-10-05", "2013-10-15", 
                        "2021-10-25", "2021-11-05", "2021-11-08", "2022-02-15", "2022-02-25",
                        "2022-04-25", "2022-10-08", "2022-10-17", "2022-10-18", "2025-01-27",
                        "2025-02-21", "2025-03-21", "2025-03-24")) ~ 5,
    date %in% as.Date(c("2021-09-28", "2021-09-29", "2021-10-08", "2021-10-09", "2021-10-19",
                        "2021-10-29", "2021-11-09", "2021-11-12", "2022-02-19", "2022-03-01",
                        "2022-04-29", "2022-10-12", "2022-10-21", "2022-10-22", "2025-01-31", 
                        "2025-02-25", "2025-03-18", "2025-03-25", "2025-03-28")) ~ 9
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
write.csv(ps_combo, file = "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ps_combo_2021_2025.csv")



#_______________________________________________________________________________
# 2023 DATA_____________________DON'T USE_______________________________________
#_______________________________________________________________________________

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

# Also not in use but worked so hard on this stuff__________________________________
# data merge to add id, fv_fm, and temp columns from the key dataset "ids_all.csv"
# Rename to your actual data frame names:
ps_2021_key  <- read_csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/input/ids_all.csv")   # has date, time, id, fv_fm, temp

# Ensure types & select the join keys/columns we want from the key table
key_clean <- ps_2021_key %>%
  transmute(
    date = as.Date(as.character(date)),   # force factor → character → Date
    time = as_hms(as.character(time)),    # ensure "HH:MM:SS" format
    id    = tolower(as.character(id)),    # ensure character and lowercase
    fv_fm = as.numeric(fv_fm),
    temp  = as.numeric(temp)
  )
# 1) Normalize types in the main table and define block boundaries
ps_2021_prepped <- ps_2021 %>%
  mutate(
    date = as.Date(as.character(date)),   # ensure Date
    time = as_hms(as.character(time))     # ensure HMS
  ) %>%
  arrange(date, time) %>%    # make sure rows are in time order within date
  group_by(date) %>%
  mutate(block = cumsum(epar == 0)) %>%    # each start row (epar==0) begins a new block
  ungroup()


# 2) Join only at the start rows (but easiest is to join all rows on date+time —
#    only the start rows will actually match because only those times exist in key_clean)
ps_2021_joined <- ps_2021_prepped %>%
  left_join(key_clean, by = c("date", "time")) %>%
  group_by(date, block) %>%
  tidyr::fill(id.y, fv_fm.y, temp, .direction = "down") %>%
  ungroup() 

# 3) Carry the id/fv_fm/temp DOWN within each block so all 9 rows get them
ps_2021 <- ps_2021_joined %>%
  group_by(date, block) %>%                      # block from your earlier mutate(cumsum(epar==0))
  tidyr::fill(id.y, fv_fm.y, temp, .direction = "downup") %>%  # fill across all rows in the block
  ungroup() %>%
  select(-block)

#Merge id.x and id.y columns into single id column
ps_2021 <- ps_2021 %>%
  mutate(
    # normalize blanks to NA and trim
    id.x = na_if(trimws(as.character(id.x)), ""),
    id.y = na_if(trimws(as.character(id.y)), ""),
    # merge without losing info
    id = case_when(
      !is.na(id.x) & !is.na(id.y) & id.x != id.y ~ paste(id.x, id.y, sep = "|"),
      !is.na(id.x)                                ~ id.x,
      !is.na(id.y)                                ~ id.y,
      TRUE                                        ~ NA_character_
    )
  ) %>%
  select(-id.x, -id.y)                    # drop the old columns

#Merge fv_fm.x and fv_fm.y columns into single fv_fm column
ps_2021 <- ps_2021 %>%
  mutate(
    fv_fm.x = as.numeric(fv_fm.x),
    fv_fm.y = as.numeric(fv_fm.y),
    fv_fm = case_when(
      !is.na(fv_fm.x) & !is.na(fv_fm.y) & fv_fm.x != fv_fm.y ~ (fv_fm.x + fv_fm.y) / 2,  # average if different
      !is.na(fv_fm.x)                                         ~ fv_fm.x,
      !is.na(fv_fm.y)                                         ~ fv_fm.y,
      TRUE                                                     ~ NA_real_
    )
  ) %>%
  select(-fv_fm.x, -fv_fm.y)            # drop the old columns

#

ps_2021 <- ps_2021 %>%
  mutate(unique_id = paste(date, tolower(id), sep = "_"))  

#More joins with other datasets to get rlc_order, plant_id, treatment, and some temp values
growth_2021 <- read_csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/input/all_runs_growth_011723.csv")

growth_2021 <- clean_names(growth_2021)

growth_2021 <- growth_2021 %>%
  mutate(date = mdy(date))

growth_2021 <- growth_2021 %>%
  mutate(unique_id = paste(date, tolower(id), sep = "_")) %>%
  mutate(temp = str_sub(temperature, 1, 2))

# column for run
growth_2021 <- growth_2021 %>%
  mutate(run = case_when(
    date == "2021-09-28" ~ 1,
    date == "2021-09-29" ~ 1,
    date == "2021-10-08" ~ 2,
    date == "2021-10-09" ~ 2,
    date == "2021-10-19" ~ 3,
    date == "2021-10-29" ~ 4,
    date == "2021-11-09" ~ 5,
    date == "2021-11-12" ~ 5,
    date == "2022-02-19" ~ 6,
    date == "2022-03-01" ~ 7,
    date == "2022-04-29" ~ 8,
    date == "2022-10-12" ~ 9,
    date == "2022-10-21" ~ 10,
    date == "2022-10-22" ~ 10,
    TRUE ~ NA_real_   # fallback for all other dates
  ))

# Count unique_ids in each dataset
n_ps  <- n_distinct(ps_2021$unique_id)
n_gr  <- n_distinct(growth_2021$unique_id)

# Overlap
n_overlap <- ps_2021 %>%
  distinct(unique_id) %>%
  semi_join(growth_2021 %>% distinct(unique_id),
            by = "unique_id") %>%
  nrow()

cat("Unique IDs in ps_2021:   ", n_ps, "\n")
cat("Unique IDs in growth_2021:", n_gr, "\n")
cat("Overlap:                  ", n_overlap, "\n")

# IDs in growth_2021 but NOT in ps_2021
missing_in_ps <- growth_2021 %>%
  distinct(unique_id) %>%
  anti_join(ps_2021 %>% distinct(unique_id),
            by = "unique_id")

missing_in_ps

#join growth_2021 to ps_2021
# 1) Create a key from growth_2021 with (run, id) → (plant_id, temp, treatment, rlc_order)
growth_key <- growth_2021 %>%
  select(run, unique_id, plant_id, temp, treatment, rlc_order) %>%
  distinct()

# (Optional) sanity checks on the key
# - exactly one row per (run, id)?
dups_in_key <- growth_key %>%
  count(run, unique_id) %>% filter(n > 1)
if (nrow(dups_in_key) > 0) {
  warning("growth_key has duplicate (run, unique_id). Resolve before joining.")
}

# 2) Left-join into ps_2021 by (run, unique_id)
ps_2021_aug <- ps_2021 %>%
  left_join(growth_key, by = c("run", "unique_id"))

# 3) Quick QA
stopifnot(nrow(ps_2021_aug) == nrow(ps_2021))  # row count unchanged

# Which (run, id) in ps_2021 didn’t find a match in growth_2021?
missing_pairs <- ps_2021_aug %>%
  distinct(run, unique_id) %>%
  anti_join(growth_key %>% distinct(run, unique_id), by = c("run","unique_id"))

# How many rows have missing joins for each joined column?
colSums(is.na(ps_2021_aug[, c("plant_id", "temp.y", "treatment", "rlc_order")]))

## Optional sanity checks
stopifnot(nrow(ps_2021_aug) == nrow(ps_2021))  # row count unchanged

# Any (run,id) in ps_2021 that didn’t find a match in growth_2021?
missing_in_key <- ps_2021 %>%
  anti_join(growth_key, by = c("run","unique_id")) %>%
  distinct(run, unique_id)
nrow(missing_in_key)

# Are there duplicate (run,id) in the key (should be 1)?
dups_in_key <- growth_key %>%
  count(run, unique_id) %>%
  filter(n > 1)
nrow(dups_in_key)

#Merge temp.x and temp.y columns into single temp column add unique_id column)
ps_2021 <- ps_2021_aug %>%
  mutate(
    temp.x = as.numeric(temp.x),
    temp.y = as.numeric(temp.y),
    temp = case_when(
      !is.na(temp.x) & !is.na(temp.y) & temp.x != temp.y ~ (temp.x + temp.y) / 2,  # average if different
      !is.na(temp.x)                                       ~ temp.x,
      !is.na(temp.y)                                       ~ temp.y,
      TRUE                                                 ~ NA_real_
    )
  ) %>%
  select(-temp.x, -temp.y)            # drop the old columns

colSums(is.na(ps_2021[, c("plant_id", "temp", "treatment", "rlc_order")]))
