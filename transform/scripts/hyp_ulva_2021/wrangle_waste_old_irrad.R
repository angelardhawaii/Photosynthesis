#NEW code for irradiance used in 2026 for use with combined dataset with wastewater
# By Angela Richards Dona
# Created February 20, 2026
# This code will prepare the two older irradiance datasets for use in hsat script and will save them as new files
# ------------------------------------------------------------
# wrangle_irradiance.R
# Build a clean irradiance dataset at 1-minute resolution.
#
# Inputs (CSV):
#   - irrad2021.csv  (date, time, epar, lanai_side)
#   - irrad2023.csv  (date, time, epar)
#   - waste_irrad.csv (date, time, epar)
#
# Output:
#   - ./data_irradiance/processed/irrad_1min.rds
#   - ./data_irradiance/processed/irrad_1min.csv
# ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(hms)
library(janitor)

TZ <- "Pacific/Honolulu"

# ---- file paths ----
paths <- tibble(
  source = c("irrad2021", "irrad2023", "waste_irrad"),
  path   = c("/Users/angela/src/Photosynthesis/data/hyp_ulv_old/2021_2022_irrad.csv",
             "/Users/angela/src/Photosynthesis/data/hyp_ulv_old/2023_irrad.csv",
             "/Users/angela/src/Photosynthesis/data/wastewater/input/waste_irrad.csv")
)

irrad_raw <- paths %>%
  mutate(
    data = map(path, ~ read_csv(
      .x,
      show_col_types = FALSE,
      col_types = cols(
        date = col_character(),
        time = col_character(),
        epar = col_double(),
        .default = col_guess()
      )
    ) %>%
      clean_names())
  ) %>%
  unnest(data) %>%
  select(-any_of("chk")) %>%
  mutate(
    lanai_side = if ("lanai_side" %in% names(.)) as.character(lanai_side) else NA_character_,
    
    # parse date correctly by source (source comes from `paths`)
    date = lubridate::parse_date_time(date, orders = 
                                        c("ymd", "mdy", "Y-m-d", "m/d/Y", "Y/m/d")) %>% as.Date(),
    
    time = stringr::str_trim(time),
    time = stringr::str_replace(time, "\\.\\d+$", ""),
    time = hms::parse_hms(time),
    
    dt   = force_tz(as.POSIXct(date) + time, TZ),
    epar = as.numeric(epar)
  ) %>%
  filter(!is.na(dt), !is.na(epar))

# ---- collapse to 1 value per minute (union of all sources/sides) ----
# choose how to handle overlaps within the same minute:
#   mean = neutral
#   median = robust
#   max = conservative for "minutes above Ek"
AGG <- "max"  

reducer <- switch(
  AGG,
  mean   = function(x) mean(x, na.rm = TRUE),
  median = function(x) median(x, na.rm = TRUE),
  max    = function(x) max(x, na.rm = TRUE)
)

irrad_1min <- irrad_raw %>%
  mutate(
    dt_min = floor_date(dt, "minute"),
    date = as.Date(dt_min, tz = TZ)
  ) %>%
  group_by(dt_min, date) %>%
  summarise(
    epar = reducer(epar),
    n_raw = n(),
    # provenance (optional; helps debug coverage)
    lanai_side = paste(sort(unique(na.omit(lanai_side))), collapse = "+"),
    source = paste(sort(unique(source)), collapse = "+"),
    .groups = "drop"
  ) %>%
  arrange(dt_min)

# ---- quick sanity prints ----
print(irrad_1min %>% summarise(min_dt = min(dt_min), max_dt = max(dt_min), n_minutes = n()))
print(irrad_1min %>% summarise(overlap_minutes = sum(n_raw > 1)))

# ---- save ----
saveRDS(irrad_1min, "/Users/angela/src/Photosynthesis/data/wastewater/transformed/irrad_1min.rds")
write_csv(irrad_1min, "/Users/angela/src/Photosynthesis/data/wastewater/transformed/irrad_1min.csv")



#____________________________
#OLD CODE SAVED JUST IN CASE
# These are time series of the irradiance measurements
# The 3rd column is the value of interest (irradiance)
ewa_irradiance_files <- dir("./data_irradiance/ewa", full.names = TRUE)
ewa_irradiance <- do.call(rbind, lapply(ewa_irradiance_files, read.csv))
ewa_irradiance$Lanai.Side <- as.factor("ewa")

dia_irradiance_files <- dir("./data_irradiance/dia", full.names = TRUE)
dia_irradiance <- do.call(rbind, lapply(dia_irradiance_files, read.csv))
dia_irradiance$Lanai.Side <- as.factor("diamond")

irradiance <- rbind(ewa_irradiance, dia_irradiance)
rm(dia_irradiance, ewa_irradiance, dia_irradiance_files, ewa_irradiance_files)

# Loads all the files in the same data frame

# Make sure the date is loaded as date
irradiance$posix_date <- as.POSIXct(irradiance$Date, format = "%m/%d/%y", tz = "")
irradiance$date_time <- as.POSIXct(paste(irradiance$Date, irradiance$Time, sep = " "), format = "%m/%d/%y %H:%M:%S")