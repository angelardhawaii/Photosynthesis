# Step 5 (bypassing populate scripts)
# This code will combine the Ek and irradiance datasets after wrangling them
# The irradiance comes from three separate datasets 2021, 2023, 2025
# The ek etc. comes from the combined dataset with wastewater data
# Tidyverse workflow; cleaner way to do the same thing as previous scripts
# P1 = days 1–3 → threshold = day1 Ek
 #P2 = days 4–6 → threshold = day5 Ek (or interpolated if missing)
 #P3 = days 7–9 → threshold = day9 Ek
 #Entire run starts at day1 rlc_time
 #Entire run ends at day9 rlc_end_time
# By Angela Richards Dona
# Created February 20, 2026
# ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(hms)
library(janitor)

# -----------------------------
# 1) Load EK (combined dataset)
# -----------------------------
ek <- read_csv(
  "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ww_combo_ek_alpha_norm.csv"
) %>%
  clean_names() %>%
  mutate(
    date = ymd(date),
    rlc_end_time = hms::as_hms(rlc_end_time),
    rlc_day = as.integer(rlc_day),
    nitrate = factor(as.character(nitrate)),
    salinity = factor(salinity),
    sal_bin = factor(sal_bin),
    plant_id = factor(plant_id),
    pmax_min = pmax * 60,
    delta_npq = as.numeric(delta_npq),
    year = year(date),
    run_combo = factor(run_combo)
  ) %>%
  group_by(date, rlc_day) %>%
  arrange(rlc_end_time, .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    rlc_order1 = factor(floor((rank - 1) / 3) + 1)
  ) %>%
  ungroup()

# Check for duplicates before we start joining and summarizing by individual keys.
ek %>%
  count(unique_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))
# -----------------------------
# 2) Load irradiance files
# -----------------------------
irrad <- readRDS("/Users/angela/src/Photosynthesis/data/wastewater/transformed/irrad_1min.rds") %>%
  mutate(
    date = ymd(date), # should already be Date, but harmless
    epar = as.numeric(epar)
  )
# Check irradiance data
irrad %>%
  summarise(min_date = min(date), max_date = max(date), n = n(),
            min_epar = min(epar, na.rm = TRUE), max_epar = max(epar, na.rm = TRUE))

# ------------------------------------------------------------
# 3) Build thresholds table per individual (day 1 base rows)
#    Join day5/day9 thresholds by individual keys
# ------------------------------------------------------------
id_key <- c("id", "run_combo", "species") # keys to join day1, day5, day9 rows by individual

day1 <- ek %>%     
  filter(rlc_day == 1) %>%
  select(all_of(id_key), plant_id,
         nitrate, salinity, sal_bin, temp, rlc_day, treat_letter, rlc_order1,
         date, rlc_time, f0, fm, fv_fm, alpha_est, ek_est, pmax, pmax_min) %>%
  rename(
    start_date = date,
    start_time = rlc_time, # this is the start time of the entire run (day1 rlc_time)
    ek_day1    = ek_est,
    pmax_min_day1 = pmax_min
  )

day5 <- ek %>%
  filter(rlc_day == 5) %>%
  select(all_of(id_key), ek_est, pmax_min) %>%
  rename(
    ek_day5    = ek_est,
    pmax_min_day5 = pmax_min
  )

day9 <- ek %>%
  filter(rlc_day == 9) %>%
  select(all_of(id_key), rlc_end_time, ek_est, pmax_min) %>%
  rename(
    end_time   = rlc_end_time, # this is the end time of the entire run (day9 rlc_end_time)
    ek_day9    = ek_est,
    pmax_min_day9 = pmax_min
  )

hsat <- day1 %>%
  left_join(day5, by = id_key) %>%
  left_join(day9, by = id_key) %>%
  mutate(
    scheme = if_else(is.na(ek_day5), "no_day5", "with_day5"),
    ek_day5 = if_else(is.na(ek_day5), (ek_day1 + ek_day9) / 2, ek_day5),
    pmax_min_day5 = if_else(is.na(pmax_min_day5),
                            (pmax_min_day1 + pmax_min_day9) / 2,
                            pmax_min_day5)
  )

calc_hsat <- function(irrad_df, start_date, days_to_consider, threshold) {
  if (is.na(threshold)) {
    return(tibble(hsat_min = NA_real_, daylight_min = NA_real_, relhsat = NA_real_))
  }
  
  dates <- start_date + lubridate::days(days_to_consider - 1)
  
  irrad_df %>%
    filter(date %in% dates) %>%
    summarise(
      hsat_min = sum(epar >= threshold, na.rm = TRUE),
      daylight_min = sum(epar >= 1, na.rm = TRUE),
      relhsat = if_else(daylight_min > 0, hsat_min / daylight_min, NA_real_),
      .groups = "drop"
    )
}

hsat_periods <- hsat %>%
  rowwise() %>%
  mutate(
    p1 = list(calc_hsat(irrad, start_date,
                        if (scheme == "no_day5") 1:4 else 1:3,
                        ek_day1)),
    p2 = list(calc_hsat(irrad, start_date,
                        if (scheme == "no_day5") 5:9 else 4:6,
                        if (scheme == "no_day5") ek_day9 else ek_day5)),
    p3 = list(calc_hsat(irrad, start_date,
                        7:9,
                        if (scheme == "with_day5") ek_day9 else NA_real_))
  ) %>%
  tidyr::unnest_wider(p1, names_sep = "_") %>%  # -> p1_hsat_min, p1_daylight_min, p1_relhsat
  tidyr::unnest_wider(p2, names_sep = "_") %>%
  tidyr::unnest_wider(p3, names_sep = "_") %>%
  ungroup()

hsat_periods <- hsat_periods %>%
  mutate(
    p3_hsat_min     = coalesce(p3_hsat_min, 0),
    p3_daylight_min = coalesce(p3_daylight_min, 0),
    
    hsat_total     = p1_hsat_min + p2_hsat_min + p3_hsat_min,
    daylight_total = p1_daylight_min + p2_daylight_min + p3_daylight_min,
    relhsat_total  = if_else(daylight_total > 0, hsat_total / daylight_total, NA_real_)
  )

# ------------------------------------------------------------
# Calculate DSPI
# ------------------------------------------------------------

hsat_periods <- hsat_periods %>%
  mutate(
    # DSPI per period
    p1_dspi = p1_hsat_min * pmax_min_day1,
    
    p2_dspi = if_else(
      scheme == "no_day5",
      p2_hsat_min * pmax_min_day9,
      p2_hsat_min * pmax_min_day5
    ),
    
    p3_dspi = if_else(
      scheme == "with_day5",
      p3_hsat_min * pmax_min_day9,
      0
    ),
    
    # total DSPI
    dspi_total = p1_dspi + p2_dspi + p3_dspi
  )
# Convert to mol/m-2
hsat_periods <- hsat_periods %>%
  mutate(
    p1_dspi_mol = round(p1_dspi / 1e6, 2),
    p2_dspi_mol = round(p2_dspi / 1e6, 2),
    p3_dspi_mol = round(p3_dspi / 1e6, 2),
    dspi_total_mol = round(dspi_total / 1e6, 2)
  )
# Save data
write_csv(hsat_periods, "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ww_combo_hsat_dspi.csv")
