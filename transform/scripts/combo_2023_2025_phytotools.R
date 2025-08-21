#This script takes the cleaned data from PAM output and runs the Phytotools package to
#produce a final dataset that can be used for analysis.
#get ek and alpha etc

#To be used with data that already has all unique IDs
#add column for treatment in next step from this output

#open appropriate libraries
library(phytotools)
library(hash)
library(dplyr)
library(lubridate)
library(hms)
library(stringr)
library(tidyr)
library(purrr)
library(tidyverse)

# 1) Read
ps_combo_pt <- readr::read_csv(
  "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ps_combo_2023_2025.csv",
  show_col_types = FALSE
)

# 2) Clean types & core derived columns
ps_combo_pt <- ps_combo_pt %>%
  mutate(
    # Make sure these are numeric; coerce safely
    delta_npq = as.numeric(delta_npq),
    fm        = as.numeric(fm),
    f         = as.numeric(f),
    epar      = as.numeric(epar),
    yield     = as.numeric(yield),
    retr      = as.numeric(retr),
    npq       = as.numeric(npq),
    # If `time` is a time string, keep it as character or parse if you want ordering
    # time    = hms::as_hms(time),
    
    # Effective quantum yield (numeric, not factor)
    quan_yield = round((fm - f) / fm, 3),
    
    # Optional: Y(II) at Epar == 0 (dark step); stays NA otherwise
    quan_yield_at0 = if_else(epar == 0, yield, NA_real_)
  )

# 3) Keep only positive rETR
selected_df <- ps_combo_pt %>%
  filter(retr > 0)

# 4) Per-curve (per unique_id) summaries
per_curve <- selected_df %>%
  group_by(unique_id) %>%
  summarise(
    retr_max = max(retr, na.rm = TRUE),
    retr_max_ypoint1 = {
      vals <- retr[yield > 0.1]
      if (length(vals) == 0) NA_real_ else max(vals, na.rm = TRUE)
    },
    npq_max = {
      vals <- npq[yield > 0.1]
      if (length(vals) == 0) NA_real_ else max(vals, na.rm = TRUE)
    },
    rlc_end_time = max(time, na.rm = TRUE),
    .groups = "drop"
  )

# 5) Fit Webb per light curve (vectorised by groups via list-columns)
#    Assumes fitWebb(epar, yield, normalize=TRUE) returns a list with numeric vectors $alpha and $ek (length 4 each)
webb_params <- selected_df %>%
  group_by(unique_id) %>%
  reframe(
    fit = list(fitWebb(epar, yield, normalize = TRUE))
  ) %>%
  mutate(
    alpha = list(fit[[1]]$alpha),
    ek    = list(fit[[1]]$ek)
  ) %>%
  select(-fit) %>%
  tidyr::unnest_wider(alpha, names_sep = "_") %>%  # alpha_1 ... alpha_4
  tidyr::unnest_wider(ek,    names_sep = "_")      # ek_1 ... ek_4

per_curve_with_webb <- per_curve %>%
  left_join(webb_params, by = "unique_id") %>%
  mutate(
    # Pmax from the first alpha/ek values
    pmax = round(alpha_1 * ek_1, 2)
  )
# 6) Join summaries back to rows (row-level dataset with curve-level columns filled)
selected_df_aug <- selected_df %>%
  left_join(per_curve, by = "unique_id")

# 7) (Optional) a per-curve table with both summaries and Webb params
per_curve_with_webb <- per_curve %>%
  left_join(webb_params, by = "unique_id")


#rlc_day_assign <- read.csv("/Users/Angela/src/work/limu/phytotools_alpha_ek/data_input/date_day_assignment.csv")

#rlc_days_by_date = array(dim = length(dates))
#hash_of_rlc_days_by_date <- hash(rlc_day_assign$Date, rlc_day_assign$RLC.Day)
#for (i in 1: length(dates)) {
#        rlc_days_by_date[i] = hash_of_rlc_days_by_date[[dates[i]]]
#}

#lanai_side_by_date = array(dim = length(dates))
#hash_of_lanai_side_by_date <- hash(rlc_day_assign$Date, tolower(rlc_day_assign$Lanai.side))
#for (i in 1:length(dates)) {
#        lanai_side_by_date[i] = as.character(hash_of_lanai_side_by_date[[dates[i]]])
#}

first_row_of_rlc <- subset(ps_combo_pt, Epar == 0 & NPQ == "-")
last_row_rlc <- subset(ps_combo_pt, deltaNPQ > 0)

# build the final data frame
final_df <- data.frame("Date" = first_row_of_rlc$posix_date, 
                       "rlc_end_time" = rlc_end_times,
                       "Specimen ID" = first_row_of_rlc$ID,
                       #"Plant ID" = first_row_of_rlc$plant.ID,
                       "deltaNPQ" = last_row_rlc$deltaNPQ,
                       "Species" = first_row_of_rlc$species,
                       "Treatment" = first_row_of_rlc$treatment,
                       "Temp (°C)" = first_row_of_rlc$temp,
                       "RLC Time" = first_row_of_rlc$Time,
                       "RLC Day" = first_row_of_rlc$rlc_day,
                       "Run" = first_row_of_rlc$run,
                       "pmax" = pmax,
                       "rETRmax" = rETRMaxes,
                       "rETRmaxYpoint1" = rETRmaxYpoint1,
                       "NPQmax" = NPQmax,
                       "alpha" = round(alpha, digits = 3),
                       "ek" = round(ek, digits = 1)
)


# save to file
write.csv(final_df, "data_output/wastewater/hyp_ulv_wastewater_ek_alpha_normalized.csv")
write.csv (final_df, "../irradiance_ek/data_ek/wastewater/hyp_ulv_wastewater_ek_alpha_normalized.csv")
write.csv(final_df, "../algal_growth_photosynthesis/data_input/wastewater/hyp_ulv_wastewater_ek_alpha_normalized.csv")

