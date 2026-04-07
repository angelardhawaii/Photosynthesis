# Step 4
#This script takes the cleaned and combined data from PAM output and runs the Phytotools package to
#produce a final dataset that can be used for analysis.
# Step 3 is from script waste_photosyn_combine_sgd.R

#To be used with data that already has all unique IDs
#add column for treatment in next step from this output
# last modified in november 2025
# By Angela Richards Dona

#open appropriate libraries
library("phytotools")
library("hash")
library("dplyr")
library(lubridate)
library(hms)
library(stringr)
library(tidyr)
library(purrr)
library(tidyverse)

ww_df <- read_csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/ps_combo_allyears.csv")
glimpse(ww_df)

# Clean types & core derived columns
ww_df <- ww_df %>%
  mutate(
    # Make sure these are numeric; coerce safely
    delta_npq = as.numeric(delta_npq),
    fm        = as.numeric(fm),
    f         = as.numeric(f),
    epar      = as.numeric(epar),
    yield     = as.numeric(yield),
    retr      = as.numeric(retr),
    npq       = as.numeric(npq),
    
    # Effective quantum yield (numeric, not factor)
    quan_yield = round((fm - f) / fm, 3),
    
    # Optional: Y(II) at Epar == 0 (dark step); stays NA otherwise
    quan_yield_at0 = if_else(epar == 0, yield, NA_real_)
  )

# 3) Keep only positive rETR
selected_df <- ww_df %>%
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
    rlc_end_time = hms::as_hms(max(time, na.rm = TRUE)),
    .groups = "drop"
  )

# 5) Fit Webb per light curve (vectorised by groups via list-columns)
#    Assumes fitWebb(epar, yield, normalize=TRUE) returns a list with numeric vectors $alpha and $ek (length 4 each)
# Helper function to extract elements from named vectors safely
safe_extract <- function(x, name) {
  if (!is.null(x) && name %in% names(x)) {
    return(x[[name]])
  } else {
    return(NA_real_)
  }
}

webb_params <- selected_df %>%
  group_by(unique_id) %>%
  group_modify(~ {
    fit <- tryCatch(
      fitWebb(.x$epar, .x$yield, normalize = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      return(tibble(
        alpha_est = NA_real_,
        alpha_se = NA_real_,
        alpha_t = NA_real_,
        alpha_pvalue = NA_real_,
        ek_est = NA_real_,
        ek_se = NA_real_,
        ek_t = NA_real_,
        ek_pvalue = NA_real_
      ))
    }
    
    tibble(
      alpha_est    = safe_extract(fit$alpha, "Estimate"),
      alpha_se     = safe_extract(fit$alpha, "Std. Error"),
      alpha_t      = safe_extract(fit$alpha, "t value"),
      alpha_pvalue = safe_extract(fit$alpha, "Pr(>|t|)"),
      ek_est       = safe_extract(fit$ek, "Estimate"),
      ek_se        = safe_extract(fit$ek, "Std. Error"),
      ek_t         = safe_extract(fit$ek, "t value"),
      ek_pvalue    = safe_extract(fit$ek, "Pr(>|t|)")
    )
  }) %>%
  ungroup()

per_curve_with_webb <- per_curve %>%
  left_join(webb_params, by = "unique_id") %>%
  mutate(
    # Pmax from the first alpha/ek values
    pmax = round(alpha_est * ek_est, 2)
  )

# Populate delta_npq column using first positive value of npq as min
# Require a minimum of 3 values for delta_npq calculation
ww_df <- ww_df %>%
  group_by(unique_id) %>%
  mutate(
    # Clean NPQ: ignore negatives and NAs when we compute min/max
    npq_clean = if_else(npq > 0, as.numeric(npq), NA_real_),
    
    # Count how many valid NPQ values we have
    n_valid_npq = sum(!is.na(npq_clean)),
    
    # Find the smallest *positive* epar value where NPQ is valid
    min_epar_pos = suppressWarnings(
      min(if_else(epar > 0 & !is.na(npq_clean), epar, NA_real_), na.rm = TRUE)
    ),
    
    # NPQ at that first valid positive epar
    npq_min = {
      idx <- which(epar == min_epar_pos & !is.na(npq_clean))
      val <- if (length(idx) == 0) NA_real_ else npq_clean[idx[1]]
      val
    },
    
    # Maximum valid NPQ over the curve
    npq_max = if_else(
      n_valid_npq > 0,
      max(npq_clean, na.rm = TRUE),
      NA_real_
    ),
    
    row_num = row_number(),
    
    # Only compute delta_npq if we have enough valid NPQ points
    delta_npq = if_else(
      row_num == 1 & n_valid_npq >= 3 & !is.na(npq_min) & !is.na(npq_max),
      npq_max - npq_min,
      # First row but not enough data → NA; other rows → 0 (as before)
      if_else(row_num == 1, NA_real_, 0)
    )
  ) %>%
  select(-npq_clean, -n_valid_npq, -min_epar_pos, -npq_min, -npq_max, -row_num) %>%
  ungroup()

# 6) Join summaries back to rows (row-level dataset with curve-level columns filled)
selected_df_aug <- selected_df %>%
  left_join(per_curve_with_webb, by = "unique_id")

first_row_of_rlc <- ww_df %>%
  filter(epar == 0 & is.na(npq))
last_row_rlc <- subset(ww_df, delta_npq > 0)

# Join the per-curve summaries into first-row-only subset
first_row_joined <- first_row_of_rlc %>%
  left_join(per_curve_with_webb, by = "unique_id")  # Adjust join key if needed

# Now safely build your final df
ww_pmax_ek_alpha <- data.frame(
  id               = first_row_joined$id,
  unique_id        = first_row_joined$unique_id,
  plant_id         = first_row_joined$plant_id,
  date             = as_date(first_row_joined$date),
  species          = first_row_joined$species,
  rlc_time         = as_hms(first_row_joined$time),
  treatment        = first_row_joined$treatment,
  temp             = first_row_joined$temp,
  nitrate          = first_row_joined$nitrate,
  salinity         = first_row_joined$salinity,
  sal_bin          = first_row_joined$sal_bin,
  n_g_m3_2d        = first_row_joined$n_g_m3_2d,
  p_g_m3_2d        = first_row_joined$p_g_m3_2d,
  treat_letter     = first_row_joined$treat_letter,
  run_combo        = first_row_joined$run_combo,
  rlc_end_time     = as_hms(first_row_joined$rlc_end_time),
  rlc_day          = first_row_joined$rlc_day,
  f0               = first_row_joined$f0,
  fm               = first_row_joined$fm,
  f                = first_row_joined$f,
  fm_prime         = first_row_joined$fm_prime,
  fv_fm            = first_row_joined$fv_fm,
  delta_npq        = first_row_joined$delta_npq,
  pmax             = first_row_joined$pmax,
  retr_max         = first_row_joined$retr_max,
  retr_max_ypoint1 = first_row_joined$retr_max_ypoint1,
  npq_max          = first_row_joined$npq_max,
  alpha_est        = round(first_row_joined$alpha_est, 3),
  alpha_se         = round(first_row_joined$alpha_se, 3),
  alpha_t          = round(first_row_joined$alpha_t, 2),
  alpha_pvalue     = round(first_row_joined$alpha_pvalue, 3),
  ek_est           = round(first_row_joined$ek_est, 1),
  ek_se            = round(first_row_joined$ek_se, 1),
  ek_t             = round(first_row_joined$ek_t, 2),
  ek_pvalue        = round(first_row_joined$ek_pvalue, 3),
  lunar_phase      = first_row_joined$lunar_phase,
  illumination     = round(first_row_joined$illumination, 3)
)


#________________________________
# save to file
write.csv(ww_pmax_ek_alpha, "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ww_combo_ek_alpha_norm.csv")
