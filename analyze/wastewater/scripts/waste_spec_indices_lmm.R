# Spectral Reflectance Indices and LMM Analysis
# By Angela Richards Dona
# Wastewater SGD 2025 experiment — Ulva lactuca and Hypnea musciformis
#
# Spectral data: day 9, 5 of 6 run_combos (22-26), ~350-750 nm
# Joins to ww_combo_hsat_dspi via `id`
# Tests whether pigment-sensitive indices respond to nitrate/salinity
# where Hsat and DSPI did not.

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(stringr)
library(ggplot2)
library(signal)
library(lme4)
library(lmerTest)
library(MuMIn)

# ── 1. Load and smooth spectral data ─────────────────────────────────────────

spec_files <- list.files(
  path       = "/Users/angela/src/Photosynthesis/data/wastewater/input/spec/",
  pattern    = "\\.csv$",
  full.names = TRUE
)

spec_raw <- spec_files %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE, name_repair = "minimal"))

# Identify wavelength columns (numeric decimal names)
wl_cols <- grep("^\\d+\\.\\d+$", names(spec_raw), value = TRUE)
wl_num  <- as.numeric(wl_cols)

# Ensure reflectance columns are numeric
spec_raw <- spec_raw %>%
  mutate(across(all_of(wl_cols), as.numeric))

# Savitzky-Golay smoothing (window = 11, polynomial = 2)
smooth_sg <- function(x, w = 11, p = 2) signal::sgolayfilt(x, p = p, n = w)

refl_mat  <- as.matrix(spec_raw[, wl_cols])
sm_mat    <- t(apply(refl_mat, 1, smooth_sg))
colnames(sm_mat) <- wl_cols

spec_smoothed <- bind_cols(
  tibble(id_full = spec_raw$id),
  as_tibble(sm_mat)
)

# ── 2. Compute spectral indices from reflectance ──────────────────────────────
# Target wavelengths (use nearest available):
#   Blue peak absorption : ~450 nm
#   Green reflectance peak: ~550 nm
#   Red absorption trough : ~670 nm
#   Red edge base         : ~700 nm
#   Red edge top / NIR    : ~740 nm  (max available ~750 nm)

nearest_wl <- function(target) {
  wl_cols[which.min(abs(wl_num - target))]
}

wl_blue  <- nearest_wl(450)
wl_green <- nearest_wl(550)
wl_red   <- nearest_wl(670)
wl_re700 <- nearest_wl(700)
wl_nir   <- nearest_wl(740)

message("Wavelength assignments:")
message("  Blue  : ", wl_blue,  " nm")
message("  Green : ", wl_green, " nm")
message("  Red   : ", wl_red,   " nm")
message("  RE700 : ", wl_re700, " nm")
message("  NIR   : ", wl_nir,   " nm")

spec_idx <- spec_smoothed %>%
  rename(
    R_blue  = all_of(wl_blue),
    R_green = all_of(wl_green),
    R_red   = all_of(wl_red),
    R_re700 = all_of(wl_re700),
    R_nir   = all_of(wl_nir)
  ) %>%
  select(id_full, R_blue, R_green, R_red, R_re700, R_nir) %>%
  mutate(
    # NDVI: (NIR - Red) / (NIR + Red)  — overall greenness / chlorophyll proxy
    ndvi     = (R_nir - R_red)   / (R_nir + R_red),
    # Simple Ratio: NIR / Red  — linear chlorophyll index
    sr       = R_nir / R_red,
    # Chlorophyll Index Green: (NIR / Green) - 1  — sensitive to Chl a+b
    ci_green = (R_nir / R_green) - 1,
    # NDRE: (NIR - RedEdge) / (NIR + RedEdge)  — red edge chlorophyll index
    ndre     = (R_nir - R_re700) / (R_nir + R_re700),
    # Red:Blue ratio  — shifts with accessory pigment ratios
    rb_ratio = R_red / R_blue
  )

# ── 3. Average replicates, parse ID to match ww_combo format ─────────────────
# id_full = "h125_1", "h125_2", "h125_3"  → id = "h125"

spec_mean_idx <- spec_idx %>%
  mutate(id = str_sub(id_full, 1, 4)) %>%
  group_by(id) %>%
  summarise(
    across(c(ndvi, sr, ci_green, ndre, rb_ratio), mean, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )

# ── 4. Load ww_combo and join ─────────────────────────────────────────────────

ww <- read_csv(
  "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ww_combo_hsat_dspi.csv",
  show_col_types = FALSE
) %>%
  filter(temp == 22) %>%  # 2025 data only
  mutate(
    id         = tolower(trimws(id)),
    nitrate    = as.factor(nitrate),
    salinity   = as.factor(salinity),
    run_combo  = as.factor(run_combo),
    plant_id   = as.factor(plant_id),
    rlc_order1 = as.factor(rlc_order1)
  )

# Join: spec was measured on day 9; one row per id in both tables
ww_spec <- ww %>%
  inner_join(spec_mean_idx, by = "id")

message("\nRows after join: ", nrow(ww_spec),
        "  (expected ~", nrow(spec_mean_idx), " per rlc_day)")

# Check coverage
message("run_combos in joined data: ",
        paste(sort(unique(ww_spec$run_combo)), collapse = ", "))

# ── 5. Species subsets ────────────────────────────────────────────────────────

ulva_spec   <- filter(ww_spec, species == "u")
hypnea_spec <- filter(ww_spec, species == "h")

# ── 6. Visual check: mean spectra by nitrate treatment ───────────────────────

# Pivot smoothed spectra back to long for plotting (subset to 400-750 nm)
spec_plot_long <- spec_smoothed %>%
  mutate(id = str_sub(id_full, 1, 4)) %>%
  group_by(id) %>%
  summarise(across(all_of(wl_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  inner_join(select(ww, id, species, nitrate, salinity), by = "id") %>%
  pivot_longer(
    cols      = all_of(wl_cols),
    names_to  = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(wavelength)) %>%
  filter(wavelength >= 400, wavelength <= 750)

# Mean ± SE by nitrate for each species
spec_summ <- spec_plot_long %>%
  group_by(species, nitrate, wavelength) %>%
  summarise(
    mean_refl = mean(reflectance, na.rm = TRUE),
    se_refl   = sd(reflectance, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  )

plot_mean_spectra <- function(spp_label, title) {
  spec_summ %>%
    filter(species == spp_label) %>%
    ggplot(aes(x = wavelength, y = mean_refl,
               color = nitrate, fill = nitrate)) +
    geom_line(linewidth = 0.6) +
    geom_ribbon(aes(ymin = mean_refl - se_refl, ymax = mean_refl + se_refl),
                alpha = 0.15, color = NA) +
    scale_color_viridis_d(option = "plasma", direction = -1) +
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    theme_bw() +
    labs(
      x     = "Wavelength (nm)",
      y     = "Mean reflectance",
      color = "Nitrate (µmol)",
      fill  = "Nitrate (µmol)",
      title = title
    )
}

plot_mean_spectra("u", "Ulva lactuca — Mean Reflectance by Nitrate (Day 9)")
plot_mean_spectra("h", "Hypnea musciformis — Mean Reflectance by Nitrate (Day 9)")

# ── 7. LMM: spectral indices as response variables ───────────────────────────
# Fixed: nitrate + salinity (2025-only, consistent with DSPI 2025 models)
# Random: run_combo only (plant_id collapses in DSPI models; test here too)

fit_spec_lmm <- function(data, response) {
  f <- as.formula(
    paste(response, "~ nitrate + salinity + (1 | run_combo)")
  )
  m <- lmer(f, data = data,
            control = lmerControl(optimizer = "bobyqa"),
            REML = TRUE)
  list(
    model   = m,
    r2      = r.squaredGLMM(m),
    summary = summary(m)
  )
}

indices <- c("ndvi", "sr", "ci_green", "ndre", "rb_ratio")

ulva_spec_results   <- map(indices, ~ fit_spec_lmm(ulva_spec,   .x)) %>% set_names(indices)
hypnea_spec_results <- map(indices, ~ fit_spec_lmm(hypnea_spec, .x)) %>% set_names(indices)

# Print results
for (idx in indices) {
  cat("\n────────────────────────────────────\n")
  cat("Ulva —", idx, "\n")
  cat("R2m:", round(ulva_spec_results[[idx]]$r2[1], 4),
      " R2c:", round(ulva_spec_results[[idx]]$r2[2], 4), "\n")
  print(ulva_spec_results[[idx]]$summary)

  cat("\nHypnea —", idx, "\n")
  cat("R2m:", round(hypnea_spec_results[[idx]]$r2[1], 4),
      " R2c:", round(hypnea_spec_results[[idx]]$r2[2], 4), "\n")
  print(hypnea_spec_results[[idx]]$summary)
}

# ── 8. Spectral index vs. DSPI scatter ───────────────────────────────────────
# Tests decoupling hypothesis: high pigment (CI_green) + flat DSPI

# Use rlc_day == 9 for DSPI day 9 values
ww_spec_d9 <- ww_spec %>%
  filter(rlc_day == 9) %>%
  select(id, species, nitrate, salinity, run_combo,
         dspi_total, ci_green, ndvi, sr)

ggplot(ww_spec_d9, aes(x = ci_green, y = dspi_total,
                        color = nitrate, shape = salinity)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~ species, labeller = as_labeller(c(u = "Ulva", h = "Hypnea"))) +
  scale_color_viridis_d(option = "plasma", direction = -1) +
  theme_bw() +
  labs(
    x     = "Chlorophyll Index Green (CI_green)",
    y     = "DSPI total (µmol O₂ m⁻² day⁻¹)",
    color = "Nitrate (µmol)",
    shape = "Salinity (ppt)",
    title = "Pigment index vs. photosynthetic output by nitrate treatment"
  )

# ── 9. Index distributions by nitrate (box plots) ────────────────────────────

ww_spec %>%
  select(species, nitrate, salinity, run_combo, all_of(indices)) %>%
  pivot_longer(all_of(indices), names_to = "index", values_to = "value") %>%
  ggplot(aes(x = nitrate, y = value, fill = salinity)) +
  geom_boxplot(outlier.size = 1, alpha = 0.7) +
  facet_grid(index ~ species, scales = "free_y",
             labeller = labeller(
               species = c(u = "Ulva", h = "Hypnea")
             )) +
  scale_fill_manual(values = c("18" = "#4DBBEE", "22" = "#E8844A")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x     = "Nitrate (µmol)",
    y     = "Index value",
    fill  = "Salinity (ppt)",
    title = "Spectral indices by nitrate treatment and salinity"
  )
