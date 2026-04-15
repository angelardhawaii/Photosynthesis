# Spectral Reflectance Indices and LMM Analysis
# By Angela Richards Dona
# Wastewater SGD 2025 experiment — Ulva lactuca and Hypnea musciformis
#
# Loads wrangled spectral data from spec_wrangle.R output (spec_norm_long.csv)
# Computes algae-appropriate spectral indices, joins to ww_combo_hsat_dspi,
# runs LMM per index.
#
# Wavelengths chosen from observed absorptance peaks in processed spectra and literature:
#
#   BOTH species:
#     435 nm — Chl a/b + carotenoid Soret absorption band
#     500 nm — Ulva: carotenoid/accessory pigment absorption;
#               Hypnea: B-phycoerythrin Soret peak (observed in spectra)
#     572 nm — Ulva: green window trough (reflectance peak);
#               Hypnea: R-phycoerythrin peak (observed in spectra; shifted from textbook ~560)
#     675 nm — Chl a red absorption maximum (both; less prominent in Ulva)
#     750 nm — NIR reference baseline (instrument max ~750.7 nm)
#
#   Ulva-specific:
#     650 nm — Chl b red absorption peak (dominant red peak in these specimens;
#               observed minimum reflectance ~650–655 nm in raw data)
#
#   Hypnea-specific:
#     625 nm — Phycocyanin absorption (observed ~620–625 nm; shifted from textbook 620)
#
# NOTE on Ulva baseline:
#   The spectrometer range ends at ~750.7 nm. Ulva's red edge reflectance peaks at ~688–691 nm
#   and then DECLINES to 749 nm — the true NIR plateau is not captured. This makes the
#   749 nm baseline correction push the 680–700 nm region negative for Ulva. The processed
#   Ulva absorptance_zero peaks at ~650–655 nm (Chl b) and is near-zero or negative by
#   675 nm. This is real spectral behavior given the limited wavelength range, not a
#   wrangling error. Ratio-based indices using R_750 as a relative reference remain valid.

library(tidyverse)#includes ggplot2, readr, forcats, tibble, purrr, tidyr, dplyr, stringr
library(lme4)
library(lmerTest)
library(MuMIn)

# ── 1. Load wrangled spectral data ────────────────────────────────────────────
# Run spec_wrangle.R first to generate this file

spec_long <- read_csv(
  "/Users/angela/src/Photosynthesis/data/wastewater/transformed/spec_norm_long.csv",
  show_col_types = FALSE
) %>%
  mutate(
    nitrate  = as.factor(nitrate),
    salinity = as.factor(salinity),
    species  = as.character(species)
  )

# ── 2. Compute spectral indices from reflectance ──────────────────────────────
# Indices computed on raw reflectance at nearest available wavelength.
# Both-species indices: chl_ratio, car_chl, rg_ratio
# Ulva-only: acc500, chlb_ratio, chl_ab
# Hypnea-only: pe_b_index (500 nm), pe_index (572 nm), pc_index (625 nm)

spec_wide <- spec_long %>%
  select(id, replicate, wavelength, reflectance, nitrate, salinity, species, run) %>%
  group_by(id, replicate, nitrate, salinity, species, run) %>%
  summarise(
    R_435 = reflectance[which.min(abs(wavelength - 435))],  # Chl a/b + carotenoid (both)
    R_500 = reflectance[which.min(abs(wavelength - 500))],  # Ulva: carotenoid; Hypnea: B-PE
    R_572 = reflectance[which.min(abs(wavelength - 572))],  # Ulva: green trough; Hypnea: R-PE
    R_625 = reflectance[which.min(abs(wavelength - 625))],  # Phycocyanin (Hypnea)
    R_650 = reflectance[which.min(abs(wavelength - 650))],  # Chl b red peak (Ulva ~650-655)
    R_675 = reflectance[which.min(abs(wavelength - 675))],  # Chl a red max (both)
    R_750 = reflectance[which.min(abs(wavelength - 750))],  # NIR baseline
    .groups = "drop"
  ) %>%
  mutate(
    # ── Both-species indices ──────────────────────────────────────────────────
    # Chlorophyll a ratio — primary Chl a proxy; R750/R675 (Gitelson et al.)
    chl_ratio = R_750 / R_675,

    # Carotenoid:Chlorophyll a ratio — increases under nutrient/light stress
    car_chl   = R_675 / R_435,

    # Red:Green ratio — broad pigment load indicator; R675/R572
    rg_ratio  = R_675 / R_572,

    # ── Ulva-specific indices ─────────────────────────────────────────────────
    # 500 nm normalized difference — carotenoid/accessory pigment absorption in Ulva
    # (same formula as pe_b_index; different pigment target in Chlorophyta)
    acc500     = (R_750 - R_500) / (R_750 + R_500),

    # Chl b ratio — observed peak ~650 nm in these specimens; Chl b absent in red algae
    chlb_ratio = R_750 / R_650,

    # Chl a:b ratio — stress indicator; rises when Chl b degrades faster than Chl a
    chl_ab     = R_675 / R_650,

    # ── Hypnea-specific indices (phycobiliproteins) ───────────────────────────
    # B-phycoerythrin index — Soret peak at ~500 nm; observed in Hypnea spectra
    pe_b_index = (R_750 - R_500) / (R_750 + R_500),

    # R-phycoerythrin index — dominant phycobiliprotein; observed peak ~572 nm
    pe_index   = (R_750 - R_572) / (R_750 + R_572),

    # Phycocyanin index — secondary phycobiliprotein; observed ~625 nm
    pc_index   = (R_750 - R_625) / (R_750 + R_625)
  )

# Average across replicates — one row per individual
indices_both   <- c("chl_ratio", "car_chl", "rg_ratio")
indices_ulva   <- c("acc500", "chlb_ratio", "chl_ab")      # Ulva-only: 500, 650 nm peaks
indices_hypnea <- c("pe_b_index", "pe_index", "pc_index")  # Hypnea-only: 500, 572, 625 nm
indices_all    <- c(indices_both, indices_ulva, indices_hypnea)

spec_idx <- spec_wide %>%
  group_by(id, nitrate, salinity, species, run) %>%
  summarise(
    across(all_of(indices_all), mean, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )

# ── 3. Join to ww_combo_hsat_dspi (2025 data only) ───────────────────────────

ww <- read_csv(
  "/Users/angela/src/Photosynthesis/data/wastewater/transformed/ww_combo_hsat_dspi.csv",
  show_col_types = FALSE
) %>%
  filter(temp == 22) %>%
  mutate(
    id        = tolower(trimws(id)),
    nitrate   = as.factor(nitrate),
    salinity  = as.factor(salinity),
    run_combo = as.factor(run_combo),
    plant_id  = as.factor(plant_id)
  )

ww_spec <- ww %>%
  inner_join(spec_idx, by = "id")

message("Rows after join: ", nrow(ww_spec))
message("run_combos in joined data: ",
        paste(sort(unique(as.character(ww_spec$run_combo))), collapse = ", "))

# ── 4. Species subsets ────────────────────────────────────────────────────────

ulva_spec   <- filter(ww_spec, species.x == "u")
hypnea_spec <- filter(ww_spec, species.x == "h")

# ── 5. Mean absorptance spectra by nitrate (visual check) ────────────────────
# Annotate key pigment wavelengths for reference

# Species-specific peak annotations (observed peak wavelengths)
pigment_lines_ulva <- tibble(
  wavelength = c(435, 500, 572, 650, 675),
  label      = c("Chl a/b + Car", "Carotenoid", "Green trough", "Chl b", "Chl a")
)

pigment_lines_hypnea <- tibble(
  wavelength = c(435, 500, 572, 625, 675),
  label      = c("Chl a/b + Car", "B-PE", "R-PE", "PC", "Chl a")
)

spec_summ <- spec_long %>%
  group_by(species, nitrate, wavelength) %>%
  summarise(
    mean_abs = mean(absorptance_zero, na.rm = TRUE),
    se_abs   = sd(absorptance_zero,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  filter(wavelength >= 400, wavelength <= 750)

plot_mean_spectra <- function(spp_label, title, peak_lines, caption_text) {
  spec_summ %>%
    filter(species == spp_label) %>%
    ggplot(aes(x = wavelength, y = mean_abs,
               color = nitrate, fill = nitrate)) +
    geom_vline(data = peak_lines,
               aes(xintercept = wavelength),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.6) +
    geom_ribbon(aes(ymin = mean_abs - se_abs, ymax = mean_abs + se_abs),
                alpha = 0.15, color = NA) +
    scale_color_viridis_d(option = "plasma", direction = -1) +
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    theme_bw() +
    labs(
      x       = "Wavelength (nm)",
      y       = "Mean zero-baseline absorptance (SNV)",
      color   = "Nitrate (µmol)",
      fill    = "Nitrate (µmol)",
      title   = title,
      caption = caption_text
    )
}

plot_mean_spectra(
  "u", "Ulva lactuca — Mean Absorptance by Nitrate (Day 9)",
  pigment_lines_ulva,
  "Dashed lines: 435 nm Chl a/b+Car | 500 nm Carotenoid | 572 nm Green trough | 650 nm Chl b | 675 nm Chl a"
)
plot_mean_spectra(
  "h", "Hypnea musciformis — Mean Absorptance by Nitrate (Day 9)",
  pigment_lines_hypnea,
  "Dashed lines: 435 nm Chl a/b+Car | 500 nm B-PE | 572 nm R-PE | 625 nm PC | 675 nm Chl a"
)

# ── 6. LMM: both-species indices ─────────────────────────────────────────────
# Fixed: nitrate + salinity  |  Random: run_combo

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

ulva_spec_results   <- map(indices_both, ~ fit_spec_lmm(ulva_spec,   .x)) %>%
  set_names(indices_both)
hypnea_spec_results <- map(indices_both, ~ fit_spec_lmm(hypnea_spec, .x)) %>%
  set_names(indices_both)

# Print both-species results
for (idx in indices_both) {
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

# ── 7a. LMM: Ulva-only indices (500 nm + Chl b 650 nm) ──────────────────────
# acc500: carotenoid/accessory at 500 nm; chlb_ratio and chl_ab: Chl b at 650 nm

ulva_chlb_results <- map(indices_ulva, ~ fit_spec_lmm(ulva_spec, .x)) %>%
  set_names(indices_ulva)

for (idx in indices_ulva) {
  cat("\n────────────────────────────────────\n")
  cat("Ulva —", idx, "\n")
  cat("R2m:", round(ulva_chlb_results[[idx]]$r2[1], 4),
      " R2c:", round(ulva_chlb_results[[idx]]$r2[2], 4), "\n")
  print(ulva_chlb_results[[idx]]$summary)
}

# ── 7b. LMM: Hypnea-only phycobiliprotein indices ────────────────────────────
# pe_b_index (500 nm), pe_index (572 nm), pc_index (625 nm): Hypnea only

hypnea_phyco_results <- map(indices_hypnea, ~ fit_spec_lmm(hypnea_spec, .x)) %>%
  set_names(indices_hypnea)

for (idx in indices_hypnea) {
  cat("\n────────────────────────────────────\n")
  cat("Hypnea —", idx, "\n")
  cat("R2m:", round(hypnea_phyco_results[[idx]]$r2[1], 4),
      " R2c:", round(hypnea_phyco_results[[idx]]$r2[2], 4), "\n")
  print(hypnea_phyco_results[[idx]]$summary)
}

# ── 8. Index distributions by nitrate (box plots) ────────────────────────────

ww_spec %>%
  select(species.x, nitrate.x, salinity.x, run_combo, all_of(indices_all)) %>%
  rename(species = species.x, nitrate = nitrate.x, salinity = salinity.x) %>%
  pivot_longer(all_of(indices_all), names_to = "index", values_to = "value") %>%
  filter(
    !(index %in% indices_hypnea & species == "u"),  # phycobiliprotein indices: Hypnea only
    !(index %in% indices_ulva   & species == "h")   # Chl b indices: Ulva only
  ) %>%
  ggplot(aes(x = nitrate, y = value, fill = salinity)) +
  geom_boxplot(outlier.size = 1, alpha = 0.7) +
  facet_grid(index ~ species, scales = "free_y",
             labeller = labeller(species = c(u = "Ulva", h = "Hypnea"))) +
  scale_fill_manual(values = c("18" = "#4DBBEE", "22" = "#E8844A")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x        = "Nitrate (µmol)",
    y        = "Index value",
    fill     = "Salinity (ppt)",
    title    = "Spectral indices by nitrate treatment and salinity",
    caption  = paste("acc500/chlb_ratio/chl_ab: Ulva only (500 nm carotenoid, 650 nm Chl b);",
                     "pe_b/pe/pc indices: Hypnea only (phycobiliproteins absent in Ulva)")
  )

# ── 9. Chl ratio vs. DSPI scatter (pigment-photosynthesis decoupling) ─────────

ww_spec %>%
  filter(rlc_day == 9) %>%
  ggplot(aes(x = chl_ratio, y = dspi_total,
             color = nitrate.x, shape = salinity.x)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~ species.x,
             labeller = as_labeller(c(u = "Ulva", h = "Hypnea"))) +
  scale_color_viridis_d(option = "plasma", direction = -1) +
  theme_bw() +
  labs(
    x     = "Chlorophyll ratio (R750/R675)",
    y     = "DSPI total (µmol O₂ m⁻² day⁻¹)",
    color = "Nitrate (µmol)",
    shape = "Salinity (ppt)",
    title = "Chlorophyll index vs. photosynthetic output by nitrate treatment"
  )
