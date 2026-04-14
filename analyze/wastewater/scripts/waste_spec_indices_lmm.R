# Spectral Reflectance Indices and LMM Analysis
# By Angela Richards Dona
# Wastewater SGD 2025 experiment — Ulva lactuca and Hypnea musciformis
#
# Loads wrangled spectral data from spec_wrangle.R output (spec_norm_long.csv)
# Computes algae-appropriate spectral indices, joins to ww_combo_hsat_dspi,
# runs LMM per index.
#
# Wavelengths chosen for marine macroalgae pigments:
#   435 nm — Chl a + Chl b + carotenoid blue absorption peak
#   560 nm — Phycoerythrin peak (Hypnea); green reflectance peak (Ulva)
#   620 nm — Phycocyanin absorption (Hypnea)
#   675 nm — Chl a red absorption maximum
#   750 nm — NIR reference baseline

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
# Hypnea-only indices: pe_index (phycoerythrin), pc_index (phycocyanin)

spec_wide <- spec_long %>%
  select(id, replicate, wavelength, reflectance, nitrate, salinity, species, run) %>%
  group_by(id, replicate, nitrate, salinity, species, run) %>%
  summarise(
    R_435 = reflectance[which.min(abs(wavelength - 435))],  # Chl a/b + carotenoid
    R_560 = reflectance[which.min(abs(wavelength - 560))],  # Phycoerythrin / green peak
    R_620 = reflectance[which.min(abs(wavelength - 620))],  # Phycocyanin
    R_675 = reflectance[which.min(abs(wavelength - 675))],  # Chl a red absorption max
    R_750 = reflectance[which.min(abs(wavelength - 750))],  # NIR baseline
    .groups = "drop"
  ) %>%
  mutate(
    # ── Both-species indices ──────────────────────────────────────────────────
    # Chlorophyll a ratio — primary Chl a proxy for macroalgae
    # Gitelson et al.; widely used for aquatic macrophytes
    chl_ratio = R_750 / R_675,

    # Carotenoid:Chlorophyll ratio — increases under nutrient/light stress
    # High value = relatively more carotenoid than Chl a
    car_chl   = R_675 / R_435,

    # Red:Green ratio — broad pigment load indicator
    # Reflects overall absorptance balance between red and green regions
    rg_ratio  = R_675 / R_560,

    # ── Hypnea-specific indices (phycobiliproteins) ───────────────────────────
    # Phycoerythrin index — dominant accessory pigment in red algae
    # Absorbs strongly at ~560 nm; not present in Ulva
    pe_index  = (R_750 - R_560) / (R_750 + R_560),

    # Phycocyanin index — secondary phycobiliprotein in red algae
    pc_index  = (R_750 - R_620) / (R_750 + R_620)
  )

# Average across replicates — one row per individual
indices_both  <- c("chl_ratio", "car_chl", "rg_ratio")
indices_hypnea <- c("pe_index", "pc_index")
indices_all   <- c(indices_both, indices_hypnea)

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

pigment_lines <- tibble(
  wavelength = c(435, 560, 620, 675),
  label      = c("Chl a/b + Car", "Phycoerythrin", "Phycocyanin", "Chl a")
)

spec_summ <- spec_long %>%
  group_by(species, nitrate, wavelength) %>%
  summarise(
    mean_abs = mean(absorptance_zero, na.rm = TRUE),
    se_abs   = sd(absorptance_zero,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  filter(wavelength >= 400, wavelength <= 750)

plot_mean_spectra <- function(spp_label, title) {
  spec_summ %>%
    filter(species == spp_label) %>%
    ggplot(aes(x = wavelength, y = mean_abs,
               color = nitrate, fill = nitrate)) +
    geom_vline(data = pigment_lines,
               aes(xintercept = wavelength),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.6) +
    geom_ribbon(aes(ymin = mean_abs - se_abs, ymax = mean_abs + se_abs),
                alpha = 0.15, color = NA) +
    scale_color_viridis_d(option = "plasma", direction = -1) +
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    theme_bw() +
    labs(
      x     = "Wavelength (nm)",
      y     = "Mean zero-baseline absorptance (SNV)",
      color = "Nitrate (µmol)",
      fill  = "Nitrate (µmol)",
      title = title,
      caption = "Dashed lines: 435 nm Chl a/b+Car | 560 nm PE | 620 nm PC | 675 nm Chl a"
    )
}

plot_mean_spectra("u", "Ulva lactuca — Mean Absorptance by Nitrate (Day 9)")
plot_mean_spectra("h", "Hypnea musciformis — Mean Absorptance by Nitrate (Day 9)")

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

# ── 7. LMM: Hypnea-only phycobiliprotein indices ─────────────────────────────
# pe_index and pc_index are only meaningful for Hypnea musciformis

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
  mutate(
    algae_group = if_else(index %in% indices_hypnea & species == "u",
                          "Hypnea only", "both")
  ) %>%
  filter(!(index %in% indices_hypnea & species == "u")) %>%  # drop PE/PC for Ulva
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
    caption  = "pe_index and pc_index shown for Hypnea only (phycobiliproteins absent in Ulva)"
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
