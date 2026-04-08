# Comparison script: DSPI computed using the method from Richards Dona et al. (2023)
# Purpose: curiosity check — applies the 2023 paper's per-period average-daily Hsat × Pmax
# approach to the current combined 2021/2025 dataset. No models — plots only.
# By Angela Richards Donà
# April 2026
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(lubridate)
library(hms)
library(ggpubr)
library(paletteer)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Load data (same source as waste_sgd_hsat_dspi_lmm.R)
# ─────────────────────────────────────────────────────────────────────────────
dat <- read_csv("/Users/angela/src/Photosynthesis/data/wastewater/transformed/ww_combo_hsat_dspi.csv")

# Remove individual that died during experiment
dat <- dat %>% filter(id != "ul09")

# Basic type fixes
dat <- dat %>%
  mutate(
    start_date = as.Date(start_date),
    year       = lubridate::year(start_date),
    across(c(sal_bin, temp, nitrate, run_combo, salinity), as.factor)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 2) Compute DSPI the 2023-paper way:
#    For each of the three irradiance periods, divide total-period Hsat by the
#    number of days in that period to get an average-daily Hsat, then multiply
#    by the corresponding Pmax (in µmol O2 m⁻² min⁻¹).  Sum the three periods.
#    This mirrors: dspi_pmax_total = Σ (avg_daily_hsat_pi × pmax_pi) in the
#    old mean_supersaturation_by_period.csv pipeline.
#
#    Period day counts (matches waste_hsat_dspi_build.R):
#      with_day5 scheme → P1 = days 1–3 (3 d), P2 = days 4–6 (3 d), P3 = days 7–9 (3 d)
#      no_day5  scheme  → P1 = days 1–4 (4 d), P2 = days 5–9 (5 d), P3 = skipped (0 d)
# ─────────────────────────────────────────────────────────────────────────────
dat <- dat %>%
  mutate(
    p1_days = if_else(scheme == "no_day5", 4L, 3L),
    p2_days = if_else(scheme == "no_day5", 5L, 3L),
    p3_days = if_else(scheme == "with_day5", 3L, 0L),

    # Average daily Hsat (min day⁻¹) for each period
    avg_hsat_p1 = p1_hsat_min / p1_days,
    avg_hsat_p2 = p2_hsat_min / p2_days,
    avg_hsat_p3 = if_else(p3_days > 0, p3_hsat_min / p3_days, 0),

    # Pmax to use per period (µmol O₂ m⁻² min⁻¹); mirrors build script logic
    pmax_p2 = if_else(scheme == "no_day5", pmax_min_day9, pmax_min_day5),

    # Per-period contribution (µmol O₂ m⁻²) — one representative day × Pmax
    dspi_old_p1 = avg_hsat_p1 * pmax_min_day1,
    dspi_old_p2 = avg_hsat_p2 * pmax_p2,
    dspi_old_p3 = avg_hsat_p3 * pmax_min_day9,

    # Sum → convert µmol to mol m⁻²
    dspi_old_mol = round((dspi_old_p1 + dspi_old_p2 + dspi_old_p3) / 1e6, 3)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 3) Apply the same dataset filter as the LMM script
#    Keep: 2021–2022 data at nitrate == 80 µmol only; all 2025 data.
# ─────────────────────────────────────────────────────────────────────────────
dat_filtered <- dat %>%
  filter(
    (year %in% c(2021, 2022) & nitrate == "80") |
      year == 2025
  )

# Treatment-graph labels (same order as LMM script)
dat_filtered <- dat_filtered %>%
  mutate(treatment_graph = case_when(
    treat_letter == "a" ~ "1) 0.5 μmol",
    treat_letter == "b" ~ "3) 14 μmol",
    treat_letter == "c" ~ "4) 27 μmol",
    treat_letter == "d" ~ "6) 53 μmol",
    treat_letter == "e" ~ "7) 53 μmol",
    treat_letter == "f" ~ "9) 80 μmol",
    treat_letter == "g" ~ "2) 0.5 μmol",
    treat_letter == "h" ~ "5) 37 μmol",
    treat_letter == "i" ~ "8) 53 μmol",
    treat_letter == "j" ~ "10) 86 μmol",
    treat_letter == "k" ~ "11) 80 μmol",
    treat_letter == "l" ~ "12) 245 μmol",
    treat_letter == "m" ~ "13) 748 μmol",
    treat_letter == "n" ~ "14) 2287 μmol"
  ))

# Species subsets
ulva  <- dat_filtered %>% filter(species == "u")
hypnea <- dat_filtered %>% filter(species == "h")

# ─────────────────────────────────────────────────────────────────────────────
# 4) Quick summary: treatment-level means (matches Table output in old script)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Ulva: mean DSPI (2023 method) by nitrate treatment ──\n")
ulva %>%
  group_by(nitrate) %>%
  summarise(mean_dspi_old = round(mean(dspi_old_mol, na.rm = TRUE), 3),
            n = n(), .groups = "drop") %>%
  print()

cat("\n── Hypnea: mean DSPI (2023 method) by nitrate treatment ──\n")
hypnea %>%
  group_by(nitrate) %>%
  summarise(mean_dspi_old = round(mean(dspi_old_mol, na.rm = TRUE), 3),
            n = n(), .groups = "drop") %>%
  print()

# ─────────────────────────────────────────────────────────────────────────────
# 5) Plots — styled to match the 2023 paper aesthetic
#    (coloured by temperature, reference line at 5.0 mol m⁻², y-axis 0–12)
# ─────────────────────────────────────────────────────────────────────────────

# -- 5a) Ulva ----------------------------------------------------------------
plot_ulva_old <- ulva %>%
  ggplot(aes(x = treatment_graph, y = dspi_old_mol, color = temp)) +
  geom_boxplot(linewidth = 0.4, outlier.shape = NA) +
  geom_point(alpha = 0.75, size = 2,
             position = position_jitter(width = 0.25), show.legend = TRUE) +
  geom_hline(yintercept = 5.0, color = "red", linewidth = 0.5, alpha = 0.6,
             linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  labs(
    x     = "nitrate treatment",
    y     = "DSPI (mol O₂ m⁻²)",
    title = "Ulva lactuca — DSPI (2023 paper method)",
    subtitle = "Red dashed line = 5.0 mol published threshold",
    color = "temp (°C)"
  ) +
  scale_x_discrete(labels = function(x) sub("^\\d+\\) ", "", x)) +
  ylim(0, 12) +
  scale_color_paletteer_d("MetBrewer::Cross") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 11),
    axis.title   = element_text(size = 13),
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )
plot(plot_ulva_old)

# -- 5b) Hypnea --------------------------------------------------------------
plot_hypnea_old <- hypnea %>%
  ggplot(aes(x = treatment_graph, y = dspi_old_mol, color = temp)) +
  geom_boxplot(linewidth = 0.4, outlier.shape = NA) +
  geom_point(alpha = 0.75, size = 2,
             position = position_jitter(width = 0.25), show.legend = TRUE) +
  geom_hline(yintercept = 5.0, color = "red", linewidth = 0.5, alpha = 0.6,
             linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  labs(
    x     = "nitrate treatment",
    y     = "DSPI (mol O₂ m⁻²)",
    title = "Hypnea musciformis — DSPI (2023 paper method)",
    subtitle = "Red dashed line = 5.0 mol published threshold",
    color = "temp (°C)"
  ) +
  scale_x_discrete(labels = function(x) sub("^\\d+\\) ", "", x)) +
  ylim(0, 12) +
  scale_color_paletteer_d("MetBrewer::Cross") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 11),
    axis.title   = element_text(size = 13),
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )
plot(plot_hypnea_old)

# -- 5c) Both species side-by-side (faceted) ---------------------------------
plot_both_old <- dat_filtered %>%
  ggplot(aes(x = treatment_graph, y = dspi_old_mol, color = temp)) +
  geom_boxplot(linewidth = 0.4, outlier.shape = NA) +
  geom_point(alpha = 0.75, size = 2,
             position = position_jitter(width = 0.25), show.legend = TRUE) +
  geom_hline(yintercept = 5.0, color = "red", linewidth = 0.5, alpha = 0.6,
             linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  facet_wrap(~ species,
             labeller = labeller(species = c("h" = "Hypnea musciformis",
                                             "u" = "Ulva lactuca"))) +
  labs(
    x     = "nitrate treatment",
    y     = "DSPI (mol O₂ m⁻²)",
    title = "DSPI — 2023 paper method applied to 2021 + 2025 dataset",
    subtitle = "Diamonds = treatment means  |  Red dashed line = 5.0 mol",
    color = "temp (°C)"
  ) +
  scale_x_discrete(labels = function(x) sub("^\\d+\\) ", "", x)) +
  ylim(0, 12) +
  scale_color_paletteer_d("MetBrewer::Cross") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 10),
    strip.text   = element_text(face = "italic", size = 12),
    axis.title   = element_text(size = 13),
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )
plot(plot_both_old)

# -- 5d) Side-by-side comparison: old method vs. new dspi_day_mean -----------
#    Requires computing dspi_day_mean on the filtered data for comparison.
#    run_days ≈ 8 elapsed days (day 1 to day 9).
dat_filtered <- dat_filtered %>%
  mutate(dspi_day_mean = round((dspi_total / 1e6) / 8, 3))

comparison_long <- dat_filtered %>%
  select(species, treatment_graph, temp, dspi_old_mol, dspi_day_mean) %>%
  pivot_longer(
    cols      = c(dspi_old_mol, dspi_day_mean),
    names_to  = "method",
    values_to = "dspi"
  ) %>%
  mutate(method = recode(method,
                         dspi_old_mol   = "2023 paper method\n(avg-daily Hsat × Pmax)",
                         dspi_day_mean  = "Current method\n(total DSPI ÷ run days)"))

plot_comparison <- comparison_long %>%
  ggplot(aes(x = treatment_graph, y = dspi, color = temp)) +
  geom_boxplot(linewidth = 0.35, outlier.shape = NA) +
  geom_point(alpha = 0.6, size = 1.5,
             position = position_jitter(width = 0.2), show.legend = TRUE) +
  geom_hline(yintercept = 5.0, color = "red", linewidth = 0.4, alpha = 0.5,
             linetype = "dashed") +
  facet_grid(method ~ species,
             labeller = labeller(
               species = c("h" = "Hypnea musciformis", "u" = "Ulva lactuca")
             )) +
  labs(
    x     = "nitrate treatment",
    y     = "DSPI (mol O₂ m⁻²)",
    title = "DSPI method comparison: 2023 paper vs. current approach",
    color = "temp (°C)"
  ) +
  scale_x_discrete(labels = function(x) sub("^\\d+\\) ", "", x)) +
  scale_color_paletteer_d("MetBrewer::Cross") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 9),
    strip.text   = element_text(face = "italic", size = 11),
    axis.title   = element_text(size = 12),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )
plot(plot_comparison)

# ─────────────────────────────────────────────────────────────────────────────
# Optional: save plots
# ─────────────────────────────────────────────────────────────────────────────
# ggsave("dspi_2023method_ulva.png",     plot = plot_ulva_old,    path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots", width = 8, height = 5, dpi = 300)
# ggsave("dspi_2023method_hypnea.png",   plot = plot_hypnea_old,  path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots", width = 8, height = 5, dpi = 300)
# ggsave("dspi_2023method_both.png",     plot = plot_both_old,    path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots", width = 10, height = 5, dpi = 300)
# ggsave("dspi_method_comparison.png",   plot = plot_comparison,  path = "/Users/angela/src/Photosynthesis/analyze/wastewater/plots", width = 12, height = 7, dpi = 300)
