# Wrangle spectral reflectance data
# By Angela Richards Dona
# 07/21/25
#The script is meant to open, clean, and plot some spectral reflectance data in prep for doing some analysis

#load some helpful libraries
library(dplyr)
library(tidyr)
library(prospectr)
library(hyperSpec)
library(pavo)
library(spectrolab)
library(tidyverse)
library(purrr)
library (tibble)
library(janitor)
library(readr)
library(signal)
library(scales)
library(grid)

# read in files, renames the columns (ids) to include replicate numbers, appends _4 to any in need of replicate
#stacks the files
#these files have been modified in script clean_script.R to remove unwanted wavelengths and to clean names

files <- list.files(
  path       = "/Users/angela/src/Photosynthesis/data/wastewater/input/spec/",
  pattern    = "\\.csv$",
  full.names = TRUE
)

# map & bind into one tibble
spec_waste_all <- files %>%
  map_dfr(
    ~ read_csv(.x, show_col_types = FALSE, name_repair = "minimal")
  )

#________________________________________________________________________________

# Make sure reflectance columns are numeric
wl_cols <- grep("^\\d+\\.\\d+$", names(spec_waste_all), value = TRUE)

spec_all <- spec_waste_all %>%
  mutate(across(all_of(wl_cols), as.numeric))

# Extract reflectance matrix
refl_mat <- spec_all %>%
  select(all_of(wl_cols)) %>%
  as.matrix()

smooth_sg <- function(x, w = 11, p = 2, m = 0) {
  sgolayfilt(x, p = p, n = w, m = m)
}

#Smooth each row, keeping the same dimensions
sm_mat <- t(apply(refl_mat,
                  1,
                  smooth_sg,
                  w = w,    # window
                  p = 2,    # polynomial order
                  m = 0))   # derivative = 0 (pure smoothing)

smoothed_wl <- as.numeric(wl_cols)

# sanity check: these must be equal
stopifnot(length(smoothed_wl) == ncol(sm_mat))

# reassemble your smoothed data.frame with the correct names
spec_smoothed <- bind_cols(
  tibble(id = spec_all$id),
  as_tibble(sm_mat)
) %>%
  set_names(c("id", as.character(smoothed_wl)))


# pivot back to long for plotting
spec_long <- spec_smoothed %>%
  pivot_longer(
    -id,
    names_to  = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(wavelength))

#add columns for absorptance (%), species, and replicates
spec_long <- spec_long %>%
  mutate(absorptance = 100 - reflectance) %>%
  mutate(species = str_sub(id, 1,1)) %>%
  mutate(replicate = str_sub(id, 6,6)) %>%
  mutate(id = str_sub(id, 1, 4)) %>%
  mutate(run = str_sub(id, 2, 2))
  

#SNV takes each spectrum, subtracts its own mean and divides by its own standard deviation:
#\text{SNV}(x_i) = \frac{x_i - \bar{x}}{\mathrm{sd}(x)}
#so all spectra have zero mean and unit variance, making them directly comparable.
spec_snv <- spec_long %>%
  group_by(id) %>%
  mutate(
    # compute each spectrum’s mean and sd
    mean_abs   = mean(absorptance),
    sd_abs     = sd(absorptance),
    # SNV: subtract the mean, then divide by the sd
    absorptance_snv = (absorptance - mean_abs) / sd_abs
  ) %>%
  ungroup() %>%
  select(-mean_abs, -sd_abs)

spec_mean <- spec_snv %>%
  group_by(id, species, wavelength) %>%
  summarise(
    absorptance_snv = mean(absorptance_snv),
    run = first(run),
    replicate = first(replicate),
    .groups = "drop"
  )

spec_norm <- spec_mean %>%
  group_by(id) %>%
  mutate(
    baseline749 = absorptance_snv[which.min(abs(wavelength - 700.44)) ],
    absorptance_zero = absorptance_snv - baseline749
  ) %>%
  ungroup()
####
#make a column for replicate using the 3rd and 4th digits from id column
spec_norm$replicate <- as.factor(spec_norm$replicate)
spec_norm <-spec_norm %>%
  mutate(id_num = as.numeric(str_sub(id, 3, 4)))

# Create treatment and salinity columns based on id number
treatment_pattern <- c(7, 8, 9, 10) # 7, 8, 9, 10 are the treatment numbers

spec_norm <- spec_norm %>%
  mutate(
    treatment = as.factor(treatment_pattern[(id_num - 1) %% 4 + 1])
  )

spec_norm <- spec_norm %>%
  mutate(
    salinity = as.factor(if_else(((id_num - 1) %/% 4) %% 2 == 0, 18, 22))
  )

#Subset for the two species to plot
spec_u <- spec_norm %>%
  dplyr::filter(species == "u")
spec_h <- spec_norm %>%
  dplyr::filter(species == "h")


# plot raw data
ggplot(spec_u, aes(wavelength, absorptance_zero, color = id)) +
  geom_line(linewidth = 0.3) +
  xlim(400,700) +
  ylim(0, 4) +
  facet_wrap(salinity~treatment, nrow = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "Wavelength (nm)",
    y = "Mean zero-baseline absorptance",
    title = "Ulva Spectrum"
  )
ggplot(spec_h, aes(wavelength, absorptance_zero, color = id)) +
  geom_line(linewidth = 0.3) +
  xlim(400,700) +
  ylim(0, 4) +
  facet_wrap(treatment ~ salinity, nrow = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "Wavelength (nm)",
    y = "Mean zero-baseline absorptance",
    title = "Hypnea Spectrum"
  )

#Get mean +- standard dev to see group level treands with variability
spec_u %>%
  group_by(wavelength, treatment, salinity) %>%
  summarise(
    mean_abs = mean(absorptance_zero, na.rm = TRUE),
    sd_abs = sd(absorptance_zero, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = wavelength, y = mean_abs)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs),
              alpha = 0.3) +
  facet_grid(salinity ~ treatment) +
  theme_bw()

spec_h %>%
  group_by(wavelength, treatment, salinity) %>%
  summarise(
    mean_abs = mean(absorptance_zero, na.rm = TRUE),
    sd_abs = sd(absorptance_zero, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = wavelength, y = mean_abs)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs),
              alpha = 0.3) +
  facet_grid(salinity ~ treatment) +
  theme_bw()

#One graph for means of each group of treatment salinity combo (8 total)
spec_means_u <- spec_u %>%
  group_by(wavelength, treatment, salinity) %>%
  summarise(mean_abs = mean(absorptance_zero, na.rm = TRUE), .groups = "drop")

spec_means_u <- spec_means_u %>%
  mutate(group = paste0("T", treatment, "_S", salinity))  # create combined group label

ggplot(spec_means_u, aes(x = wavelength, y = mean_abs, color = group)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(
    x = "Wavelength (nm)",
    y = "Mean zero-baseline absorptance",
    color = "treatment/salinity",
    title = "Ulva lactuca Mean Spectra by Treatment and Salinity"
  ) +
  xlim(400, 700) +
  ylim(0, 3)

spec_means_h <- spec_h %>%
  group_by(wavelength, treatment, salinity) %>%
  summarise(mean_abs = mean(absorptance_zero, na.rm = TRUE), .groups = "drop")

spec_means_h <- spec_means_h %>%
  mutate(group = paste0("T", treatment, "_S", salinity))  # create combined group label

ggplot(spec_means_h, aes(x = wavelength, y = mean_abs, color = group)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(
    x = "Wavelength (nm)",
    y = "Mean zero-baseline absorptance",
    color = "treatment/salinity",
    title = "Hypnea musciformis Mean Spectra by Treatment and Salinity"
  ) +
  xlim(400, 700) +
  ylim(0, 1.75)

# add background to mean spectrum
spec_means_u <- spec_u %>%
  group_by(wavelength, treatment, salinity) %>%
  summarise(mean_abs = mean(absorptance_zero, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = paste0("T", treatment, "_S", salinity))

# Step 2: Create visible spectrum gradient (violet to red)
n_colors <- 200
hues <- seq(270/360, 0, length.out = n_colors)  # violet to red in hue space
vis_colors <- grDevices::hcl(h = hues * 360, c = 100, l = 65, alpha = 0.4)  # 40% opacity
vis_spec_gradient <- matrix(vis_colors, nrow = 1)

# Step 3: Convert to raster graphic
g <- rasterGrob(vis_spec_gradient, width = unit(1, "npc"), height = unit(1, "npc"),
                interpolate = TRUE)

# Step 4: Plot
ggplot(spec_means_u, aes(x = wavelength, y = mean_abs, linetype = group)) +
  annotation_custom(g, xmin = 400, xmax = 700, ymin = -Inf, ymax = Inf) +
  geom_line(color = "black", linewidth = 0.5) +
  scale_linetype_manual(values = c(
    "T7_S18" = "solid",
    "T8_S18" = "dashed",
    "T9_S18" = "dotted",
    "T10_S18" = "dotdash",
    "T7_S22" = "longdash",
    "T8_S22" = "twodash",
    "T9_S22" = "44",          # encode custom dash pattern as string
    "T10_S22" = "11"      # also as string (must be valid length: 2, 4, 6, or 8)
  )) +
  coord_cartesian(xlim = c(400, 700), ylim = c(0, 4), expand = FALSE) +
  theme_bw() +
  labs(
    x = "Wavelength (nm)",
    y = "Mean zero-baseline absorptance",
    linetype = "Treatment / Salinity",
    title = "Mean Spectra with Visible Spectrum Background"
  )

