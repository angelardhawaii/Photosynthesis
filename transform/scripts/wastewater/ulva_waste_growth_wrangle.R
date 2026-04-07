#Wrangling Ulva growth data from Wastewater experiments
#By Angela Richards Donà
#February 2,2026

library(dplyr)
library(tidyverse)
ulva_growth <- read_csv("/Users/angela/src/Photosynthesis/data/wastewater/input/growth_2025/ulva_waste_growth.csv")

#Reframe data as growth periods or intervals
ulva_growth <- ulva_growth %>%
  mutate(
    g_1_5 = d5_weight - i_weight,
    g_5_9 = f_weight  - d5_weight,
    g_1_9 = f_weight  - i_weight
  )

#Classify allocation state by interval
ulva_growth <- ulva_growth %>%
  mutate(
    alloc_1_5 = case_when(
      g_1_5 < 0 ~ "reproduction",
      g_1_5 > 0 ~ "growth",
      TRUE      ~ "neutral"
    ),
    alloc_5_9 = case_when(
      g_5_9 < 0 ~ "reproduction",
      g_5_9 > 0 ~ "growth",
      TRUE      ~ "neutral"
    )
  )

#Identify recovery trajectories
ulva_growth <- ulva_growth %>%
  mutate(
    trajectory = paste(alloc_1_5, "→", alloc_5_9)
  )


#Estimate growth capacity (how much should they have grown?)
interval_growth <- ulva_growth %>%
  select(run, treatment, plant_id,
         g_1_5, g_5_9,
         alloc_1_5, alloc_5_9) %>%
  pivot_longer(
    cols = c(g_1_5, g_5_9),
    names_to = "interval",
    values_to = "growth"
  ) %>%
  pivot_longer(
    cols = c(alloc_1_5, alloc_5_9),
    names_to = "alloc_interval",
    values_to = "allocation"
  ) %>%
  filter(
    (interval == "g_1_5" & alloc_interval == "alloc_1_5") |
      (interval == "g_5_9" & alloc_interval == "alloc_5_9")
  )
#Estimate expected growth per interval
expected_interval_growth <- interval_growth %>%
  filter(allocation == "growth") %>%
  group_by(run, treatment, interval) %>%
  summarise(
    mean_growth = mean(growth, na.rm = TRUE),
    .groups = "drop"
  )
ulva_growth <- ulva_growth %>%
  left_join(
    expected_interval_growth %>%
      pivot_wider(names_from = interval, values_from = mean_growth) %>%
      rename(
        exp_g_1_5 = g_1_5,
        exp_g_5_9 = g_5_9
      ),
    by = c("run", "treatment")
  ) %>%
  mutate(
    loss_1_5 = if_else(g_1_5 < 0, abs(g_1_5), 0),
    loss_5_9 = if_else(g_5_9 < 0, abs(g_5_9), 0),
    
    gross_est_1_5 = g_1_5 + loss_1_5,
    gross_est_5_9 = g_5_9 + loss_5_9,
    
    # if exp_g_1_5 is NA (e.g., nobody grew in that run/treatment/interval),
    # this will propagate NA; that's honest, but you can decide later how to handle it.
    unrealized_1_5 = if_else(g_1_5 < 0, exp_g_1_5 - g_1_5, 0),
    unrealized_5_9 = if_else(g_5_9 < 0, exp_g_5_9 - g_5_9, 0)
  )
ulva_growth %>% 
  select(run, treatment, plant_id, i_weight, d5_weight, f_weight,
         g_1_5, g_5_9, alloc_1_5, alloc_5_9, trajectory,
         loss_1_5, loss_5_9, gross_est_1_5, gross_est_5_9,
         exp_g_1_5, exp_g_5_9, unrealized_1_5, unrealized_5_9) %>% 
  glimpse()

#Plot spaghetti
library(ggplot2)

ulva_long <- ulva_growth %>%
  select(run, treatment, plant_id, i_weight, d5_weight, f_weight) %>%
  pivot_longer(
    cols = c(i_weight, d5_weight, f_weight),
    names_to = "day",
    values_to = "weight"
  ) %>%
  mutate(
    day = recode(day,
                 i_weight = "day_1",
                 d5_weight = "day_5",
                 f_weight  = "day_9"),
    day = factor(day, levels = c("day_1", "day_5", "day_9"))
  )

ggplot(ulva_long, aes(day, weight, group = plant_id)) +
  geom_line(alpha = 0.35) +
  geom_point(size = 1.2, alpha = 0.6) +
  facet_grid(run ~ treatment, scales = "free_y") +
  labs(x = NULL, y = "Weight (g)", title = "Ulva weights through time by run and treatment")

#plot mean trajectory with uncertainty (run x treatment)
ulva_summary <- ulva_long %>%
  group_by(run, treatment, day) %>%
  summarise(
    mean_w = mean(weight, na.rm = TRUE),
    se_w   = sd(weight, na.rm = TRUE) / sqrt(sum(!is.na(weight))),
    n      = sum(!is.na(weight)),
    .groups = "drop"
  )

ggplot(ulva_summary,
       aes(day, mean_w,
           group = treatment,
           color = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean_w - se_w,
        ymax = mean_w + se_w),
    width = 0.08
  ) +
  facet_wrap(~ run, scales = "free_y") +
  labs(
    x = NULL,
    y = "Mean weight (g) ± SE",
    color = "Treatment",
    title = "Mean Ulva trajectory by treatment within each run"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold") +
    scale_color_brewer(palette = "Dark2")
  )
