# Script to clean PAM data for midday depression
# Only effective quantum yield and F values needed for each ID
# Three readings were taken for each individual every hour from
# 08:00 to 18:00
# By Angela Richards Donà
# June 18, 2025


#load packages
library(dplyr)
library(tidyverse)
library(hms)
library(lubridate)

midday <- read.csv("../data/midday_sad_chondria.csv")

cleaned_data <- midday %>%
  filter(Type %in% c("SLCS", "FO")) %>%       # Find and keep only rows with SLCS and FO for column Type
  select(-Datetime, -PAR, -ETR, -FvFm, -no)   # Get rid of columns indicated

cleaned_data <- cleaned_data %>%
  mutate(row_num = row_number()) %>%
  mutate(ID = lag(FmPrime, default = NA)) %>%                # get ID from the previous row
  filter(Type != "SLCS" & F != "Light Curve start") %>%       # keep only data rows
  select(Date, Time, Type, F, FmPrime, YII, ID) 

cleaned_data <- cleaned_data %>%
  select(-Type) %>%
  rename("F0" = "F", "Fm" = "FmPrime")

# Convert Time to hms
cleaned_data <- cleaned_data %>%
  mutate(Time = as_hms(Time))
cleaned_data <- cleaned_data %>%
  mutate(Time = as_hms(as.numeric(Time) - 3600))  # subtract 1 hour = 3600 seconds

# Create target hourly time points
target_hours <- sprintf("%02d:00:00", 8:18) # 08:00 to 18:00

# Create time windows using hms-compatible arithmetic
time_windows <- tibble(
  target_time = as_hms(target_hours),
  window_start = as_hms(target_time) - 5 * 60,  # subtract 5 minutes
  window_end   = as_hms(target_time) + 23 * 60  # add 23 minutes
)

# Join each row with its matching time window
hourly_data <- cleaned_data %>%
  rowwise() %>%
  mutate(hour_group = time_windows$target_time[
    which(Time >= time_windows$window_start & Time <= time_windows$window_end)][1]
  ) %>%
  ungroup()

hourly_data$Date <- mdy(hourly_data$Date) 
hourly_data$ID <- as.factor(hourly_data$ID)

hourly_data <- hourly_data %>%
  mutate(hour_label = format(as.POSIXct(hour_group, format = "%H:%M:%S"), "%H:%M"))

hourly_data$hour_label <- factor(hourly_data$hour_label, levels = sprintf("%02d:00", 8:18))


#plot the data
midday_plot <- hourly_data %>%
  ggplot(aes(x = hour_label, y = YII, color = ID)) +
  theme_bw() +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(vars(ID), nrow = 2) +
  labs(x = "Hour") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
midday_plot
