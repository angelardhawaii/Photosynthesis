# Irradiance data cleaner, removes unnecessary columns,  datapoints with zero irradiance removed secondarily
# By Angela Richards Donà
# September 2024

# load libraries
library(dplyr)
library(tidyr)
library(janitor)
library(readr)
library(lubridate)
library(ggplot2)
library(hms)
library(zoo)

#load files and concatenate
files <- list.files("../data/midway_2024/input/licor", pattern = "\\.csv$", full.names = TRUE)

e_all <- files %>%
  lapply(read_csv) %>%
  bind_rows()
glimpse(e_all)
head(paste(e_all$Date, as.character(e_all$Time)))

#make one datetime column from the Date and Time columns
epar_data <- e_all %>%
  mutate(
    datetime = mdy_hms(paste(Date, as.character(Time)))
  )

epar_long <- epar_data %>%
  pivot_longer(cols = c(Epar_canopy, Epar_under),
               names_to = "plant_part",
               values_to = "irradiance")

ggplot(epar_long, aes(x = datetime, y = irradiance, color = plant_part)) +
  geom_line() +
  labs(x = "Time", y = "Irradiance (Epar, μmols photons m-2 s-1)", title = "Canopy vs Understory Irradiance") +
  theme_bw() +
  scale_color_manual(values = c("Epar_canopy" = "goldenrod", "Epar_under" = "deeppink4"))
ggsave("epar_canopy_under.png", path = "midway_2024/plots/", 
       width = 13, height = 6, units = "in", dpi = 300, scale = 1)

#Get the means for the maximum irradiance by day and plant part
mean_max_epar <- epar_long %>%
  group_by(Date, plant_part) %>%
  summarise(max_Epar = max(irradiance, na.rm = TRUE)) %>%  # get daily max per part
  ungroup() %>%
  group_by(plant_part) %>%
  summarise(mean_max_Epar = mean(max_Epar, na.rm = TRUE))  # get mean of those max values
mean_max_epar
