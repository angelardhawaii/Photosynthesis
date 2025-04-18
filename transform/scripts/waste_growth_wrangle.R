#Make minor changes and additions to the growth dataset that has already been cleaned for use in Google Sheets
#By Angela Richards Dona
#April 4, 2025

#Load libraries
library(lubridate)
library(hms)
library(dplyr)
library(janitor)
library(lunar)
library(suncalc)

#open weight dataset
waste_p2_growth <- read.csv("../data/wastewater/input/growth_2025/p2_growth_all.csv")

#clean those column titles
clean_names(waste_p2_growth)


#make sure date is properly formatted
waste_p2_growth <- waste_p2_growth %>%
  mutate(date = mdy(date)) 

#assigns treatment as characters from integers then to factors
waste_p2_growth$treatment <- as.factor(as.character(waste_p2_growth$treatment))
waste_p2_growth$salinity <- as.factor(as.character(waste_p2_growth$salinity))

#add column for lunar phase to be used for Ulva only
waste_p2_growth <- waste_p2_growth %>%
  mutate(lunar_phase = lunar.phase(date, name = TRUE))

#for more precision on lunar illumination
moon_data <- getMoonIllumination(waste_p2_growth$date)

waste_p2_growth <- waste_p2_growth %>%
  mutate(illumination = moon_data$fraction,
         lunar_illumination = case_when(
           illumination < 0.1 ~ "new",
           illumination >= 0.1 & illumination < 0.4 ~ "crescent",
           illumination >= 0.4 & illumination < 0.6 ~ "first_q / last_q",
           illumination >= 0.6 & illumination < 0.9 ~ "gibbous",
           illumination >= 0.9 ~ "full"
         ))
write.csv(waste_p2_growth, file = "../data/wastewater/input/growth_2025/waste_p2_growth.csv", 
          row.names = FALSE)
