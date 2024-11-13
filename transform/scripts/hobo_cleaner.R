#Clean Hobo data files by getting rid of unneeded columns, cleaning headers,
#and changing datetime column to two columns for date and time
#By Angela Richards Donà
#October 31, 2024 SPOOKY DAY!
#Updated November 5, 2024 - another spooky day

#load libraries
library(tidyverse)
library(lubridate)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(hms)

#Automate this for all files
hobo_files <- list.files(path = "../data/midway_2024/input/hobo/", pattern = "*.csv", full.names = TRUE)

#read all files and assign new columns based on the bin and run numbers from the file names
hobo_all <- map_dfr(hobo_files, ~ {
  data <- read_csv(.x, show_col_types = FALSE) # Read the CSV file
  bin <- str_extract(basename(.x), "bin\\d+") #extract the bin number
  bin <- as.numeric(str_remove(bin, "bin")) #remove the term bin
  run <- str_extract(basename(.x), "run\\d+") #extract the run number
  run <- as.numeric(str_remove(run, "run")) #remove the term run
  data %>%
mutate(bin = bin, run = run) #add new columns to new concatenated hobo data set       
})      

#Now clean up the data
hobo_all_clean <- hobo_all %>%
  drop_na(`Temperature   (°C)`) %>%
  select(starts_with("Date-Time"), starts_with("Temperature"), "bin", "run") %>% #keep these and get rid of junk columns
  rename(date_time = "Date-Time (U.S. Outlying Islands Standard Time)", #clean these dirty column names
         temp = "Temperature   (°C)") %>%
  separate(col = date_time, #separate datetime column into two columns
           into = c("date", "time"), #name the new columns date and time
           sep = " ", remove = FALSE) %>% #date and time are separated by a space
  mutate(date = mdy(date), #lubridate them! This had to be last to avoid other issues
         time = as_hms(time), #if converted first with lubridate then separated, 00:00:00 was not seen as midnight and was assigned NA
         temp = round(temp, 2)) %>% #round the decimal point to 2 digits
  mutate(plant_part = ifelse(bin %in% 1:6, "canopy", "understory"))
glimpse(hobo_all_clean)

#Get summaries for the data
temp_summary <- hobo_all_clean %>%
  group_by(plant_part, bin, run, date) %>% #summary based on whether canopy or understory, which bin, and run
  summarise(max_temp = max(temp), #get max temps
            min_temp = min(temp), #get min temps
            mean_temp = mean(temp)) #get mean temps
print(temp_summary)

maxtemp_summary <- temp_summary %>%
  ggplot(aes(x = as.factor(bin), y = max_temp, color = as.factor(run))) + #plot it
  geom_boxplot(alpha = 0.75) +
  geom_point(position = position_dodge(width = 0.5), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~ plant_part, scales = "free_x", nrow = 2) +
  labs(x= "Bin",
       y = "Max temperature (°C)",
       title = "Daily Maximum Bin Temperatures by Run",
       color = "Run") +
  scale_color_manual(values = c("darkslateblue", "purple2", "darkslategray3"))
plot(maxtemp_summary)
ggsave("maxtemp_summary.png", path = "plots/")

#plot mean similarly
mean_temp <- temp_summary %>%
  ggplot(aes(x = as.factor(bin), y = mean_temp, color = as.factor(run))) + #plot it
  geom_boxplot(alpha = 0.75) +
  geom_point(position = position_dodge(width = 0.5), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~ plant_part, scales = "free_x", nrow = 2) +
  labs(x= "Bin",
       y = "Mean temperature (°C)",
       title = "Daily Mean Bin Temperatures by Run",
       color = "Run") +
  scale_color_manual(values = c("cornflowerblue", "firebrick3", "cyan3"))
plot(mean_temp)
ggsave("mean_temp_bin_run.png", path = "plots/")

#save dataframe of summary for analysis
summary_temps_midway <- write.csv(temp_summary, file = "../data/midway_2024/transformed/temp_meanmaxmin.csv")
