#Clean Hobo data files by getting rid of unneeded columns, cleaning headers,
#and changing datetime column to two columns for date and time
#Summary dataset created and plots made from this
#Total time temperatures exceeded 28 C extracted for further analysis by bin
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
library(zoo) #for rolling means

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
  mutate(plant_part = ifelse(bin %in% 1:6, "canopy", "under"))
glimpse(hobo_all_clean)

#Get summaries for the data by date
temp_summary <- hobo_all_clean %>%
  group_by(plant_part, bin, run, date) %>% #summary based on whether canopy or understory, which bin, and run
  summarise(max_temp = max(temp), #get max temps
            min_temp = min(temp), #get min temps
            mean_temp = mean(temp)) #get mean temps
print(temp_summary)

#Get summaries for the data based on time
temp_summary_time <- hobo_all_clean %>%
  group_by(plant_part, run, date, time) %>% #summary based on whether canopy or understory, which bin, and run
  summarise(max_temp = max(temp), #get max temps
            min_temp = min(temp), #get min temps
            mean_temp = mean(temp)) #get mean temps
print(temp_summary_time)


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

#Extract time in minutes that temperature exceeded 28C
run1 <- as.Date(c("2024-07-25", "2024-07-26", "2024-07-27", "2024-07-28", "2024-07-29", 
                  "2024-07-30", "2024-07-31", "2024-08-01", "2024-08-02")) #make this vecture so to exclude run 1
time_over_28 <- hobo_all_clean %>%
  filter(temp >= 28, !date %in% run1 ) %>% #filter rows where temperature is above 28 and keep only those from runs 2 and 3
  group_by(plant_part, run, bin, date) %>%
  summarise(daily_min_over_28 = n() * 10) #count rows and multiply by 10 minutes per instance
glimpse(time_over_28)

#plot it
hot_plot <- time_over_28 %>%
  ggplot(aes(x = as.factor(bin), y = daily_min_over_28, color = as.factor(run))) +
  geom_boxplot(alpha = 0.75) +
  geom_point(position = position_dodge(width = 0.5), alpha = 0.5) +
  theme_bw() + 
  facet_wrap(~ plant_part + bin, scales = "free_x", nrow = 2) +
  labs(x = "Bin",
       y = "Minutes over 28°C",
       color = "Run") +
  scale_color_manual(values = c("orange", "firebrick3"))
plot(hot_plot)
ggsave("min_over_28.png", path = "plots/")
write.csv(time_over_28, file = "../data/midway_2024/transformed/time_over_28.csv")

#Take sum of time over 28 for all bins
time_over_28_sum <- time_over_28 %>%
  group_by(plant_part, bin, run) %>%
  summarise(sum_over_28 = round(sum(daily_min_over_28), digits = 1))
glimpse(time_over_28_sum)
write.csv(time_over_28_sum, file = "../data/midway_2024/transformed/time_over_28_sum.csv")

#Get mean of all bins in canopy and understory sections
mins_over_28_plant_part <- time_over_28_sum %>%
  group_by(plant_part, run) %>%
  summarise(mean_mins28_plant_part = round(mean(sum_over_28), digits = 0))
glimpse(mins_over_28_plant_part)
write.csv(mins_over_28_plant_part, file = "../data/midway_2024/transformed/mins_28_plant_part.csv")

#Meh, plot daily raw temperatures by bin
daily_raw_bin2 <- temp_summary %>%
  filter(run == 2) %>%
  ggplot(aes(x = date,
             y = mean_temp,
             color = as.factor(run))) +
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~plant_part + bin, nrow = 4) +
  labs(x = "Date",
       y = "Mean Temperature (°C) by Day",
       color = "Run") +
  scale_color_manual(values = c("darkslateblue", "cyan3")) +
  theme_bw()
  
plot(daily_raw_bin2)




