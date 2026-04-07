#Clean Hobo data files by getting rid of unneeded columns, cleaning headers,
#and changing datetime column to two columns for date and time
#Summary dataset created and plots made from this
#Total time temperatures exceeded 28 C extracted for further analysis by bin
#By Angela Richards Donà
#October 31, 2024 SPOOKY DAY!
#Updated for wastewater data November 13, 2025

#load libraries
library(tidyverse)
library(lubridate)
library(hms)
library(stringr)
library(glue)

# 1) Directories
input_dir  <- "/Users/angela/src/Photosynthesis/data/wastewater/input/hobo_data"
output_dir <- "/Users/angela/src/Photosynthesis/data/wastewater/input/hobo_data/hobo_run" 
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2) Define run date windows
runs <- tibble(
  run   = c("run1",        "run2",        "run3",        "run4"),
  start = ymd(c("2025-01-23", "2025-02-17", "2025-03-06", "2025-03-17")),
  end   = ymd(c("2025-01-31", "2025-02-25", "2025-03-18", "2025-03-28"))
)

# Helper to turn a date into a run label
which_run <- function(d) {
  case_when(
    d >= ymd("2025-01-23") & d <= ymd("2025-01-31") ~ "run1",
    d >= ymd("2025-02-17") & d <= ymd("2025-02-25") ~ "run2",
    d >= ymd("2025-03-06") & d <= ymd("2025-03-18") ~ "run3",
    d >= ymd("2025-03-17") & d <= ymd("2025-03-28") ~ "run4",
    TRUE ~ NA_character_
  )
}

# 3) Read and standardise all HOBO files into ONE data frame: hobo_all
hobo_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

hobo_all <- map_dfr(hobo_files, function(f) {
  hobo_id <- str_sub(basename(f), 1, 2)   # first two chars = HOBO number
  
  # Read everything as character so readr won't guess wrong types
  df_raw <- read_csv(
    f,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    # ⬇️ KEEP ONLY the two columns you actually need
    select(
      date_time_raw = `Date-Time (HST)`,
      temp_raw      = `Temperature   (°C)`
    )
  
  df_raw %>%
    mutate(
      # Parse datetime with flexible formats (with or without seconds)
      datetime = parse_date_time(
        date_time_raw,
        orders = c("mdy HMS", "mdy HM"),
        quiet  = TRUE
      ),
      date = as_date(datetime),
      time = as_hms(datetime),
      
      # Parse temperature safely
      temp = parse_double(
        temp_raw,
        na = c("", "NA", "NaN", "---", "LOW BAT")
      ),
      
      run  = which_run(date),
      hobo = hobo_id
    ) %>%
    # Remove rows that clearly aren't real data (unparsed datetime)
    filter(!is.na(datetime))
})
#Drop rows that are not in a run
hobo_all <- hobo_all %>%
  filter(!is.na(run))

#Write out one .csv file for each hobo and run combo
hobo_all %>%
  group_by(hobo, run) %>%
  group_walk(~ {
    out_name <- paste0(.y$hobo, "_", .y$run, ".csv")
    out_path <- file.path(output_dir, out_name)
    write_csv(.x, out_path)
    message("Wrote: ", out_path)
  })

hobo_all <- hobo_all %>%
  drop_na(temp) %>%                             # same as drop_na(`Temperature   (°C)`)
  # We already selected/renamed earlier, so no need for select()/rename()/separate()
  mutate(
    temp = round(temp, 2),
    )
glimpse(hobo_all)


#______________
#Plots and Summaries
temp_summary <- hobo_all %>%
  group_by(hobo, run, date) %>%
  summarise(
    max_temp  = max(temp),
    min_temp  = min(temp),
    mean_temp = mean(temp),
    .groups   = "drop"
  )

print(temp_summary)

temp_summary_time <- hobo_all %>%
  group_by(hobo, run, date, time) %>%
  summarise(
    max_temp  = max(temp),
    min_temp  = min(temp),
    mean_temp = mean(temp),
    .groups   = "drop"
  )

print(temp_summary_time)

# Daily max by bin & run
ww_maxtemp_summary <- temp_summary %>%
  ggplot(aes(x = as.factor(hobo), y = max_temp, color = as.factor(run))) +
  geom_boxplot(alpha = 0.75) +
  geom_point(position = position_dodge(width = 0.5), alpha = 0.5) +
  theme_bw() +
  #facet_wrap(~ plant_part, scales = "free_x", nrow = 2) +
  labs(
    x     = "Hobo",
    y     = "Max temperature (°C)",
    title = "Daily Maximum Bin Temperatures by Run",
    color = "Run"
  ) +
  scale_color_manual(values = c("darkslateblue", "purple2", "darkslategray3", "firebrick4"))

plot(ww_maxtemp_summary)
ggsave("ww_maxtemp_summary.png", path = "plots/")

# Daily mean by bin & run
ww_mean_temp_plot <- temp_summary %>%
  ggplot(aes(x = as.factor(hobo), y = mean_temp, color = as.factor(run))) +
  geom_boxplot(alpha = 0.75) +
  geom_point(position = position_dodge(width = 0.5), alpha = 0.5) +
  theme_bw() +
  #facet_wrap(~ plant_part, scales = "free_x", nrow = 2) +
  labs(
    x     = "Hobo",
    y     = "Mean temperature (°C)",
    title = "Daily Mean Bin Temperatures by Run",
    color = "Run"
  ) +
  scale_color_manual(values = c("cornflowerblue", "firebrick3", "cyan3", "orange3"))

plot(ww_mean_temp_plot)
ggsave("ww_mean_temp_bin_run.png", path = "plots/")
