# Calculate mean supersaturation time by period 

# Period 1 is days 1-3 and uses Ek of day 1 as threshold
# Period 2 is days 4-5 and uses Ek of day 5 as threshold
# Period 3 is days 6-8 and uses Ek of day 9 as threshold

# Desired output columns
# specimen_id run lanai_side treatment temp rlc_order plant_id  species supersat_time_p1 supersat_time_p2 supersat_time_p3 supersat_time_total

# execute the R program that populates the ek data frame
source("acan_populate_ek.R")

# execute the R program that populates the irradiance data frame
source("acan_populate_irradiance.R")

# Make a new dataframe by selecting the named columns from the first day of each run in the dataset ek 
supersat <- ek[ek$rlc_day == 1, c(
  "specimenID", "run", "treatment", "temperature",
  "rlc_order", "plantID", "ek.est", "posix_date", "rlc_end_time", "pmax", "pmax_min")]

#rename the last three columns
names(supersat)[9] <- "day1_rlc_time" 
names(supersat)[10] <- "day1_pmax"
names(supersat)[11] <- "day1_pmax_min"

# Make another new dataframe as above and populate a new column with the hour when the RLC was done, and Ek, on day 9
supersat9 <- ek[ek$rlc_day == 9, c("specimenID", "run", "rlc_end_time", "ek.est", "pmax")]
r <- mapply(function(specimenID, run) {
  d9 <- supersat9[supersat9$specimenID == specimenID
                  & supersat9$run == run, ]
  return(c(d9$rlc_end_time, d9$ek.est, d9$pmax))
}, supersat$specimenID, supersat$run)
rm(supersat9)
supersat$day9_rlc_time <- unlist(r[1, ])
supersat$day9_ek <- as.numeric(unlist(r[2, ]))
supersat$day9_pmax <- as.numeric(unlist(r[3,]))
supersat$day9_pmax_min <- supersat$day9_pmax * 60

# Add new column for Ek on day 5
supersat$day5_ek <- NA
supersat$day5_pmax <- NA
supersat$day5_rlc_end_time <- NA
supersat5 <- ek[ek$rlc_day == 5, c("specimenID", "run", "rlc_end_time", "ek.est",  "pmax")]
for (i in 1:nrow(supersat5)) {
  supersat[supersat$specimenID == supersat5[i, "specimenID"]
           & supersat$run == supersat5[i, "run"], ]$day5_ek <- supersat5[i, "ek.est"]
  supersat[supersat$specimenID == supersat5[i, "specimenID"]
           & supersat$run == supersat5[i, "run"], ]$day5_rlc_end_time <- supersat5[i, "rlc_end_time"]
  supersat[supersat$specimenID == supersat5[i, "specimenID"]
           & supersat$run == supersat5[i, "run"], ]$day5_pmax <- supersat5[i, "pmax"]
}
supersat$day5_pmax_min <- supersat$day5_pmax * 60
rm(supersat5)

calc_supersat_by_day <- function(run_irradiance, day1_date, days_to_consider, threshold) {
  supersat_by_day <- rep(NA, length(days_to_consider))
  daylight_minutes <- rep(NA, length(days_to_consider))
  rel_supersat_by_day <- rep(NA, length(days_to_consider))
  n = 0
  for (day in days_to_consider) {
    start <- as.POSIXct(paste(day1_date + 86400 * day, "00:00:01", sep = " "))
    end <- as.POSIXct(paste(day1_date + 86400 * day, "23:59:59", sep = " "))
    n = n + 1
    supersat_by_day[n] = nrow(run_irradiance[run_irradiance$date_time > start 
                                             & run_irradiance$date_time < end 
                                             & run_irradiance$Epar > threshold, ])
    daylight_minutes[n] = nrow(run_irradiance[run_irradiance$date_time > start 
                                              & run_irradiance$date_time < end 
                                              & run_irradiance$Epar >= 1, ])
    rel_supersat_by_day[n] = supersat_by_day[n] / daylight_minutes[n]
  }
  return(list(supersat_by_day=supersat_by_day, daylight_minutes=daylight_minutes, rel_supersat_by_day=rel_supersat_by_day))
}

# Test the function above
x = calc_supersat_by_day(irradiance, supersat[1, "posix_date"], 0:2, supersat[1, "ek.est"])

calculate_supersat <- function(day1_date, day1_rlc_time, day9_rlc_time, day1_ek, 
                               day5_ek, day9_ek, run) {
  run_start <- as.POSIXct(paste(day1_date, day1_rlc_time, sep = " "))
  run_end <- as.POSIXct(paste(day1_date + 86400 * 8, day9_rlc_time, sep = " "))
  run_irradiance <- irradiance[irradiance$date_time > run_start & irradiance$date_time < run_end, ]

    r = calc_supersat_by_day(run_irradiance, day1_date, 0:2, day1_ek)
    supersat_p1 <- r$supersat_by_day
    day_length_p1 <- r$daylight_minutes
    supersat_rel_p1 <- r$rel_supersat_by_day
    r = calc_supersat_by_day(run_irradiance, day1_date, 3:4, day5_ek)
    supersat_p2 <- r$supersat_by_day
    day_length_p2 <- r$daylight_minutes
    supersat_rel_p2 <- r$rel_supersat_by_day
    r = calc_supersat_by_day(run_irradiance, day1_date, 5:8, day9_ek)
    supersat_p3 <- r$supersat_by_day
    day_length_p3 <- r$daylight_minutes
    supersat_rel_p3 <- r$rel_supersat_by_day

  return(list(mean(supersat_p1), mean(supersat_p2), mean(supersat_p3), 
              mean(c(supersat_p1, supersat_p2, supersat_p3)), 
              mean(c(day_length_p1, day_length_p2, day_length_p3), na.rm = TRUE),
              mean(c(supersat_rel_p1, supersat_rel_p2, supersat_rel_p3, na.rm = TRUE))
  ))
}

r <- mapply(calculate_supersat, 
            supersat$posix_date,
            supersat$day1_rlc_time, 
            supersat$day9_rlc_time, 
            supersat$ek.est, 
            supersat$day5_ek,
            supersat$day9_ek,
            supersat$run)
supersat$supersat_p1 <- unlist(r[1,])
supersat$supersat_p2 <- unlist(r[2,])
supersat$supersat_p3 <- unlist(r[3,])
supersat$supersat_avg <- unlist(r[4,])
supersat$day_length_avg <- unlist(r[5,])
supersat$supersat_rel <- round(unlist(r[6,]), 3) * 100
rm(r)

supersat$dspi_pmax1 <- supersat$supersat_p1 * supersat$day1_pmax_min
supersat$dspi_pmax5 <- supersat$supersat_p2 * supersat$day5_pmax_min
supersat$dspi_pmax9 <- supersat$supersat_p3 * supersat$day9_pmax_min
supersat$dspi_pmax_total <- supersat$dspi_pmax1 + supersat$dspi_pmax5 + supersat$dspi_pmax9

write.csv(supersat, "../../../../data/limu/acan_2024/transformed/acan_hsat_dspi.csv")

