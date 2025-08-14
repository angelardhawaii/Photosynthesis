# Calculate mean supersaturation time by period 

# Period 1 is days 1-3 and uses Ek of day 1 as threshold
# Period 2 is days 4-5 and uses Ek of day 5 as threshold
# Period 3 is days 6-8 and uses Ek of day 9 as threshold

# Desired output columns
# specimenID run nitrate rlc_order plantID  plant_part supersat_time_p1 supersat_time_p2 supersat_time_p3 supersat_time_total

#By Angela Richards Donà
#Modified script for Midway data September-November 2024


# execute the R script that populates the ek data frame
source("midway_populate_ek.R")

# execute the R program that populates the irradiance data frame
source("midway_populate_irradiance.R")

# Make a new dataframe by selecting the named columns from the first day of each run in the dataset ek 
hsat <- ek[ek$rlc_day == 1, c(
  "specimenID", "run", "nitrate", "salinity",
  "rlc_order", "plantID", "ek.est", "date", "rlc_end_time", "pmax", "pmax_min")]

#rename the last three columns
names(hsat)[9] <- "day1_rlc_time" 
names(hsat)[10] <- "day1_pmax"
names(hsat)[11] <- "day1_pmax_min"

# Make another new dataframe as above and populate a new column with the hour when the RLC was done, and Ek, on day 9
hsat9 <- ek[ek$rlc_day == 9, c("specimenID", "run", "rlc_end_time", "ek.est", "pmax")]
r <- mapply(function(specimenID, run) {
  d9 <- hsat9[hsat9$specimenID == specimenID
                  & hsat9$run == run, ]
  return(c(d9$rlc_end_time, d9$ek.est, d9$pmax))
}, hsat$specimenID, hsat$run)
rm(hsat9)
hsat$day9_rlc_time <- unlist(r[1, ])
hsat$day9_ek <- as.numeric(unlist(r[2, ]))
hsat$day9_pmax <- as.numeric(unlist(r[3,]))
hsat$day9_pmax_min <- hsat$day9_pmax * 60

# Add new column for Ek on day 5
hsat$day5_ek <- NA
hsat$day5_pmax <- NA
hsat$day5_rlc_end_time <- NA
hsat5 <- ek[ek$rlc_day == 5, c("specimenID", "run", "rlc_end_time", "ek.est",  "pmax")]
for (i in 1:nrow(hsat5)) {
  hsat[hsat$specimenID == hsat5[i, "specimenID"]
           & hsat$run == hsat5[i, "run"], ]$day5_ek <- hsat5[i, "ek.est"]
  hsat[hsat$specimenID == hsat5[i, "specimenID"]
           & hsat$run == hsat5[i, "run"], ]$day5_rlc_end_time <- hsat5[i, "rlc_end_time"]
  hsat[hsat$specimenID == hsat5[i, "specimenID"]
           & hsat$run == hsat5[i, "run"], ]$day5_pmax <- hsat5[i, "pmax"]
}
hsat$day5_pmax_min <- hsat$day5_pmax * 60
rm(hsat5)

calc_hsat_by_day <- function(run_irradiance, day1_date, days_to_consider, threshold) {
  hsat_by_day <- rep(NA, length(days_to_consider))
  daylight_minutes <- rep(NA, length(days_to_consider))
  rel_hsat_by_day <- rep(NA, length(days_to_consider))
  n = 0
  for (day in days_to_consider) {
    start <- as.POSIXct(paste(day1_date + 86400 * day, "00:00:01", sep = " "))
    end <- as.POSIXct(paste(day1_date + 86400 * day, "23:59:59", sep = " "))
    n = n + 1
    hsat_by_day[n] = nrow(run_irradiance[run_irradiance$date_time > start 
                                             & run_irradiance$date_time < end 
                                             & run_irradiance$Epar > threshold, ])
    daylight_minutes[n] = nrow(run_irradiance[run_irradiance$date_time > start 
                                              & run_irradiance$date_time < end 
                                              & run_irradiance$Epar >= 1, ])
    rel_hsat_by_day[n] = hsat_by_day[n] / daylight_minutes[n]
  }
  return(list(hsat_by_day=hsat_by_day, daylight_minutes=daylight_minutes, rel_hsat_by_day=rel_hsat_by_day))
}

# Test the function above
x = calc_hsat_by_day(irradiance, hsat[1, "posix_date"], 0:2, hsat[1, "ek.est"])

calculate_hsat <- function(day1_date, day1_rlc_time, day9_rlc_time, day1_ek, 
                               day5_ek, day9_ek, run) {
  run_start <- as.POSIXct(paste(day1_date, day1_rlc_time, sep = " "))
  run_end <- as.POSIXct(paste(day1_date + 86400 * 8, day9_rlc_time, sep = " "))
  run_irradiance <- irradiance[irradiance$date_time > run_start & irradiance$date_time < run_end, ]
  
  r = calc_hsat_by_day(run_irradiance, day1_date, 0:2, day1_ek)
  hsat_p1 <- r$hsat_by_day
  day_length_p1 <- r$daylight_minutes
  hsat_rel_p1 <- r$rel_hsat_by_day
  r = calc_hsat_by_day(run_irradiance, day1_date, 3:4, day5_ek)
  hsat_p2 <- r$hsat_by_day
  day_length_p2 <- r$daylight_minutes
  hsat_rel_p2 <- r$rel_hsat_by_day
  r = calc_hsat_by_day(run_irradiance, day1_date, 5:8, day9_ek)
  hsat_p3 <- r$hsat_by_day
  day_length_p3 <- r$daylight_minutes
  hsat_rel_p3 <- r$rel_hsat_by_day
  
  return(list(mean(hsat_p1), mean(hsat_p2), mean(hsat_p3), 
              mean(c(hsat_p1, hsat_p2, hsat_p3)), 
              mean(c(day_length_p1, day_length_p2, day_length_p3), na.rm = TRUE),
              mean(c(hsat_rel_p1, hsat_rel_p2, hsat_rel_p3, na.rm = TRUE))
  ))
}

r <- mapply(calculate_hsat, 
            hsat$posix_date,
            hsat$day1_rlc_time, 
            hsat$day9_rlc_time, 
            hsat$ek.est, 
            hsat$day5_ek,
            hsat$day9_ek,
            hsat$run)
hsat$hsat_p1 <- unlist(r[1,])
hsat$hsat_p2 <- unlist(r[2,])
hsat$hsat_p3 <- unlist(r[3,])
hsat$hsat_avg <- unlist(r[4,])
hsat$day_length_avg <- unlist(r[5,])
hsat$hsat_rel <- round(unlist(r[6,]), 3) * 100
rm(r)

hsat$dspi_pmax1 <- hsat$hsat_p1 * hsat$day1_pmax_min
hsat$dspi_pmax5 <- hsat$hsat_p2 * hsat$day5_pmax_min
hsat$dspi_pmax9 <- hsat$hsat_p3 * hsat$day9_pmax_min
hsat$dspi_pmax_total <- hsat$dspi_pmax1 + hsat$dspi_pmax5 + hsat$dspi_pmax9

write.csv(hsat, "../../../../data/limu/acan_2024/transformed/acan_hsat_dspi.csv")

