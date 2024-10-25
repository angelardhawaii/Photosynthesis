# Calculate mean supersaturation time by period 

# For runs 1-7
# Period 1 is days 1-3 and uses Ek of day 1 as threshold
# Period 2 is days 4-5 and uses Ek of day 5 as threshold
# Period 3 is days 6-8 and uses Ek of day 9 as threshold

# For runs 8-9 which don't have a day 5 measurement
# Period 1 is days 1-4 and uses Ek of day 1 as threshold
# Period 2 is days 5-8 and uses Ek of day 9 as threshold

# Desired output columns
# specimen_id run lanai_side treatment temp rlc_order plant_id  species supersat_time_p1 supersat_time_p2 supersat_time_p3 supersat_time_total

# execute the R program that populates the ek data frame
source("./populate_ek.R")

# execute the R program that populates the irradiance data frame
source("./populate_irradiance.R")

# Consider the first day of each run 
supersat <- ek[ek$RLC.Day == 1, c(
  "Specimen.ID", "Run", "Lanai.Side", "Treatment", "Temp...C.",
  "RLC.Order", "Plant.ID", "Species", "rETRmaxYpoint1", "rETRmaxYpoint1_min", "ek.1", "posix_date", "rlc_end_time", "pmax", "pmax_min")]
names(supersat)[13] <- "day1_rlc_time" 
names(supersat)[14] <- "day1_pmax"
names(supersat)[15] <- "day1_pmax_min"

# Populate the column with the hour when the RLC was done, and Ek, on day 9
supersat9 <- ek[ek$RLC.Day == 9, c("Specimen.ID", "Run", "Lanai.Side", "rlc_end_time", "ek.1", "rETRmaxYpoint1_min", "pmax")]
r <- mapply(function(specimen_id, run, lanai_side) {
  d9 <- supersat9[supersat9$Specimen.ID == specimen_id 
                   & supersat9$Run == run 
                   & supersat9$Lanai.Side == lanai_side, ]
  return(c(d9$rlc_end_time, d9$ek.1, d9$rETRmaxYpoint1_min, d9$pmax))
}, supersat$Specimen.ID, supersat$Run, supersat$Lanai.Side)
rm(supersat9)
supersat$day9_rlc_time <- unlist(r[1, ])
supersat$day9_ek <- as.numeric(unlist(r[2, ]))
supersat$day9_rETRmaxYpoint1_min <- as.numeric(unlist(r[3, ]))
supersat$day9_pmax <- as.numeric(unlist(r[4,]))
supersat$day9_pmax_min <- supersat$day9_pmax * 60

# Add new column for Ek on day 5
supersat$day5_ek <- NA
supersat$day5_rETRmaxYpoint1_min <- NA
supersat$day5_pmax <- NA
supersat5 <- ek[ek$RLC.Day == 5, c("Specimen.ID", "Run", "Lanai.Side", "ek.1", "rETRmaxYpoint1_min", "pmax")]
for (i in 1:nrow(supersat5)) {
  supersat[supersat$Specimen.ID == supersat5[i, "Specimen.ID"]
           & supersat$Run == supersat5[i, "Run"]
           & supersat$Lanai.Side == supersat5[i, "Lanai.Side"], ]$day5_ek <- supersat5[i, "ek.1"]
  supersat[supersat$Specimen.ID == supersat5[i, "Specimen.ID"]
           & supersat$Run == supersat5[i, "Run"]
           & supersat$Lanai.Side == supersat5[i, "Lanai.Side"], ]$day5_rETRmaxYpoint1_min <- supersat5[i, "rETRmaxYpoint1_min"]
  supersat[supersat$Specimen.ID == supersat5[i, "Specimen.ID"]
           & supersat$Run == supersat5[i, "Run"]
           & supersat$Lanai.Side == supersat5[i, "Lanai.Side"], ]$day5_pmax <- supersat5[i, "pmax"]
}

rm(supersat5)
# For runs where day5 is missing, use the mean of day 1 and day 9 rETRmaxYpoint1 
find_day5_values <- function(day1_rETRmaxYpoint1, day5_rETRmaxYpoint1, day9_rETRmaxYpoint1, day1_pmax, day5_pmax, day9_pmax) {
  if (is.na(day5_rETRmaxYpoint1)) {
    return(c(mean(c(day1_rETRmaxYpoint1, day9_rETRmaxYpoint1)), mean(c(day1_pmax, day9_pmax))))
  } 
  return(c(day5_rETRmaxYpoint1, day5_pmax))
}
r <- mapply(find_day5_values, 
            supersat$rETRmaxYpoint1_min, 
            supersat$day5_rETRmaxYpoint1_min, 
            supersat$day9_rETRmaxYpoint1_min, 
            supersat$day1_pmax, 
            supersat$day5_pmax, 
            supersat$day9_pmax)
supersat$day5_rETRmaxYpoint1_min <- unlist(r[1, ])
supersat$day5_pmax <- unlist(r[2, ])
rm(r)
supersat$day5_pmax_min <- supersat$day5_pmax * 60


calc_supersat_by_day <- function(run_irradiance, day1_date, days_to_consider, threshold, lanai_side) {
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
                                        & run_irradiance$Lanai.Side == lanai_side
                                        & run_irradiance$Epar > threshold, ])
    daylight_minutes[n] = nrow(run_irradiance[run_irradiance$date_time > start 
                                           & run_irradiance$date_time < end 
                                           & run_irradiance$Lanai.Side == lanai_side
                                           & run_irradiance$Epar >= 1, ])
    rel_supersat_by_day[n] = supersat_by_day[n] / daylight_minutes[n]
  }
  return(list(supersat_by_day=supersat_by_day, daylight_minutes=daylight_minutes, rel_supersat_by_day=rel_supersat_by_day))
}

# Test the function above
#x = calc_supersat_by_day(irradiance, supersat[1, "posix_date"], 0:2, supersat[1, "ek.1"], "ewa")

calculate_supersat <- function(day1_date, day1_rlc_time, day9_rlc_time, day1_ek, 
                               day5_ek, day9_ek, run, lanai_side) {
  run_start <- as.POSIXct(paste(day1_date, day1_rlc_time, sep = " "))
  run_end <- as.POSIXct(paste(day1_date + 86400 * 8, day9_rlc_time, sep = " "))
  run_irradiance <- irradiance[irradiance$date_time > run_start & irradiance$date_time < run_end, ]
  if (run == 8 | run == 9) {
    r = calc_supersat_by_day(run_irradiance, day1_date, 0:2, day1_ek, lanai_side)
    supersat_p1 <- r$supersat_by_day
    day_length_p1 <- r$daylight_minutes
    supersat_rel_p1 <- r$rel_supersat_by_day
    r = calc_supersat_by_day(run_irradiance, day1_date, 3:4, day1_ek, lanai_side)
    supersat_p2 <- r$supersat_by_day
    day_length_p2 <- r$daylight_minutes
    supersat_rel_p2 <- r$rel_supersat_by_day
    r = calc_supersat_by_day(run_irradiance, day1_date, 5:8, day9_ek, lanai_side)
    supersat_p3 <- r$supersat_by_day
    day_length_p3 <- r$daylight_minutes
    supersat_rel_p3 <- r$rel_supersat_by_day
  } else {
    r = calc_supersat_by_day(run_irradiance, day1_date, 0:2, day1_ek, lanai_side)
    supersat_p1 <- r$supersat_by_day
    day_length_p1 <- r$daylight_minutes
    supersat_rel_p1 <- r$rel_supersat_by_day
    r = calc_supersat_by_day(run_irradiance, day1_date, 3:4, day5_ek, lanai_side)
    supersat_p2 <- r$supersat_by_day
    day_length_p2 <- r$daylight_minutes
    supersat_rel_p2 <- r$rel_supersat_by_day
    r = calc_supersat_by_day(run_irradiance, day1_date, 5:8, day9_ek, lanai_side)
    supersat_p3 <- r$supersat_by_day
    day_length_p3 <- r$daylight_minutes
    supersat_rel_p3 <- r$rel_supersat_by_day
  }
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
            supersat$ek.1, 
            supersat$day5_ek,
            supersat$day9_ek,
            supersat$Run,
            supersat$Lanai.Side)
supersat$supersat_p1 <- unlist(r[1,])
supersat$supersat_p2 <- unlist(r[2,])
supersat$supersat_p3 <- unlist(r[3,])
supersat$supersat_avg <- unlist(r[4,])
supersat$day_length_avg <- unlist(r[5,])
supersat$supersat_rel <- round(unlist(r[6,]), 3) * 100
rm(r)

supersat$carbon1 <- supersat$supersat_p1 * supersat$rETRmaxYpoint1_min
supersat$carbon5 <- supersat$supersat_p2 * supersat$day5_rETRmaxYpoint1_min
supersat$carbon9 <- supersat$supersat_p3 * supersat$day9_rETRmaxYpoint1_min
supersat$carbon_total <- supersat$carbon1 + supersat$carbon5 + supersat$carbon9

supersat$dspi_pmax1 <- supersat$supersat_p1 * supersat$day1_pmax_min
supersat$dspi_pmax5 <- supersat$supersat_p2 * supersat$day5_pmax_min
supersat$dspi_pmax9 <- supersat$supersat_p3 * supersat$day9_pmax_min
supersat$dspi_pmax_total <- supersat$dspi_pmax1 + supersat$dspi_pmax5 + supersat$dspi_pmax9

write.csv(supersat, "./output/mean_supersaturation_by_period.csv")
write.csv(supersat, "../ek_irrad_model/input_data/mean_supersaturation_by_period.csv")
