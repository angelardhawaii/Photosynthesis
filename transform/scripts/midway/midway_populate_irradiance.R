# These are time series of the irradiance measurements for Midway Runs 2 and 3
# The 3rd column is the value of interest (irradiance)
irradiance <- read.csv("../data/midway_2024/input/licor/E_run2_midway.csv")

# Make sure the datread.csv()# Make sure the date is loaded as date
irradiance$posix_date <- as.POSIXct(irradiance$date, format = "%m/%d/%y", tz = "")
irradiance$date_time <- as.POSIXct(paste(irradiance$date, irradiance$time, sep = " "), format = "%m/%d/%y %H:%M:%S")