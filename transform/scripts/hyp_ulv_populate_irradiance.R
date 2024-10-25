# These are time series of the irradiance measurements
# The 3rd column is the value of interest (irradiance)
ewa_irradiance_files <- dir("./data_irradiance/ewa", full.names = TRUE)
ewa_irradiance <- do.call(rbind, lapply(ewa_irradiance_files, read.csv))
ewa_irradiance$Lanai.Side <- as.factor("ewa")

dia_irradiance_files <- dir("./data_irradiance/dia", full.names = TRUE)
dia_irradiance <- do.call(rbind, lapply(dia_irradiance_files, read.csv))
dia_irradiance$Lanai.Side <- as.factor("diamond")

irradiance <- rbind(ewa_irradiance, dia_irradiance)
rm(dia_irradiance, ewa_irradiance, dia_irradiance_files, ewa_irradiance_files)

# Loads all the files in the same data frame

# Make sure the date is loaded as date
irradiance$posix_date <- as.POSIXct(irradiance$Date, format = "%m/%d/%y", tz = "")
irradiance$date_time <- as.POSIXct(paste(irradiance$Date, irradiance$Time, sep = " "), format = "%m/%d/%y %H:%M:%S")