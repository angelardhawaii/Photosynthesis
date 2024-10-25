#This script takes the cleaned data from PAM output and runs the Phytotools package to
#produce a final dataset that can be used for analysis.
#take cleaned data from Python script and run through phytotools package to get ek and alpha etc
#open appropriate libraries
library("phytotools")
library("hash")
library("dplyr")
library("hms")
library("lubridate")

mid_ptools <- read.csv("../../../../data/limu/midway/transformed/midway_pam_clean.csv", sep = ",")

# add a column that turns date format into POSIXct
mid_ptools$posix_date <- as.POSIXct(mid_ptools$Date, format = "%Y-%m-%d")

#same for time
mid_ptools$hms_time <- as_hms(mid_ptools$Time)

#add a new column that gets rid of characters in ID column
mid_ptools$IDnumber <- as.numeric(substr(mid_ptools$ID, 3, 5))

#get only last two digits from IDnumber for plantID to correspond to individuals
mid_ptools$plantID <- mid_ptools$IDnumber

#add a new column that assigns treatment based on the remainder of integer division by 8 (modulus)
mid_ptools$treatment <- ifelse(mid_ptools$plantID %% 8 == 0, "4",(((mid_ptools$plantID %% 8) +1) %/% 2))

#add a new column for salinity based on plantID
mid_ptools$salinity <- ifelse(mid_ptools$plantID %% 2 == 0, "28 ppt", "35 ppt")

#add a column for canopy or understory
mid_ptools$plant_part <- ifelse(between(mid_ptools$plantID, 01, 48), "canopy", "under")

#add a column for run
mid_ptools$run <- as.character(ifelse(mid_ptools$Date == "2024-08-03" | mid_ptools$Date == "2024-08-07" | mid_ptools$Date == "2024-08-11", 1, 2))

#add a new column that bins the measurement periods by date and AM/PM
bin_times <- list(
        bin1 = c(hms("07:00:00"), hms("12:00:00")),
        bin2 = c(hms("13:00:00"), hms("18:00:00"))        
)

mid_ptools <- mid_ptools %>%
        mutate(day_group = case_when(
                posix_date == as.Date("2024-08-03") & hms_time > bin_times$bin1[1] & hms_time < bin_times$bin1[2] ~ "A",
                posix_date == as.Date("2024-08-03") & hms_time > bin_times$bin2[1] & hms_time < bin_times$bin2[2] ~ "B",
                posix_date == as.Date("2024-08-07") & hms_time > bin_times$bin1[1] & hms_time < bin_times$bin1[2] ~ "C",
                posix_date == as.Date("2024-08-07") & hms_time > bin_times$bin2[1] & hms_time < bin_times$bin2[2] ~ "D",
                posix_date == as.Date("2024-08-11") & hms_time > bin_times$bin1[1] & hms_time < bin_times$bin1[2] ~ "E",
                posix_date == as.Date("2024-08-11") & hms_time > bin_times$bin2[1] & hms_time < bin_times$bin2[2] ~ "F",
                posix_date == as.Date("2024-08-12") & hms_time > bin_times$bin1[1] & hms_time < bin_times$bin1[2] ~ "G",
                posix_date == as.Date("2024-08-12") & hms_time > bin_times$bin2[1] & hms_time < bin_times$bin2[2] ~ "H",
                posix_date == as.Date("2024-08-16") & hms_time > bin_times$bin1[1] & hms_time < bin_times$bin1[2] ~ "I",
                posix_date == as.Date("2024-08-16") & hms_time > bin_times$bin2[1] & hms_time < bin_times$bin2[2] ~ "J",
                posix_date == as.Date("2024-08-20") & hms_time > bin_times$bin1[1] & hms_time < bin_times$bin1[2] ~ "K",
                posix_date == as.Date("2024-08-20") & hms_time > bin_times$bin2[1] & hms_time < bin_times$bin2[2] ~ "L"
))

##RLC-order is added to dataset in model script

# add a new column with date and specimen ID as unique key
mid_ptools$uid <- paste(mid_ptools$posix_date, mid_ptools$ID, sep = "_")

mid_ptools$NPQ <- as.numeric(mid_ptools$NPQ)

# create a new column based on Y(II) at Epar 0 (Effective quantum yield)
mid_ptools <- transform(mid_ptools, QuanYield = ifelse(Epar == "0", Y.II., NA))

# remove all rows where rETR is null or negative
selectedData <- subset(mid_ptools, rETR > 0) 

# unique function eliminates duplicates returns all unique IDs into a vector
uniqueIds <- unique(selectedData$uid)

# store the number of unique IDs in the variable n
n <- length(uniqueIds)

# create a matrix full of NAs with n rows
rETRMaxes = array(NA,c(n,1))
rlc_end_times = array(NA, c(n,1))
rETRmaxYpoint1 = array(NA, c(n,1))
maxNPQ_ypoint1 = array(NA, c(n,1))
deltaNPQ_ypoint1 = array(NA, c(n,1))
minNPQ = array(NA, c(n,1))


# create the rETRmax column
#also create a rETRmaxYpoint1 column per Beer & Axelsson where we take the max value but only if Y>0.1
for (i in 1:n){
        light_curve <- subset(selectedData, uid == uniqueIds[i])
        subMaxETR <- max(light_curve$rETR)
        rETRmaxYpoint1[i] <- max(subset(light_curve, Y.II. > .1)$rETR, na.rm = TRUE)
        maxNPQ_ypoint1[i] <- as.numeric(max(subset(light_curve, Y.II. > .1)$NPQ, na.rm = TRUE))
        minNPQ[i] <- as.numeric(max(subset(light_curve, Epar == 65)$NPQ))
        #store the subMaxETR in the new column but only in rows where uid is same as uniqueIds of i
        selectedData$rETRmax[selectedData$uid == uniqueIds[i]] <- subMaxETR
        #also store subMaxETR in the matrix rETRMaxes created previously, to later calculate ETRmax
        rETRMaxes[i] = subMaxETR
        #temperatures[i] = max(sub$Temp)
        rlc_end_times[i] <- max(light_curve$Time)
        
}

# prepare empty matrices to hold output from fitWebb
alpha <- array(NA,c(n,4))
colnames(alpha) <- c("est", "st_err", "t", "p")
ek <- array(NA,c(n,4))
colnames(ek) <- c("est", "st_err", "t", "p")

#create a treatmemt vector from SelecteData treatment column
treatment_names = c("0.5uM", "2uM", "4uM", "8uM")
treatment = array(NA,c(n,1))



for (i in 1:n){
        #Get ith data
        Epar <- selectedData$Epar[selectedData$uid==uniqueIds[i]]
        rETR <- selectedData$rETR[selectedData$uid==uniqueIds[i]]
        y.II <- selectedData$Y.II.[selectedData$uid==uniqueIds[i]]
        
        #Call function (using y.II and adding normalize = TRUE will yield more accurate PE parameters)
        myfit <- fitWebb(Epar, y.II, normalize = TRUE)
        #store the four values outputs in the matrix
        alpha[i,] <- myfit$alpha
        ek[i,] <- myfit$ek
}

#extracting the date and the specimen ID from the uniqueIds
dates = substr(uniqueIds, 1, 10)
specimens = substr(uniqueIds, 12, 16)

pmax <- array(NA, c(n,1))
pmax = round(alpha[,1] * ek[,1], digits = 2)


rlc_day_assign <- read.csv("../../../../data/limu/midway/input/date_rlcday/midway_date_rlcday.csv")

rlc_days_by_date = array(dim = length(dates))
hash_of_rlc_days_by_date <- hash(rlc_day_assign$date, rlc_day_assign$rlc_day)
for (i in 1: length(dates)) {
        rlc_days_by_date[i] = hash_of_rlc_days_by_date[[dates[i]]]
}

first_row_of_rlc <- subset(mid_ptools, Epar == 0)
last_row_rlc <- subset(mid_ptools, Epar == 820)

deltaNPQ_ypoint1 <-  (maxNPQ_ypoint1 - minNPQ)

# build the result data frame
result_df <- data.frame(Date = substr(uniqueIds, 1, 10), 
                        "rlc_end_time" = rlc_end_times,
                        "specimenID" = substr(uniqueIds, 12, 16),
                        "uid" = uniqueIds, 
                        "plant_part" = first_row_of_rlc$plant_part,
                        "plantID" = first_row_of_rlc$plantID,
                        "salinity" = first_row_of_rlc$salinity,
                        "treatment" = first_row_of_rlc$treatment,
                        "day_group" = first_row_of_rlc$day_group,
                        "rlc_day" = rlc_days_by_date,
                        "pmax" = pmax,
                        "rETRmax" = rETRMaxes,
                        "rETRmaxYpoint1" = rETRmaxYpoint1,
                        "maxNPQ_Ypoint1" = maxNPQ_ypoint1,
                        "deltaNPQ" = deltaNPQ_ypoint1,
                        "alpha" = round(alpha, digits = 3),
                        "ek" = round(ek, digits = 3),
                        "run" = first_row_of_rlc$run
)


# save to file
write.csv(result_df, "../../../../data/limu/midway/transformed/midway_ek_alpha_normalized_2024.csv")
