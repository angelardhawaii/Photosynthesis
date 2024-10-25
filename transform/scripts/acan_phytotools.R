#This script takes the cleaned data from PAM output and runs the Phytotools package to
#produce a final dataset that can be used for analysis.
#take cleaned data from Python script and run through phytotools package to get ek and alpha etc
#open appropriate libraries
library("phytotools")
library("hash")
library("dplyr")

acan_phyto <- read.csv("../../../../data/limu/acan_2024/transformed/acanthophora_clean.csv", sep = ",")

# add a column that turns date format into POSIXct
acan_phyto$posix_date <- as.POSIXct(acan_phyto$Date, format = "%m/%d/%y")

#add a new column that gets rid of Ac in ID column
acan_phyto$IDnumber <- as.numeric(substr(acan_phyto$ID, 3, 5))

#add a new column that assigns treatment based on the modulus of ID number
acan_phyto$treatment <- ifelse(acan_phyto$IDnumber %% 2 == 0, 1, 0)

#add a new column for temperature based on IDnumber values
acan_phyto$temp <- ifelse(between(acan_phyto$IDnumber, 1, 8) | between(acan_phyto$IDnumber, 17, 24) | 
                                  between(acan_phyto$IDnumber, 33, 40) | 
                                between(acan_phyto$IDnumber, 101, 108) | 
                                  between(acan_phyto$IDnumber, 117, 124) | 
                                  between(acan_phyto$IDnumber, 133, 140), 28, 24)

# add a new column with date and specimen ID as unique key
acan_phyto$uid <- paste(acan_phyto$posix_date, acan_phyto$ID, sep = "_")

acan_phyto$NPQ <- as.numeric(acan_phyto$NPQ)


# add a new column for the effective quantum yield
#acan_phyto$QuanYield <- as.factor(round((acan_phyto$Fm. - acan_phyto$F) / acan_phyto$Fm., digits = 3))

# create a new column based on Y(II) at Epar 0 (Effective quantum yield)
acan_phyto <- transform(acan_phyto, QuanYield = ifelse(Epar == "0", Y.II., NA))



# remove all rows where rETR is null or negative
selectedData <- subset(acan_phyto, rETR > 0) 

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
        minNPQ[i] <- as.numeric(max(subset(light_curve, Epar == 66)$NPQ))
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
treatment_names = c("35ppt/0.5uM", "28ppt/14uM")
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
specimens = substr(uniqueIds, 12, 18)

pmax <- array(NA, c(n,1))
pmax = round(alpha[,1] * ek[,1], digits = 2)


rlc_day_assign <- read.csv("../../../../data/limu/acan_2024/input/acan_date_rlcday.csv")

rlc_days_by_date = array(dim = length(dates))
hash_of_rlc_days_by_date <- hash(rlc_day_assign$Date, rlc_day_assign$RLC_Day)
for (i in 1: length(dates)) {
        rlc_days_by_date[i] = hash_of_rlc_days_by_date[[dates[i]]]
}

first_row_of_rlc <- subset(acan_phyto, Epar == 0)
last_row_rlc <- subset(acan_phyto, Epar == 820)
 
deltaNPQ_ypoint1 <-  (maxNPQ_ypoint1 - minNPQ)

# build the result data frame
result_df <- data.frame(Date = substr(uniqueIds, 1, 10), 
                        "rlc_end_time" = rlc_end_times,
                        "specimenID" = substr(uniqueIds, 12, 18),
                        "uid" = uniqueIds, 
                        "plantID" = first_row_of_rlc$plantID,
                        "treatment" = first_row_of_rlc$treatment,
                        "temperature" = first_row_of_rlc$temp,
                        "rlc_order" = first_row_of_rlc$rlc_order,
                        "rlc_day" = rlc_days_by_date,
                        "deltaNPQ" = last_row_rlc$deltaNPQ,
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
write.csv(result_df, "../../../../data/limu/acan_2024/transformed/acan_ek_alpha_normalized_2024.csv")
