#This script takes the cleaned data from PAM output and runs the Phytotools package to
#produce a final dataset that can be used for analysis.
#take cleaned data from Python script and run through phytotools package to get ek and alpha etc
#open appropriate libraries
library("phytotools")
library("hash")
library("dplyr")

alpha_ek_alga <- read.csv("data_input/hyp_ulva_all_runs_clean.csv", sep = ",")

# add a column that turns date format into POSIXct
alpha_ek_alga$posix_date <- as.POSIXct(alpha_ek_alga$Date, format = "%m/%d/%y")


# add a new column with date and specimen ID as unique key
alpha_ek_alga$uid <- paste(alpha_ek_alga$posix_date, alpha_ek_alga$ID, sep = "_")

# add a new column with a number to indicate just the species
alpha_ek_alga$species <- as.factor(tolower(substr(alpha_ek_alga$ID2, 1, 2)))

#add a new column for treatment assigned from the 5th character in ID2
alpha_ek_alga$treatment <- as.factor(substr(alpha_ek_alga$ID2, 5, 7))

# add a new column for the effective quantum yield
#alpha_ek_alga$QuanYield <- as.factor(round((alpha_ek_alga$Fm. - alpha_ek_alga$F) / alpha_ek_alga$Fm., digits = 3))

# create a new column based on Y(II) at Epar 0 (Effective quantum yield)
alpha_ek_alga <- transform(alpha_ek_alga, QuanYield = ifelse(Epar == "0", Y.II., NA))  

# remove all rows where rETR is null or negative
selectedData <- subset(alpha_ek_alga, rETR > 0) 

# unique function eliminates duplicates returns all unique IDs into a vector
uniqueIds <- unique(selectedData$uid)

# store the number of unique IDs in the variable n
n <- length(uniqueIds)

# create a matrix full of NAs with n rows and 1 column to later calculate ETRmax 
rETRMaxes = array(NA,c(n,1))
rlc_end_times = array(NA, c(n,1))
rETRmaxYpoint1 = array(NA, c(n,1))
#temperatures = array(dim = n)

# create the rETRmax column
#also create a rETRmaxYpoint1 column per Beer & Axelsson where we take the max value but only if Y>0.1
for (i in 1:n){
        light_curve <- subset(selectedData, uid == uniqueIds[i])
        subMaxETR <- max(light_curve$rETR)
        rETRmaxYpoint1[i] <- max(subset(light_curve, y.II > .1)$rETR, na.rm = TRUE)
        #store the subMaxETR in the new column but only in rows where uid is same as uniqueIds of i
        selectedData$rETRmax[selectedData$uid == uniqueIds[i]] <- subMaxETR
        #also store subMaxETR in the matrix rETRMaxes created previously, to later calculate ETRmax
        rETRMaxes[i] = subMaxETR
        #temperatures[i] = max(sub$Temp)
        rlc_end_times[i] <- max(light_curve$Time)
}

# prepare empty matrices to hold output from fitWebb
alpha <- array(NA,c(n,4))
ek <- array(NA,c(n,4))

#create a treatmemt vector from SelecteData treatment column
treatment_names = c("35ppt/0.5uM", "35ppt/14uM", "28ppt/27uM", "28ppt/53uM", "18ppt/53uM", "11ppt/80uM")
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
        
#add a new column to selectedData for ETRmax
selectedData$ETRmax <- round(selectedData$rETRmax*0.5*0.84, digits = 2)

rlc_day_assign <- read.csv("/Users/Angela/src/work/limu/phytotools_alpha_ek/data_input/date_day_assignment.csv")

rlc_days_by_date = array(dim = length(dates))
hash_of_rlc_days_by_date <- hash(rlc_day_assign$Date, rlc_day_assign$RLC.Day)
for (i in 1: length(dates)) {
        rlc_days_by_date[i] = hash_of_rlc_days_by_date[[dates[i]]]
}

lanai_side_by_date = array(dim = length(dates))
hash_of_lanai_side_by_date <- hash(rlc_day_assign$Date, tolower(rlc_day_assign$Lanai.side))
for (i in 1:length(dates)) {
        lanai_side_by_date[i] = as.character(hash_of_lanai_side_by_date[[dates[i]]])
}

first_row_of_rlc <- subset(alpha_ek_alga, Epar == 0 & NPQ == "-")
last_row_rlc <- subset(alpha_ek_alga, Epar == 820)

# build the result data frame
result_df <- data.frame(Date = substr(uniqueIds, 1, 10), 
                        "rlc_end_time" = rlc_end_times,
                        "Specimen ID" = substr(uniqueIds, 12, 18),
                        uid = uniqueIds, 
                        "Plant ID" = first_row_of_rlc$plant.ID,
                        "Species" = first_row_of_rlc$species,
                        "Lanai Side" = first_row_of_rlc$lanai.side,
                        "Treatment" = first_row_of_rlc$treatment,
                        "Temp (°C)" = first_row_of_rlc$Temp,
                        "RLC Order" = first_row_of_rlc$RLC.order,
                        "RLC Day" = rlc_days_by_date,
                        "Run" = first_row_of_rlc$run,
                        "deltaNPQ" = last_row_rlc$deltaNPQ,
                        "pmax" = pmax,
                        "rETRmax" = rETRMaxes,
                        "rETRmaxYpoint1" = rETRmaxYpoint1,
                        "alpha" = round(alpha, digits = 3),
                        "ek" = round(ek, digits = 1)
                        )


# save to file
write.csv(result_df, "data_output/hyp_ulva_all_runs_ek_alpha_normalized.csv")
write.csv (result_df, "../irradiance_ek/data_ek/hyp_ulva_all_runs_ek_alpha_normalized.csv")
write.csv(result_df, "../algal_growth_photosynthesis/data_input/hyp_ulva_all_runs_ek_alpha_normalized.csv")

