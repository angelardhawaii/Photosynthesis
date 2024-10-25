#This script takes the cleaned data from PAM output and runs the Phytotools package to
#produce a final dataset that can be used for analysis.
#get ek and alpha etc

#To be used with data that already has all unique IDs
#add column for treatment in next step from this output

#open appropriate libraries
library("phytotools")
library("hash")
library("dplyr")

wastewater_df <- read.csv("data_input/wastewater_clean_sorted.csv", sep = ",")

# add a column that turns date format into POSIXct
wastewater_df$posix_date <- as.POSIXct(wastewater_df$Date, format = "%Y-%m-%d")

# add a new column with date and specimen ID as unique key
wastewater_df$uid <- paste(wastewater_df$posix_date, wastewater_df$ID, sep = "_")

# add a new column with a number to indicate just the species
wastewater_df$species <- as.factor(tolower(substr(wastewater_df$ID, 1, 2)))

#make deltaNPQ a factor
wastewater_df$deltaNPQ <- as.numeric(wastewater_df$deltaNPQ)

# add a new column for the effective quantum yield
wastewater_df$QuanYield <- as.factor(round((wastewater_df$Fm. - wastewater_df$F) / wastewater_df$Fm., digits = 3))

# create a new column based on Y(II) at Epar 0 (Effective quantum yield)
wastewater_df <- transform(wastewater_df, QuanYield = ifelse(Epar == "0", Y.II., NA))  

# remove all rows where rETR is null or negative
selected_df <- subset(wastewater_df, rETR > 0) 

# unique function eliminates duplicates returns all unique IDs into a vector
uniqueIds <- unique(selected_df$uid)

# store the number of unique IDs in the variable n
n <- length(uniqueIds)

# create a matrix full of NAs with n rows and 1 column to later calculate ETRmax 
rETRMaxes = array(NA,c(n,1))
rlc_end_times = array(NA, c(n,1))
rETRmaxYpoint1 = array(NA, c(n,1))
NPQmax = array(NA, c(n,1))
#temperatures = array(dim = n)

# create the rETRmax column
#also create a rETRmaxYpoint1 column per Beer & Axelsson where we take the max value but only if Y>0.1
for (i in 1:n){
        light_curve <- subset(selected_df, uid == uniqueIds[i])
        subMaxETR <- max(light_curve$rETR)
        rETRmaxYpoint1[i] <- max(subset(light_curve, Y.II. > .1)$rETR, na.rm = TRUE)
        subMaxNPQ <- max(light_curve$NPQ)
        NPQmax[i] <- max(subset(light_curve, Y.II. > .1)$NPQ, na.rm = TRUE)
        #store the subMaxETR in the new column but only in rows where uid is same as uniqueIds of i
        selected_df$rETRmax[selected_df$uid == uniqueIds[i]] <- subMaxETR
        #also store subMaxETR in the matrix rETRMaxes created previously, to later calculate ETRmax
        rETRMaxes[i] = subMaxETR
        #temperatures[i] = max(sub$Temp)
        rlc_end_times[i] <- max(light_curve$Time)
}

# prepare empty matrices to hold output from fitWebb
alpha <- array(NA,c(n,4))
ek <- array(NA,c(n,4))

#create a treatmemt vector from SelecteData treatment column
treatment_names = c("35ppt/0.5uM", "27ppt/86uM", "21ppt/53uM", "15ppt/17uM")
treatment = array(NA,c(n,1))


for (i in 1:n){
        #Get ith data
        Epar <- selected_df$Epar[selected_df$uid==uniqueIds[i]]
        rETR <- selected_df$rETR[selected_df$uid==uniqueIds[i]]
        y.II <- selected_df$Y.II.[selected_df$uid==uniqueIds[i]]
        
        #Call function (using y.II and adding normalize = TRUE will yield more accurate PE parameters)
        myfit <- fitWebb(Epar, y.II, normalize = TRUE)
        #store the four values outputs in the matrix
        alpha[i,] <- myfit$alpha
        ek[i,] <- myfit$ek
}

#extracting the date and the specimen ID from the uniqueIds
#dates = substr(uniqueIds, 1, 10)
#specimens = substr(uniqueIds, 12, 18)

pmax <- array(NA, c(n,1))
pmax = round(alpha[,1] * ek[,1], digits = 2)
        
#add a new column to selectedData for ETRmax
selected_df$ETRmax <- round(selected_df$rETRmax*0.5*0.84, digits = 2)

#rlc_day_assign <- read.csv("/Users/Angela/src/work/limu/phytotools_alpha_ek/data_input/date_day_assignment.csv")

#rlc_days_by_date = array(dim = length(dates))
#hash_of_rlc_days_by_date <- hash(rlc_day_assign$Date, rlc_day_assign$RLC.Day)
#for (i in 1: length(dates)) {
#        rlc_days_by_date[i] = hash_of_rlc_days_by_date[[dates[i]]]
#}

#lanai_side_by_date = array(dim = length(dates))
#hash_of_lanai_side_by_date <- hash(rlc_day_assign$Date, tolower(rlc_day_assign$Lanai.side))
#for (i in 1:length(dates)) {
#        lanai_side_by_date[i] = as.character(hash_of_lanai_side_by_date[[dates[i]]])
#}

first_row_of_rlc <- subset(wastewater_df, Epar == 0 & NPQ == "-")
last_row_rlc <- subset(wastewater_df, deltaNPQ > 0)

# build the final data frame
final_df <- data.frame("Date" = first_row_of_rlc$posix_date, 
                        "rlc_end_time" = rlc_end_times,
                        "Specimen ID" = first_row_of_rlc$ID,
                        #"Plant ID" = first_row_of_rlc$plant.ID,
                        "deltaNPQ" = last_row_rlc$deltaNPQ,
                        "Species" = first_row_of_rlc$species,
                        "Treatment" = first_row_of_rlc$treatment,
                        "Temp (°C)" = first_row_of_rlc$temp,
                        "RLC Time" = first_row_of_rlc$Time,
                        "RLC Day" = first_row_of_rlc$rlc_day,
                        "Run" = first_row_of_rlc$run,
                        "pmax" = pmax,
                        "rETRmax" = rETRMaxes,
                        "rETRmaxYpoint1" = rETRmaxYpoint1,
                        "NPQmax" = NPQmax,
                        "alpha" = round(alpha, digits = 3),
                        "ek" = round(ek, digits = 1)
                        )


# save to file
write.csv(final_df, "data_output/wastewater/hyp_ulv_wastewater_ek_alpha_normalized.csv")
write.csv (final_df, "../irradiance_ek/data_ek/wastewater/hyp_ulv_wastewater_ek_alpha_normalized.csv")
write.csv(final_df, "../algal_growth_photosynthesis/data_input/wastewater/hyp_ulv_wastewater_ek_alpha_normalized.csv")

