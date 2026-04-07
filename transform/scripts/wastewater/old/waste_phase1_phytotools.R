#This script takes the cleaned data from PAM data-cleaner output and runs the Phytotools package to
#produce a final dataset that can be used for analysis with ek and alpha
#Dataset is from phase I wastewater with 8 nitrate/phosphate levels to determine saturation points 
#for Hypnea and Ulva
#By Angela Richards Donà
#Date: 01/16/25


#open appropriate libraries
library("phytotools")
library("hash")
library("dplyr")
library("hms")
library("lubridate")

waste_p1_ptools <- read.csv("../data/wastewater/transformed/waste_p1_2025_clean.csv", sep = ",")
waste_p1_growth <- read.csv("../data/wastewater/input/growth_2025/waste_phase1_growth_all.csv")

# add a column that turns date format into POSIXct
waste_p1_ptools$Date <- ymd(waste_p1_ptools$Date)

#same for time
waste_p1_ptools$Time <- as_hms(waste_p1_ptools$Time)

#add a new column that gets assigns species based on the first letter in ID column
#then add a new column that joins the growth dataset for plant number and treatment columns
waste_p1_ptools <- waste_p1_ptools %>%
  mutate(species = as.character(substr(ID, 1, 1))) %>%
  left_join(waste_p1_growth, join_by(ID))

##RLC-order is added to dataset in model script

# add a new column with date and specimen ID as unique key
waste_p1_ptools$uid <- paste(waste_p1_ptools$Date, waste_p1_ptools$ID, sep = "_")

waste_p1_ptools$NPQ <- as.numeric(waste_p1_ptools$NPQ)

# create a new column based on Y(II) at Epar 0 (Effective quantum yield)
waste_p1_ptools <- transform(waste_p1_ptools, QuanYield = ifelse(Epar == "0", Y.II., NA))

# remove all rows where rETR is null or negative
selectedData <- subset(waste_p1_ptools, rETR > 0) 

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
treatment_names = c("80 μM", "140 μM", "245 μM", "428 μM", "748 μM", "1308 μM", "2287 μM", "4000 μM")
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

first_row_of_rlc <- subset(waste_p1_ptools, Epar == 0)
last_row_rlc <- subset(waste_p1_ptools, Epar == 420)

deltaNPQ_ypoint1 <-  (maxNPQ_ypoint1 - minNPQ)



# build the result data frame
result_df <- data.frame(Date = substr(uniqueIds, 1, 10), 
                        "rlc_end_time" = rlc_end_times,
                        "specimenID" = substr(uniqueIds, 12, 16),
                        "uid" = uniqueIds,
                        "species" = first_row_of_rlc$species,
                        "plantID" = first_row_of_rlc$plant,
                        "treatment" = first_row_of_rlc$treatment,
                        "pmax" = pmax,
                        "rETRmax" = rETRMaxes,
                        "rETRmaxYpoint1" = rETRmaxYpoint1,
                        "maxNPQ_Ypoint1" = maxNPQ_ypoint1,
                        "deltaNPQ" = deltaNPQ_ypoint1,
                        "alpha" = round(alpha, digits = 3),
                        "ek" = round(ek, digits = 3),
                        "FvFm" = first_row_of_rlc$Fv.Fm,
                        "growth" = first_row_of_rlc$growth
)


# save to file
write.csv(result_df, "../data/wastewater/transformed/waste_p1_ek_alpha_normalized_2025.csv")
