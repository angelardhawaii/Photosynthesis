# populate_ek

ek = read.csv("./data_ek/hyp_ulva_all_runs_ek_alpha_normalized.csv")
# Make  sure the date is loaded as date
ek$posix_date <- as.POSIXct(ek$Date, format = "%Y-%m-%d")
ek$RLC.Order <- as.factor(ek$RLC.Order)
ek$RLC.Day <- as.factor(ek$RLC.Day)
ek$Treatment <- as.character(ek$Treatment)
ek$Plant.ID <- as.factor(ek$Plant.ID)
ek$Run <- as.character(ek$Run)
ek$rETRmaxYpoint1_min <- as.numeric(ek$rETRmaxYpoint1 * 60)
ek$rETRmaxYpoint1 <- as.numeric(ek$rETRmaxYpoint1)
ek$pmax_min <- ek$pmax * 60