# Acanthophora spicifera populate_ek

ek = read.csv("../../../../data/limu/acan_2024/transformed/acan_ek_alpha_normalized_2024.csv")
# Make  sure the date is loaded as date
ek$posix_date <- as.POSIXct(ek$Date, format = "%Y-%m-%d")
ek$rlc_order <- as.factor(ek$rlc_order)
ek$rlc_day <- as.factor(ek$rlc_day)
ek$treatment <- as.character(ek$treatment)
ek$plantID <- as.factor(ek$plantID)
ek$run <- as.character(ek$run)
ek$pmax_min <- ek$pmax * 60
