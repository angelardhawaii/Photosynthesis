# Chondria tumulosa populate_ek

ek = read.csv("../data/midway_2024/transformed/midway_ek_alpha_normalized_2024.csv")
ek$date <- ymd(ek$Date) # Make  sure the date is loaded as date
ek$rlc_end_time <- as_hms(ek$rlc_end_time) # make sure time is hms
ek$rlc_day <- as.factor(ek$rlc_day)
ek$nitrate <- as.factor(as.character(ek$treatment))
ek$salinity <- as.factor(ek$salinity)
ek$plantID <- as.factor(ek$plantID)
ek$run <- as.character(ek$run)
ek$pmax_min <- ek$pmax * 60
ek$deltaNPQ <- as.numeric(ek$deltaNPQ)

#assign new column for chronological ranking of individuals in each day_group
ek <- ek %>%
  group_by(day_group) %>%
  mutate(rank = rank(rlc_end_time)) %>%
  ungroup()

#use ranking to make smaller groups of ~15 minutes for rlc_order
ek <- ek %>%
  group_by(day_group) %>%
  mutate(rlc_order = floor((rank -1)/3) +1)
