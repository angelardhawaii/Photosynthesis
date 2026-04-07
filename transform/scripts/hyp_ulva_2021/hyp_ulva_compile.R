# Step 2 (after pam data cleaner tidy)
#This script is to compile all of the necessary data into one file for further use
#Necessary because of the change in the pam data cleaner from Python to R script and tidy
#The 2021 dataset was done in a messy way, lots of code to clean it and compile necessary
#After compiling all needed data into the file, pull this into the combo script
# Created November 21, 2025
# By Angela Richards Donà

ps_2021 <- read_csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/transformed/pam_2021data_clean_tidy.csv")
# data merge to add id, fv_fm, and temp columns from the key dataset "ids_all.csv"

ps_2021_key  <- read_csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/input/ids_all.csv")   
# has date, time, id, fv_fm, temp for all up to end of 2021

#add a small helper to fix IDs in both datasets
#(?i) → case-insensitive (so HUl, hul, HUL all match)
#^hul → only at the start of the string
#\\s* → optional spaces after ("HUl 1-1" and "HUl1-1" both work)
clean_ids <- function(df) {
  df %>%
    mutate(
      date = as.Date(date),
      id   = as.character(id),
      id   = trimws(id),
      
      # 1) Remove parentheses (and contents) only on those Feb 2022 dates
      id = if_else(
        date %in% as.Date(c("2022-02-11", "2022-02-15", "2022-02-19")),
        str_replace_all(id, "\\s*\\([^)]*\\)", ""),
        id
      ),
      
      # 2) Normalize case
      id = tolower(id),
      
      # 3) Fix hul → h or u on specific dates
      id = case_when(
        # hul should become h (Hypnea)
        date %in% as.Date(c("2022-10-14", "2022-10-22",
                            "2022-04-21", "2022-04-25", "2022-04-29")) ~
          str_replace(id, "^hul\\s*", "h"),
        
        # hul should become u (Ulva)
        date %in% as.Date(c("2022-10-13", "2022-10-21")) ~
          str_replace(id, "^hul\\s*", "u"),
        
        TRUE ~ id
      )
    )
}

ps_2021    <- clean_ids(ps_2021)
ps_2021_key <- clean_ids(ps_2021_key)

unique(ps_2021$id)
unique(ps_2021_key$id)

# Ensure types & select the join keys/columns we want from the key table
key_clean <- ps_2021_key %>%
  transmute(
    date = as.Date(as.character(date)),   # force factor → character → Date
    time = as_hms(as.character(time)),    # ensure "HH:MM:SS" format
    id    = tolower(as.character(id)),    # ensure character and lowercase
    fv_fm = as.numeric(fv_fm),
    temp  = as.numeric(temp)
  )


# Normalize types in the main table and define block boundaries
ps_2021_prepped <- ps_2021 %>%
  mutate(
    date = as.Date(as.character(date)),   # ensure Date
    time = as_hms(as.character(time))     # ensure HMS
  ) %>%
  arrange(date, time) %>%    # make sure rows are in time order within date
  group_by(date) %>%
  mutate(block = cumsum(epar == 0)) %>%    # each start row (epar==0) begins a new block
  ungroup()


# Join only at the start rows (but easiest is to join all rows on date+time —
#    only the start rows will actually match because only those times exist in key_clean)
ps_2021_joined <- ps_2021_prepped %>%
  left_join(key_clean, by = c("date", "time")) %>%
  group_by(date, block) %>%
  tidyr::fill(id.y, fv_fm.y, temp, .direction = "down") %>%
  ungroup() 

# Carry the id/fv_fm/temp DOWN within each block so all 9 rows get them
ps_2021 <- ps_2021_joined %>%
  group_by(date, block) %>%                      # block from your earlier mutate(cumsum(epar==0))
  tidyr::fill(id.y, fv_fm.y, temp, .direction = "downup") %>%  # fill across all rows in the block
  ungroup() %>%
  select(-block)

# Merge id.x and id.y columns into single id column
ps_2021 <- ps_2021 %>%
  mutate(
    # normalize blanks to NA and trim
    id.x = na_if(trimws(as.character(id.x)), ""),
    id.y = na_if(trimws(as.character(id.y)), ""),
    # merge without losing info
    id = case_when(
      !is.na(id.x) & !is.na(id.y) & id.x != id.y ~ paste(id.x, id.y, sep = "|"),
      !is.na(id.x)                                ~ id.x,
      !is.na(id.y)                                ~ id.y,
      TRUE                                        ~ NA_character_
    )
  ) %>%
  select(-id.x, -id.y)                    # drop the old columns

# Merge fv_fm.x and fv_fm.y columns into single fv_fm column
ps_2021 <- ps_2021 %>%
  mutate(
    fv_fm.x = as.numeric(fv_fm.x),
    fv_fm.y = as.numeric(fv_fm.y),
    fv_fm = case_when(
      !is.na(fv_fm.x) & !is.na(fv_fm.y) & fv_fm.x != fv_fm.y ~ (fv_fm.x + fv_fm.y) / 2,  # average if different
      !is.na(fv_fm.x)                                         ~ fv_fm.x,
      !is.na(fv_fm.y)                                         ~ fv_fm.y,
      TRUE                                                     ~ NA_real_
    )
  ) %>%
  select(-fv_fm.x, -fv_fm.y)            # drop the old columns

ps_2021 <- ps_2021 %>%
  mutate(unique_id = paste(date, tolower(id), sep = "_"))  

# column for run
ps_2021 <- ps_2021 %>%
  mutate(run = case_when(
    date %in% c("2021-09-20", "2021-09-24", "2021-09-28") ~ 1,
    date %in% c("2021-09-21", "2021-09-25", "2021-09-29") ~ 1,
    date %in% c("2021-09-30", "2021-10-04", "2021-10-08") ~ 2,
    date %in% c("2021-10-01", "2021-10-05", "2021-10-09") ~ 2,
    date %in% c("2021-10-11", "2021-10-15", "2021-10-19") ~ 3,
    date %in% c("2021-10-21", "2021-10-25", "2021-10-29") ~ 4,
    date %in% c("2021-11-01", "2021-11-05", "2021-11-09") ~ 5,
    date %in% c("2021-11-04", "2021-11-08", "2021-11-12") ~ 5,
    date %in% c("2022-02-11", "2022-02-15", "2022-02-19") ~ 6, #ulva only
    date %in% c("2022-02-21", "2022-02-25", "2022-03-01") ~ 7, #ulva only
    date %in% c("2022-04-21", "2022-04-25", "2022-04-29") ~ 8, #hyp only
    date %in% c("2022-10-04", "2022-10-08", "2022-10-12") ~ 9, #hyp only
    date %in% c("2022-10-13", "2022-10-21") ~ 10, #ulva only
    date %in% c("2022-10-14", "2022-10-22") ~ 10, #hyp only
    TRUE ~ NA_real_   # fallback for all other dates
  ))

# Make a new run-based column to use for filling in missing plant_id and rlc_order
# for days 1 and 5 as this already exists for day 9 (from merge in hyp_ulva_compile.R)
ps_2021 <- ps_2021 %>%
  mutate(run_proxy = case_when(
    date %in% c("2021-09-20", "2021-09-24", "2021-09-28") ~ 1,
    date %in% c("2021-09-21", "2021-09-25", "2021-09-29") ~ 2,
    date %in% c("2021-09-30", "2021-10-04", "2021-10-08") ~ 3,
    date %in% c("2021-10-01", "2021-10-05", "2021-10-09") ~ 4,
    date %in% c("2021-10-11", "2021-10-15", "2021-10-19") ~ 5,
    date %in% c("2021-10-21", "2021-10-25", "2021-10-29") ~ 6,
    date %in% c("2021-11-01", "2021-11-05", "2021-11-09") ~ 7,
    date %in% c("2021-11-04", "2021-11-08", "2021-11-12") ~ 8,
    date %in% c("2022-02-11", "2022-02-15", "2022-02-19") ~ 9, #ulva only
    date %in% c("2022-02-21", "2022-02-25", "2022-03-01") ~ 10, #ulva only
    date %in% c("2022-04-21", "2022-04-25", "2022-04-29") ~ 11, #hyp only
    date %in% c("2022-10-04", "2022-10-08", "2022-10-12") ~ 12, #hyp only
    date %in% c("2022-10-13", "2022-10-21") ~ 13, #ulva only
    date %in% c("2022-10-14", "2022-10-22") ~ 14, #hyp only
    TRUE ~ NA_real_   # fallback for all other dates
  ))



#More joins with other datasets to get rlc_order, plant_id, treatment, and some temp values
growth_2021 <- read_csv("/Users/angela/src/old_limu/hyp_ulv_2022/data/input/all_runs_growth_011723.csv")

growth_2021 <- clean_names(growth_2021)

growth_2021 <- growth_2021 %>%
  mutate(date = mdy(date))

#clean messy IDs
clean_growth_2021 <- function(df) {
  
  df %>%
    mutate(
      date = as.Date(date),
      
      id = case_when(
        date == as.Date("2022-10-22") ~ str_replace(id, "(?i)^hul\\s*", "h"),
        date == as.Date("2022-04-29") ~ str_replace(id, "(?i)^hul\\s*", "h"),
        date == as.Date("2022-10-21") ~ str_replace(id, "(?i)^hul\\s*", "u"),
        TRUE ~ id
      )
    )
}
growth_2021 <- clean_growth_2021(growth_2021)

growth_2021 <- growth_2021 %>%
  mutate(unique_id = paste(date, tolower(id), sep = "_")) %>%
  mutate(temp = str_sub(temperature, 1, 2))

# remove dead before joining (not present in ps data)
growth_2021 <- growth_2021 %>%
  filter(unique_id != "2021-11-12_hm6-4")  

# column for run
growth_2021 <- growth_2021 %>%
  mutate(run = case_when(
    date == "2021-09-28" ~ 1,
    date == "2021-09-29" ~ 1,
    date == "2021-10-08" ~ 2,
    date == "2021-10-09" ~ 2,
    date == "2021-10-19" ~ 3,
    date == "2021-10-29" ~ 4,
    date == "2021-11-09" ~ 5,
    date == "2021-11-12" ~ 5,
    date == "2022-02-19" ~ 6,
    date == "2022-03-01" ~ 7,
    date == "2022-04-29" ~ 8,
    date == "2022-10-12" ~ 9,
    date == "2022-10-21" ~ 10,
    date == "2022-10-22" ~ 10,
    TRUE ~ NA_real_   # fallback for all other dates
  ))

# Count unique_ids in each dataset
n_ps  <- n_distinct(ps_2021$unique_id)
n_gr  <- n_distinct(growth_2021$unique_id)

# Overlap
n_overlap <- ps_2021 %>%
  distinct(unique_id) %>%
  semi_join(growth_2021 %>% distinct(unique_id),
            by = "unique_id") %>%
  nrow()

cat("Unique IDs in ps_2021:   ", n_ps, "\n")
cat("Unique IDs in growth_2021:", n_gr, "\n")
cat("Overlap:                  ", n_overlap, "\n")

# IDs in growth_2021 but NOT in ps_2021
missing_in_ps <- growth_2021 %>%
  distinct(unique_id) %>%
  anti_join(ps_2021 %>% distinct(unique_id),
            by = "unique_id")

missing_in_ps


#join growth_2021 to ps_2021
# 1) Create a key from growth_2021 with (run, unique_id) → (plant_id, temp, treatment, rlc_order)
growth_key <- growth_2021 %>%
  select(run, unique_id, plant_id, temp, treatment, rlc_order) %>%
  distinct()

# (Optional) sanity checks on the key
# - exactly one row per (run, unique_id)?
dups_in_key <- growth_key %>%
  count(run, unique_id) %>% filter(n > 1)
if (nrow(dups_in_key) > 0) {
  warning("growth_key has duplicate (run, unique_id). Resolve before joining.")
}

# 2) Left-join into ps_2021 by (unique_id)
ps_2021_aug <- ps_2021 %>%
  left_join(growth_key, by = c("unique_id"))

# 3) Quick QA
stopifnot(nrow(ps_2021_aug) == nrow(ps_2021))  # row count unchanged

# Which (unique_id) in ps_2021 didn’t find a match in growth_2021?
missing_pairs <- ps_2021_aug %>%
  distinct(unique_id) %>%
  anti_join(growth_key %>% distinct(run, unique_id), by = c("unique_id"))



# How many rows have missing joins for each joined column?
colSums(is.na(ps_2021_aug[, c("plant_id", "temp.y", "treatment", "rlc_order")]))

## Optional sanity checks
stopifnot(nrow(ps_2021_aug) == nrow(ps_2021))  # row count unchanged

# Any (id) in ps_2021 that didn’t find a match in growth_2021?
missing_in_key <- ps_2021 %>%
  anti_join(growth_key, by = c("unique_id")) %>%
  distinct(unique_id)
nrow(missing_in_key)

# Are there duplicate (unique_id) in the key (should be 1)?
dups_in_key <- growth_key %>%
  count(unique_id) %>%
  filter(n > 1)
nrow(dups_in_key)

#Merge temp.x and temp.y columns into single temp column add unique_id column)
ps_2021 <- ps_2021_aug %>%
  mutate(
    temp.x = as.numeric(temp.x),
    temp.y = as.numeric(temp.y),
    temp = case_when(
      !is.na(temp.x) & !is.na(temp.y) & temp.x != temp.y ~ (temp.x + temp.y) / 2,  # average if different
      !is.na(temp.x)                                       ~ temp.x,
      !is.na(temp.y)                                       ~ temp.y,
      TRUE                                                 ~ NA_real_
    ),
    
    # merge run.x and run.y into a single run column
    run = dplyr::coalesce(run.x, run.y)
  ) %>%
  select(-temp.x, -temp.y, -run.x, -run.y)    # drop the old columns

# Assign plant_id from same unique_id in same run
# rows where plant_id or rlc_order are already known (your "third date" rows)
# helper: return first non-NA, or NA if all NA
first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[1]
}

lookup_ids <- ps_2021 %>%
  group_by(run_proxy, id) %>%
  summarise(
    plant_id  = first_non_na(plant_id),
    rlc_order = first_non_na(rlc_order),
    temp      = first_non_na(temp),
    .groups   = "drop"
  ) %>%
  # keep only combos where at least one of these is known
  filter(!is.na(plant_id) | !is.na(rlc_order) | !is.na(temp))

ps_2021 <- ps_2021 %>%
  left_join(lookup_ids, by = c("run_proxy", "id"),
            suffix = c("", "_src")) %>%
  mutate(
    plant_id  = coalesce(plant_id,  plant_id_src),
    rlc_order = coalesce(rlc_order, rlc_order_src),
    temp      = coalesce(temp,      temp_src)
  ) %>%
  select(-plant_id_src, -rlc_order_src, -temp_src)

# Check: should be exactly ONE row per (run_proxy, id)
lookup_ids %>%
  count(run_proxy, id) %>%
  filter(n > 1)
# → should return 0 rows

#add species column
ps_2021 <- ps_2021 %>%
  mutate(species = str_sub(id, 1, 1))

# Assign treatments for special cases T0 and T2.5 (=2_5)
ps_2021 <- ps_2021 %>%
  mutate(
    date = as.Date(date),
    run  = as.numeric(run),                      # now run exists
    last_digit = str_extract(unique_id, "\\d$"),
    
    treatment = case_when(
      # ---- DATE OVERRIDES ----
      date %in% as.Date(c("2022-10-22", "2022-10-18", "2022-10-14")) ~ "0",
      date %in% as.Date(c("2022-10-21", "2022-10-13", "2022-10-17")) ~ "2_5",
      
      # ---- RUN OVERRIDES ----
      run %in% c(6, 7, 9) ~ "0",
      run == 8            ~ "2_5",
      
      # ---- RUNS 1–5: assign based on last digit ----
      run %in% 1:5 ~ case_when(
        last_digit == "1" ~ "1",
        last_digit == "2" ~ "2",
        last_digit == "3" ~ "3",
        last_digit == "4" ~ "4",
        TRUE ~ NA_character_
      ),
      
      TRUE ~ NA_character_
    )
  )

#missing 2021-10-29_hm4-2 (fixed) and 2021-11-12_hm6-4 (F0 too low so remove)
# but 2021-11-08 is the same individual and is still alive on this date so adding plant_id and rlc_order
ps_2021 <- ps_2021 %>%
  mutate(
    plant_id  = if_else(unique_id == "2021-11-08_hm6-4", 48, plant_id),
    rlc_order = if_else(unique_id == "2021-11-08_hm6-4", 6,  rlc_order)
  )

colSums(is.na(ps_2021[, c("plant_id", "temp", "treatment", "rlc_order")]))

ps_2021 %>%
  filter(is.na(id)) %>%
  count(date)

# Rows missing plant_id or rlc_order
missing_id_rows <- ps_2021 %>%
  filter(is.na(plant_id) | is.na(rlc_order)) %>%
  select(date, run, run_proxy, id, unique_id, plant_id, rlc_order) %>%
  arrange(run_proxy, id, date)

missing_id_rows %>% head(20)

missing_ids <- ps_2021 %>%
  filter(is.na(plant_id) | is.na(rlc_order)) %>%
  distinct(run_proxy, id, unique_id)

nrow(missing_ids)
missing_ids
write.csv(ps_2021, file = "/Users/angela/src/Photosynthesis/data/wastewater/transformed/pam_2021data_compiled.csv")
