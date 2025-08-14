# Clean spectral reflectance data files
#Specifically, make sure they all conform to the same naming conventions
# By Angela Richards Dona
# 07/23/25
#

#load some helpful libraries
library(dplyr)
library(tidyr)
library(prospectr)
library(hyperSpec)
#library(pavo)
library(spectrolab)
library(tidyverse)
library(purrr)
library (tibble)
library(janitor)
library(readr)


# Start with one file since data is so extensive
spec_waste <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/input/nearly_raw_spec/hyp_ulv_r2_dia_d9.csv")


#ids are often just repeated, need replicate numbers appended after names cleaned
spec_clean <- spec_waste %>%
  # 1) make syntactic names, turning dots into underscores
  clean_names() %>% 
  select_if( ~ any(!is.na(.)) ) %>% #get rid of empty columns
  
  # 2) for any column that is just letters+digits (no "_"), append "_4"
  #this is because some measurements may have more than three replicates
  rename_with(
    .cols = matches("^[a-z]\\d+$"),
    .fn   = ~ paste0(.x, "_4")
  )%>%
  # keep only wavelengths between 350.086 and 749.1 nm
  filter(wavelength >= 350.086,
         wavelength <= 750.7)

# basic transpose of dataframe so wavelengths are columns and ids are rows
spec_transposed <- spec_clean %>%
  # move wavelength into the rownames
  column_to_rownames("wavelength") %>%
  # transpose the matrix
  t() %>%
  # back to a data.frame
  as.data.frame(check.names = FALSE) %>%
  # make the old rownames (h125_1, h125_2…) into a column called id
  rownames_to_column("id")


write.csv(spec_transposed, file = "../data/wastewater/input/spec/hu_r2_dia_d9.csv")
#these are being written and saved into the same file folder as the og.





#check to make sure hu_r1_d9 is same as these. This dataset is a bit different from Malin
spec_waste_r1 <- read.csv("/Users/angela/src/Photosynthesis/data/wastewater/input/nearly_raw_spec/hyp_ulv_r1_dia_d9.csv")

#ids are often just repeated, need replicate numbers appended after names cleaned
spec_clean_r1 <- spec_waste_r1 %>%
  # 1) make syntactic names, turning dots into underscores
  clean_names() %>% 
  select_if( ~ any(!is.na(.)) ) %>% #get rid of empty columns
 
  # keep only wavelengths between 350.086 and 749.1 nm
  filter(wavelength >= 350.086,
         wavelength <= 750.699)

# basic transpose of dataframe so wavelengths are columns and ids are rows
spec_transposed_r1 <- spec_clean_r1 %>%
  # move wavelength into the rownames
  column_to_rownames("wavelength") %>%
  # transpose the matrix
  t() %>%
  # back to a data.frame
  as.data.frame(check.names = FALSE) %>%
  # make the old rownames (h125_1, h125_2…) into a column called id
  rownames_to_column("id")

write.csv(spec_transposed_r1, file = "../data/wastewater/input/spec/hu_r1_dia_d9.csv")
