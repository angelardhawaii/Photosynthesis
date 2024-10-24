# Photosynthesis
Projects in R and Python, data, scripts, and markdowns for photosynthesis data analysis

This project contains three subprojects that must be run in order: pam-data-cleaner, transform, and analyze. A data folder is also included and sorted by project species and/or year. More detailed README files are to be found in each.

### The Basics:
The data folder.
  This is set up by species/project and each subfolder within "data" has two more folders: input and transformed. All raw data is found in input folders, i.e., irradiance data, PAM data, temperature data, growth data. All data that results from the scripts herein get saved in the transformed folders.
  - hyp_ulv_2022 is a two (macroalgae) species dataset. The hyp_ulv scripts were the first written and may be a bit messier, but they work. Use these if you have two species or two color morphs or two of something you wish to analyze separately. 
  - acan_2024 is a single species dataset. Use the accompany

## The pam-data-cleaner folder (Python)
Step 1. pam-data-cleaner contains a python script that takes data from Win-Control software from Walz and will clean the very messy data of extraneous information. It will also concatenate all .csv files that have been put into a single folder. 
If your rapid light curves are already assigned to a sample ID (cannot do so with diving PAM while underwater), run pam-cleaner-2.0.py, otherwise, run the original version pam-cleaner.py, which will use another file you provide that matches date and time between both files and adds the sample ID to the final clean output file. The clean file gets saved in the "transformed" data folder under the appropriate species project.

## The transform folder (R)
Step 2. Phytotools. The transform project folder contains the transform.Rproj and scripts that...well, transform the data. Script titles that begin with phytotools use the phytotools package by Silsbe and Kromkamp to calculate different values of Ek and alpha, Pmax, depending on the algorithm you choose. Per their 2012 paper, the script uses the Webb et al. algorithm. The output is saved to the appropriate transformed data folder.
NOTE: the Phytotools package is no longer supported by CRAN, so you will need to use the workaround to get the package - see the README inside the transform folder for more details on how to get it working.

Step 3. Hsat. Script titles that begin with hsat pair irradiance data from LiCor irradiance loggers with the dataset that now includes Ek, alpha, and Pmax from the previous step - phytotools. See the README for this folder for more details.

## The analyze folder (R)
  The scripts contained in this folder can be run in any order. These feature mixed effects model scripts for growth, pmax_ek_npq, and hsat_dspi. Markdowns are also saved in each species/project subfolder.

  FOR MORE DETAILS ABOUT WHAT THE SCRIPTS IN EACH FOLDER DO, SEE THE README FILES.
