#R 4.3.3
#RStudio 2023.12.1+402 "Ocean Storm" Release for windows
#BDS data analysis for PumpFlow testrig (May 2024)

###NOTE###
#Before continuing, scrutinize pre-processing output plots and batch_summary to determine 
#usable data sets for import and analysis. Remove unwanted datasets from BDS100/BDS250 folders
###NOTE###

#Required libraries####

library(tidyverse)

#1. Load data ####

#set path directories for process_batch_files function.
#by default use project directory /csv

batch_meta_dir <- "./csv/"
bds100_dir <- file.path(batch_meta_dir, "BDS100")
bds250_dir <- file.path(batch_meta_dir, "BDS250")

#Run process_batch_files to consolidate data frames by sensor sample rate, remove unnecessary columns and set variable structure 
process_batch_files(batch_meta_dir, bds100_dir, bds250_dir)

#add treatment variables to 100hz and 250hz consolidated files
add_treatment(data_250hz, "data_250hz")
add_treatment(data_100hz, "data_100hz")

#subset wrangled data frames by dataset name
subset_by_long_id(data_100hz)
subset_by_long_id(data_250hz)

#processed data now ready for basic plots and ROI analysis

#2. Visualization and ROI analysis ####








#create plots and determine ROI
# Could automate determining where the nadir occurs and then plot time either side of this region 
# So use the nadir time from batch_summary to determine this
# do the same for max acceleration, but this will miss multiple acceleration events, so ROI analysis is required, and note the accleration times
# to create plots 5s either side of acceleration, for example
#filter data sets based on ROI
#determine number of impact events (e.g., 1s apart?)
#determine severity of impact events
#determine nadir, max pressure before
#determine rate of change and ratio change
#determine LPR
#plot lpr between operational scenarios
#need a normalize time step to create the final plots where nadir = 0.5 and data around = 0 - 1
# the same is required for acceleration too 