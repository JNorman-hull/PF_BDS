#R 4.3.3
#RStudio 2023.12.1+402 "Ocean Storm" Release for windows
#BDS data analysis for PumpFlow testrig (May 2024)

###NOTE###
#Before continuing, scrutinize pre-processing output plots and batch_summary to determine 
#usable data sets for import and analysis. Remove unwanted datasets from BDS100/BDS250 folders
###NOTE###

#Required libraries####

library(tidyverse)
library(plotly)
library(pracma)

#1. Load data ####

#set path directories for process_batch_files function.
#by default use project directory /csv

batch_meta_dir <- "./csv/"
bds100_dir <- file.path(batch_meta_dir, "BDS100")
bds250_dir <- file.path(batch_meta_dir, "BDS250")

#Run process_batch_files to consolidate data frames by sensor sample rate, remove unnecessary columns and set variable structure 
#requires batch_summary file to be present
#ensure previously processed batch_summary has been moved

process_batch_files(batch_meta_dir, bds100_dir, bds250_dir)

#2. BDS analysis ####

#Run prototype BDS analysis tool for: 
#nadir identification and validation, max pressure identification (1s) and validation
#Atmospheric/static pressure validation, rate pressure change calculation
#Ratio pressure change calculation, log ratio pressure change calculation
#acceleration magnitude identification (0.1s), injection and tailwater identification
#passage point assignment, time series normalization

BDSAnalysisTool(batch_summary, data_250hz, data_100hz, 250)

#After analysis, add treatment variables
add_treatment(data_250hz, "data_250hz")

#When batch complete, crop data to passage points, save and export
#If continuing with analysis, remove NA passage_points to reduce memory consumption from data files
save_data(batch_summary, data_100hz, data_250hz)

#Prototype tool to be updated with
#Calculate and add passage duration
#e.g., total passage (inj. - tail), injection to nadir, nadir to tail,
#max impact within 1s nadir
#usability features e.g., user can end process at any point, current calculations saved to batch_summary, pop-out viewer window

#subset wrangled data frames by dataset name
#subset_by_long_id(data_250hz)

#Expectation with pump
#1. Atmospheric pressure when sensor calibrates in air (~1000mbar)
#2. Static pressure when BDS exits injector pipe into water column should increase from atmospheric pressure
#   but will be minimal
#3. Static pressure will possibly remain unaffected prior to pump interaction as the pipework is a fixed diameter
#4. The impeller of the pump creates a low-pressure zone to draw water into the pump, so pressure should drop below
#   atmospheric pressure (nadir). Velocity increases = pressure decreases
#5. Pressure will then increase at discharge, and then drop to a similar static pressure to pre-suction
#6. Pressure will then increase when the sensors enter the recapture tank due to water depth