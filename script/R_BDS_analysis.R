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
library(cowplot)
library(pracma)

#1. Load data ####

#set path directories for process_batch_files function.
#by default use project directory /csv

batch_meta_dir <- "./csv/"
bds100_dir <- file.path(batch_meta_dir, "BDS100")
bds250_dir <- file.path(batch_meta_dir, "BDS250")
rapid_imp_dir <- file.path(batch_meta_dir, "RAPID_IMP")
rapid_hig_dir <- file.path(batch_meta_dir, "RAPID_HIG")

#1.1 Misc ####

#Remove a specific row
# batch_summary <- batch_summary %>%
#   slice(-2)

Filter and save specific sensor dataset
C590909160221 <- data_250hz %>%
  filter(long_id == "C590909160221" & !is.na(passage_point))

write.csv(C590909160221, "C590909160221.csv", row.names = FALSE)

#Filter by variable conditions
# vertical bar (|) = OR  & = AND
# != not equal to
# == is equal to
# batch_summary <- batch_summary %>%
#   filter(pres_processed != 'N'| acc_processed != 'N')

# Import a specific dataset
# batch_summary=read_csv("./R_output_files/4_batch_summary.csv")

#There is a bug with the BDSAnalysisTool, use this to assign passage manually if needed.
#100_imp = sensor class, not dataframe name
#assign_passage_points_manual("100_imp", batch_summary)
#time_normalization_manual("100_imp", batch_summary)

#Prototype tool to be updated with
#Calculate and add passage duration
#e.g., total passage (inj. - tail), injection to nadir, nadir to tail,
#max impact within 1s nadir
#usability features e.g., user can end process at any point, current calculations saved to batch_summary, pop-out viewer window

#Expectation with pump
#1. Atmospheric pressure when sensor calibrates in air (~1000mbar)
#2. Static pressure when BDS exits injector pipe into water column should increase from atmospheric pressure
#   but will be minimal
#3. Static pressure will possibly remain unaffected prior to pump interaction as the pipework is a fixed diameter
#4. The impeller of the pump creates a low-pressure zone to draw water into the pump, so pressure should drop below
#   atmospheric pressure (nadir). Velocity increases = pressure decreases
#5. Pressure will then increase at discharge, and then drop to a similar static pressure to pre-suction
#6. Pressure will then increase when the sensors enter the recapture tank due to water depth




#2. BDS analysis tools ####

#Run process_batch_files to consolidate data frames by sensor sample rate, remove unnecessary columns and set variable structure 
#requires batch_summary file to be present
#ensure previously processed batch_summary has been moved

process_batch_files(batch_meta_dir, bds100_dir, bds250_dir, rapid_imp_dir, rapid_hig_dir)

#Run prototype BDS analysis tool for: 
#nadir identification and validation, max pressure identification (1s) and validation
#Atmospheric/static pressure validation, rate pressure change calculation
#Ratio pressure change calculation, log ratio pressure change calculation
#acceleration magnitude identification (0.1s), injection and tailwater identification
#passage point assignment, time series normalization

BDSAnalysisTool(batch_summary, data_250hz, data_100hz, data_100_imp, data_2000_hig, "100_imp")

#After analysis, add treatment variables
add_treatment(data_250hz, "data_250hz")

#When batch complete, crop data to passage points, save and export
#If continuing with analysis, remove NA passage_points to reduce memory consumption from data files
save_data(batch_summary, data_100hz, data_250hz, data_100_imp, data_2000_hig)

#3. PumpFlow data analysis ####

##3.1 500RPM 100% BEP ####

###3.1.1 consolidate data ####

#Video metrics
# 108 validated blade strike video -> 95 recovered recovered data sets -> 79 complete data with video (82 - 3 NV)
# -> 44 clean passage with video -> 28 blade strike with video -> 7 hub contact with video
# An additional 6 blade strike videos could be correlated with incomplete sensor data to enhance blade strike comparison statistics (e.g., acceleration only)

#Data pre-processing steps performed:
#1. All data recovered from sensors and validated against deployment (e.g., non-recoverable removed).
# 113 deployed -> 95 recovered data
#2. All deployed sensors correlated with video footage, passage and strike conditions validated
# 113 deployed -> 108 validated in video
#3. Recovered sensor data (95) pre-processed using BDS_import.py. Output plots checked for sensor errors.
# 95 recovered data -> 82 complete data for analysis (pres, acceleration)
#4. batch_summary file cleaned to just contain 82 complete data for analysis
#5. process_batch_files function used to import and consolidate data sets

#Processed 500_1 consolidated files saved
# processed_500_1 <- data_100_imp
# write.csv(processed_500_1, "./PROCESSED_data/processed_500_1.csv", row.names = FALSE)

#500_1 data set batch summary saved
# processed_500_1_batch_summary <- batch_summary
# write.csv(processed_500_1_batch_summary, "./PROCESSED_data/processed_500_1_batch_summary.csv", row.names = FALSE)

###3.1.2 Process ROI ####

#Import pre-processed data sets required for ROI and pressure/accelration metrics

data_100_imp=read_csv("./PROCESSED_data/processed_500_1.csv")
batch_summary=read_csv("./PROCESSED_data/processed_500_1_batch_summary.csv")

#Process using BDSAnalysisTool

BDSAnalysisTool(batch_summary, data_250hz, data_100hz, data_100_imp, data_2000_hig, "100_imp")

#Apply manual passage/time normalisation due to bug in tool
assign_passage_points_manual("250", batch_summary)
time_normalization_manual("250", batch_summary)

#Save and exported ROI processed data frame, either for continued processing, or for further analysis
save_data(batch_summary, data_100hz, data_250hz, data_100_imp, data_2000_hig)

###3.1.3 Analysis ####

process_batch_files(batch_meta_dir, bds100_dir, bds250_dir, rapid_imp_dir, rapid_hig_dir)

BDSAnalysisTool(batch_summary, data_250hz, data_100hz, data_100_imp, data_2000_hig, "250")


write.csv(data_250hz, "./PROCESSED_data/processed_KATHA.csv", row.names = FALSE)


