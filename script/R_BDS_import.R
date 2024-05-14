#Python installation and operation within R/RStudio

#1. Python installation####
#Run only for first time installation 

#Install reticulate to use python within R
install.packages("reticulate")

#check package operation
library(reticulate)

#install miniconda
install_miniconda()

#Verify python installation
py_config()

#installing required python packages for txt > csv conversion

py_install("jupyter_client")
py_install("matplotlib")
py_install("pandas")
py_install("quaternion")

#2. Python operation####

##2.1 Run python####

#Used when working interactively with python scripts. Not needed if just running the conversion script

#Function for running python via runtasks.R
#Load the reticulate package for integrating python into R
#Start a Python REPL (Read-Eval-Print Loop) session within the R environment

run_python <- function() {
  source('./script/runtasks.R')
}

run_python()

##2.2 Restart python####

#Often necessary after making changes to objects/attributes

restart_python <- function() {
  rstudioapi::restartSession(command = "source('./script/runtasks.R')")
}

# Call the restart_python function
restart_python()

##2.3 Delete pycache####

#Delete pycache when needed
delete_pycache <- function() {
  pycache_path <- file.path("./script", "__pycache__")
  if (dir.exists(pycache_path)) {
    unlink(pycache_path, recursive = TRUE)
    cat("Deleted __pycache__ folder.\n")
  } else {
    cat("__pycache__ folder not found.\n")
  }
}

delete_pycache()

##2.4 Misc python####

##Set Python Search Path: Add the search path for the current working directory in Python. 
# import sys
# import os
# current_dir = os.getcwd()
# if current_dir not in sys.path:
#     sys.path.append(current_dir)
# for root, dirs, files in os.walk(current_dir):
#   for name in dirs:
#     dir_path = os.path.join(root, name)
#   if dir_path not in sys.path:
#     sys.path.append(dir_path)

#Check if needed
#import sys; print(sys.path)

# acc_5_ = (data['accmag'] >= 5).sum()
# acc_5g = ((data['accmag'] >= 5) & (data['time'].diff().fillna(float('inf')) >= 0.02)).sum()


# summary_csv_path = self.dir_csv / (self.filename.stem + "_summary.csv")
# summary_data.to_csv(summary_csv_path, index=False)

# def plot_acceleration_overview(self, save: bool = True, show: bool = False) -> None:
#     t = self.data["time"][::10]  # Downsampling for visualization
#     accmag = self.data["accmag"].rolling(10).mean()[::10]  # Rolling average for smoothing
# 
#     fig, ax = plt.subplots(figsize=(25, 5))
#     ax.set_xlabel("time [s]")
#     ax.set_ylabel("Acceleration magnitude [g]")
#     ax.plot(t, accmag, color="C1")
#     ax.tick_params(axis="y", labelcolor="C1")
#     fig.tight_layout()
# 
#     if save:
#         new_filename =  self.filename.stem + "_acc"
#         plt.savefig((self.dir_plots / new_filename).with_suffix(".png"))
#     if show:
#         plt.show()
#     plt.close()
#     
# def plot_magnetic_overview(self, save: bool = True, show: bool = False) -> None:
#     t = self.data["time"][::10]  # Downsampling for visualization
# 
#     fig, ax = plt.subplots(figsize=(25, 5))
#     ax.set_xlabel("time [s]")
#     ax.set_ylabel("Magnetic flux density [mT]")
#     
#     ax.plot(t, self.data["magx"].rolling(10).mean()[::10], color="C1", label="magx")
#     ax.plot(t, self.data["magy"].rolling(10).mean()[::10], color="C2", label="magy")
#     ax.plot(t, self.data["magz"].rolling(10).mean()[::10], color="C3", label="magz")
#     
#     ax.legend()
#     fig.tight_layout()
# 
#     if save:
#         new_filename = self.filename.stem + "_mag" 
#         plt.savefig((self.dir_plots / new_filename).with_suffix(".png"))
#     if show:
#         plt.show()
#     plt.close()
#     


#3. BDS txt -> CSV conversion####

#Add BDS 100hz and BDS 250hz txt files to RAW_data/BDS100 and RAW_data/BDS250
#Run BDS_import script to convert raw txt files to csv, calculate average pressure (P1, P2, P3 - mbar),
#calculate acceleration magnitude (g-force, gravity removed), check pressure sensor operation incl. warnings,
#check time series incl. warnings, plot pressure sensor checks, plot pressure and acceleration overview,
#plot magnetic flux (BDS100), create batch_summary containing sensor metadata
#output files stored in ./csv and ./plots

reticulate::source_python('./script/BDS_import.py')

#4. BDS analysis functions####

# Function to find the batch summary file
load_batch_summary <- function(batch_meta_dir) {
  cat("Loading batch summary file...\n")
  batch_summary_files <- list.files(batch_meta_dir, pattern = "_batch_summary\\.csv$", full.names = TRUE)
  if (length(batch_summary_files) == 0) {
    stop("No batch summary file found in the specified directory.")
  }
  batch_summary_path <- batch_summary_files[1]  # Assuming there's only one such file
  batch_summary <- read_csv(batch_summary_path, show_col_types = FALSE)
  cat("Batch summary loaded successfully.\n")
  return(batch_summary)
}

# Function to process 100hz data files
process_data_100hz <- function(file_path, sample_rate, long_id, short_id) {
  cat(paste("Processing 100hz data file:", file_path, "\n"))
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      long_id = factor(long_id),  # sensor_id is set to the value of 'file' from batch_summary
      short_id = factor(short_id), # short_id is set to the value of 'sensor' from batch_summary
      sample_rate = factor(sample_rate),
      time = as.numeric(time), 
      pres = as.numeric(pres), 
      accmag = as.numeric(accmag)  
    ) %>%
    select(-c(P1, P2, P3, absaccx, absaccy, absaccz))
  return(data)
}

# Function to process 250hz data files
process_data_250hz <- function(file_path, sample_rate, long_id, short_id) {
  cat(paste("Processing 250hz data file:", file_path, "\n"))
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      long_id = factor(long_id),  # sensor_id is set to the value of 'file' from batch_summary
      short_id = factor(short_id), # short_id is set to the value of 'sensor' from batch_summary
      sample_rate = factor(sample_rate),
      time = as.numeric(time), 
      pres = as.numeric(pres), 
      accmag = as.numeric(accmag)
    )%>%
    select(-c(P1, P2, P3))
  return(data)
}

#process data 
process_batch_files <- function(batch_meta_dir, bds100_dir, bds250_dir) {
  batch_summary <- load_batch_summary(batch_meta_dir)
  data_250hz <- tibble()
  data_100hz <- tibble()
  cat("Processing batch files...\n")
  for (i in seq_len(nrow(batch_summary))) {
    long_id <- batch_summary$file[i]
    short_id <- batch_summary$sensor[i]
    sample_rate <- as.integer(gsub("hz", "", batch_summary$class[i]))
    subdir <- ifelse(sample_rate == 250, bds250_dir, bds100_dir)
    file_path <- file.path(subdir, paste0(long_id, ".csv"))
    
    if (file.exists(file_path)) {
      if (sample_rate == 250) {
        data_250hz <- bind_rows(data_250hz, process_data_250hz(file_path, sample_rate, long_id, short_id))
      } else if (sample_rate == 100) {
        data_100hz <- bind_rows(data_100hz, process_data_100hz(file_path, sample_rate, long_id, short_id))
      }
    } else {
      cat(paste("File not found:", file_path, "\n"))
    }
  }
  assign("batch_summary", batch_summary, envir = .GlobalEnv)
  assign("data_250hz", data_250hz, envir = .GlobalEnv)
  assign("data_100hz", data_100hz, envir = .GlobalEnv)
  cat("Batch file processing completed.\n")
}

add_treatment <- function(data, data_name) {
  # Define the possible entries
  pump_rpm_choices <- c(400, 500, 600, 700)
  inj_pos_choices <- c('us_pump_1_c', 'us_pump_1_h', 'us_pump_1_l', 'us_pump_2_c', 'us_pump_2_h', 'us_pump_2_l', 'ds_pump')
  treatment_choices <- c(
    '400_control', '400_impact_us1c', '400_impact_us1h', '400_impact_us1l', '400_impact_us2c', '400_impact_us2h', '400_impact_us2l', 
    '500_control', '500_impact_us1c', '500_impact_us1h', '500_impact_us1l', '500_impact_us2c', '500_impact_us2h', '500_impact_us2l',
    '600_control', '600_impact_us1c', '600_impact_us1h', '600_impact_us1l', '600_impact_us2c', '600_impact_us2h', '600_impact_us2l',
    '700_control', '700_impact_us1c', '700_impact_us1h', '700_impact_us1l', '700_impact_us2c', '700_impact_us2h', '700_impact_us2l'
  )
  pump_rpm <- pump_rpm_choices[menu(pump_rpm_choices, title = "Select Pump RPM:")]
  inj_pos <- inj_pos_choices[menu(inj_pos_choices, title = "Select Injection position:")]
  treatment <- treatment_choices[menu(treatment_choices, title = "Select Treatment:")]
  cat("Selected values:\n")
  cat("pump_rpm:", pump_rpm, "\n")
  cat("inj_pos:", inj_pos, "\n")
  cat("treatment:", treatment, "\n")
  data <- data %>%
    mutate(
      pump_rpm = factor(pump_rpm),
      inj_pos = factor(inj_pos),
      treatment = factor(treatment)
    )
  assign(data_name, data, envir = .GlobalEnv)
  cat(paste("Treatment variables assigned to", data_name, "\n"))
}

# Function to subset data by unique levels in 'long_id' and assign to global environment
subset_by_long_id <- function(data) {
  unique_ids <- unique(data$long_id)
  
  for (id in unique_ids) {
    subset_data <- data %>% filter(long_id == id)
    subset_name <- as.character(id)
    assign(subset_name, subset_data, envir = .GlobalEnv)
    cat(paste("Created subset:", subset_name, "\n"))
  }
}

