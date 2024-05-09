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
