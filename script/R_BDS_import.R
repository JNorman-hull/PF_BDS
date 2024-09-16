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
py_install("opencv")
py_install("statsmodels")


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


#summary_data = data.iloc[::summary_interval].reset_index(drop=True)
#enable saving for checking summary_data if needed
#summary_csv_path = self.dir_csv / (self.filename.stem + "_summary.csv")
#summary_data.to_csv(summary_csv_path, index=False)

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

# def plot_acc_mag_overview(self, save: bool = True, show: bool = False) -> None:
#     """Plot acceleration magnitude with magnetic flux for 100hz sensors."""
#     t = self.data["time"][::10]
#     
#     fig, ax1 = plt.subplots(figsize=(25, 5))
#     ax1.set_xlabel("time [s]")
#     
#     color = "C0"
#     ax1.set_ylabel("Acceleration magnitude [g]", color=color)
#     ax1.plot(t, self.data["accmag"].rolling(10).mean()[::10], color=color, label="Acceleration")
#     ax1.tick_params(axis='y', labelcolor=color)
#     ax1.ticklabel_format(useOffset=False)
# 
#     ax2 = ax1.twinx()
#     color = "C1"
#     ax2.set_ylabel("Magnetic flux density [mT]", color=color)
#     ax2.plot(t, self.data["magx"].rolling(10).mean()[::10], color="C1", label="magx")
#     ax2.plot(t, self.data["magy"].rolling(10).mean()[::10], color="C2", label="magy")
#     ax2.plot(t, self.data["magz"].rolling(10).mean()[::10], color="C3", label="magz")
#     ax2.tick_params(axis='y', labelcolor=color)
#     
#     lines, labels = ax1.get_legend_handles_labels()
#     lines2, labels2 = ax2.get_legend_handles_labels()
#     ax2.legend(lines + lines2, labels + labels2, loc='upper right')
#     
#     fig.tight_layout()
# 
#     if save == True:
#         new_filename = self.filename.stem + "_acc_mag" 
#         plt.savefig((self.dir_plots / new_filename).with_suffix(".png"))
#     if show == True:
#         plt.show()
#     plt.close()   


#3. BDS txt -> CSV conversion####

#Add BDS 100hz and BDS 250hz txt files to RAW_data/BDS100 and RAW_data/BDS250
#Run BDS_import script to convert raw txt files to csv, calculate average pressure (P1, P2, P3 - mbar),
#calculate acceleration magnitude (g-force, gravity removed), check pressure sensor operation incl. warnings,
#check time series incl. warnings, plot pressure sensor checks, plot pressure and acceleration overview,
#plot magnetic flux (BDS100), create batch_summary containing sensor metadata
#output files stored in ./csv and ./plots

reticulate::source_python('./script/BDS_import.py')

#4. BDS analysis functions####



theme_JN <- function(base_size=10){ 
  theme_grey() %+replace%
    theme(
      axis.text = element_text(colour="black"),
      axis.title = element_text(colour="black"),
      axis.ticks = element_line(colour="black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(colour = "black",fill = NA),
      panel.spacing.x = unit(12, "pt")
    ) 
}

load_batch_summary <- function(batch_meta_dir) {
  cat("Loading batch summary file...\n")
  batch_summary_files <- list.files(batch_meta_dir, pattern = "_batch_summary\\.csv$", full.names = TRUE)
  if (length(batch_summary_files) == 0) {
    stop("No batch summary file found in the specified directory.")
  }
  batch_summary_path <- batch_summary_files[1]  # Assuming there's only one such file
  batch_summary <- read_csv(batch_summary_path, show_col_types = FALSE)
  
  cols_to_numeric <- c(
    "pres_min[mbar]", "pres_min[time]", "pres_max[mbar]",
    "pres_max[time]", "acc_max[m/s2]", "acc_max[time]",
    "HIG_max[g]", "HIG_max[time]")
  
  for (col in cols_to_numeric) {
    if (col %in% colnames(batch_summary)) {
      batch_summary[[col]] <- as.numeric(batch_summary[[col]])
    } else {
      batch_summary[[col]] <- NA_real_  # Initialize with NA if the column does not exist
    }
  }
  
  additional_cols <- c(
    "t_max_p_1s_nadir","t_nadir", "nadir",
    "max_p_1s_nadir", "rate_pc", "ratio_pc", "log_RPC")
  for (col in additional_cols) {
    batch_summary[[col]] <- as.numeric(0)
  }
  
  # Add 'processed' column as a factor with default value "N"
  batch_summary$pres_processed <- as.character("N")
  batch_summary$acc_processed <- as.character("N")
  batch_summary$badsens <- as.character("N")
  
  cat("Batch summary loaded successfully.\n")
  return(batch_summary)
}

process_data_100hz <- function(file_path, sample_rate, long_id, short_id) {
  cat(paste("Processing 100hz data file:", file_path, "\n"))
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      long_id = factor(long_id),  # long_id is set to the value of 'file' from batch_summary
      short_id = factor(short_id), # short_id is set to the value of 'sensor' from batch_summary
      sample_rate = factor(sample_rate),
      time = as.numeric(time), 
      pres = as.numeric(pres), 
      accmag = as.numeric(accmag)  
    ) %>%
    select(-c(P1, P2, P3, absaccx, absaccy, absaccz))
  return(data)
}

process_data_250hz <- function(file_path, sample_rate, long_id, short_id) {
  cat(paste("Processing 250hz data file:", file_path, "\n"))
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      long_id = factor(long_id),  # long_id is set to the value of 'file' from batch_summary in process_batch_files
      short_id = factor(short_id), # short_id is set to the value of 'sensor' from batch_summary
      sample_rate = factor(sample_rate),
      time = as.numeric(time), 
      pres = as.numeric(pres), 
      accmag = as.numeric(accmag)
    )%>%
    select(-c(P1, P2, P3))
  return(data)
}

process_data_100_imp <- function(file_path, sample_rate, long_id, short_id) {
  cat(paste("Processing RAPID 100hz IMP data file:", file_path, "\n"))
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      long_id = factor(long_id),  # long_id is set to the value of 'file' from batch_summary
      short_id = factor(short_id), # short_id is set to the value of 'sensor' from batch_summary
      sample_rate = factor(sample_rate),
      time = as.numeric(`Time (s)`), 
      pres = as.numeric(`Pressure (mbar)`), 
      accmag = as.numeric(`Accel_Mag (g)`)
    )%>%
    select(-c('Time (s)', 'Pressure (mbar)', 'Accel_Mag (g)', 'P_Temp (C)', 'Battery (V)')
           )
  return(data)
}

process_data_2000_hig <- function(file_path, sample_rate, long_id, short_id) {
  cat(paste("Processing RAPID 2000hz HIG data file:", file_path, "\n"))
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      long_id = factor(long_id),  # long_id is set to the value of 'file' from batch_summary
      short_id = factor(short_id), # short_id is set to the value of 'sensor' from batch_summary
      sample_rate = factor(sample_rate),
      time = as.numeric(`Time (s)`), 
      HIGaccmag = as.numeric(`HIGAccel_Mag (g)`)
    )%>%
    select(-c('Time (s)', 'HIGAccel_Mag (g)')
    )
  return(data)
}


process_batch_files <- function(batch_meta_dir, bds100_dir, bds250_dir, rapid_imp_dir, rapid_hig_dir) {
  batch_summary <- load_batch_summary(batch_meta_dir)
  data_250hz <- tibble()
  data_100hz <- tibble()
  data_100_imp <- tibble()
  data_2000_hig <- tibble()
  cat("Processing batch files...\n")
  
  for (i in seq_len(nrow(batch_summary))) {
    long_id <- batch_summary$file[i]
    short_id <- batch_summary$sensor[i]
    sample_rate <- batch_summary$class[i]
    
    if (sample_rate == "250hz") {
      sample_rate_value <- 250
      subdir <- bds250_dir
    } else if (sample_rate == "100hz") {
      sample_rate_value <- 100
      subdir <- bds100_dir
    } else if (sample_rate == "100_imp") {
      sample_rate_value <- "100_imp"
      subdir <- rapid_imp_dir
    } else if (sample_rate == "2000_hig") {
      sample_rate_value <- "2000_hig"
      subdir <- rapid_hig_dir
    } else {
      cat("Invalid sample rate:", sample_rate, "\n")
      next
    }
    
    file_path <- file.path(subdir, paste0(long_id, ".csv"))
    
    if (file.exists(file_path)) {
      if (sample_rate_value == 250) {
        data_250hz <- bind_rows(data_250hz, process_data_250hz(file_path, sample_rate_value, long_id, short_id))
      } else if (sample_rate_value == 100) {
        data_100hz <- bind_rows(data_100hz, process_data_100hz(file_path, sample_rate_value, long_id, short_id))
      } else if (sample_rate_value == "100_imp") {
        data_100_imp <- bind_rows(data_100_imp, process_data_100_imp(file_path, sample_rate_value, long_id, short_id))
      } else if (sample_rate_value == "2000_hig") {
        data_2000_hig <- bind_rows(data_2000_hig, process_data_2000_hig(file_path, sample_rate_value, long_id, short_id))
      }
    } else {
      cat(paste("File not found:", file_path, "\n"))
    }
  }
  
  assign("batch_summary", batch_summary, envir = .GlobalEnv)
  assign("data_250hz", data_250hz, envir = .GlobalEnv)
  assign("data_100hz", data_100hz, envir = .GlobalEnv)
  assign("data_100_imp", data_100_imp, envir = .GlobalEnv)
  assign("data_2000_hig", data_2000_hig, envir = .GlobalEnv)
  cat("Batch file processing completed.\n")
}


filter_batch_summary <- function(batch_summary, sample_rate) {
  if (sample_rate == 250) {
    data_class <- "250hz"
    num_rows <- 96
  } else if (sample_rate == 100) {
    data_class <- "100hz"
    num_rows <- 96
  } else if (sample_rate == "100_imp") {
    data_class <- "100_imp"
    num_rows <- 96
  } else if (sample_rate == "2000_hig") {
    data_class <- "2000_hig"
    num_rows <- 2000
  } else {
    stop("Invalid sample rate")
  }
  
  #check that the hig is 2k rows
  # num_rows is 96 for both sensors as 96 rows = 1s @100hz & 250hz. Only included for future redundancy 
  
  # Modify this to allow for either to be Y but not both
  sensors <- batch_summary %>%
    filter(class == data_class & (pres_processed == "N" | acc_processed == "N") & badsens == "N") %>%
    select(file) %>%
    distinct()
  
  if (nrow(sensors) == 0) {
    return(list(message = "All BDS files already processed\n", sensors = NULL, num_rows = num_rows))
  }
  
  list(message = NULL, sensors = sensors, num_rows = num_rows)
}

prompt_user_to_select_sensor <- function(sensors) {
  sensor_choices <- sensors$file
  cat("Select a sensor dataset to begin BDS analysis:\n")
  for (i in seq_along(sensor_choices)) {
    cat(i, ": ", sensor_choices[i], "\n", sep = "")
  }
  selected_sensor_index <- as.integer(readline(prompt = "Enter the number corresponding to the sensor: "))
  if (is.na(selected_sensor_index) || selected_sensor_index < 1 || selected_sensor_index > length(sensor_choices)) {
    stop("Invalid sensor selection")
  }
  sensor_choices[selected_sensor_index]
}

create_combined_plot <- function(data, sensor_summary, selected_sensor, stage, num_rows) {
  # Filter data for the selected sensor
  data <- filter(data, long_id == selected_sensor)
  
  if (stage == 2) {
    nadir_time <- sensor_summary$t_nadir
    nadir_index <- which(data$time == nadir_time)
    row_start <- max(1, nadir_index - num_rows * 2)
    row_end <- min(nrow(data), nadir_index + num_rows * 2)
    data <- data[row_start:row_end, ]
  }
  
  if (stage == 4 || stage == 5) {
    data <- filter(data, !is.na(passage_point))
  }
  
  # Base plot
  p <- plot_ly(data)
  
  # Always add pressure line
  p <- p %>%
    add_lines(x = if(stage == 5) ~time_norm else ~time, y = ~pres, name = 'Pressure', yaxis = 'y1', line = list(color = 'red', width = 1))
  
  # Add acceleration line for stages 3 and above
  if (stage >= 3) {
    p <- p %>%
      add_lines(x = if(stage == 5) ~time_norm else ~time, y = ~accmag, name = 'Acceleration', yaxis = 'y2', line = list(color = 'blue', width = 1))
  }
  
  # Add stage-specific elements
  if (stage == 2) {
    nadir_index <- which(data$time == sensor_summary$t_nadir)
    time_start <- data$time[max(1, nadir_index - num_rows)]
    time_end <- data$time[min(nrow(data), nadir_index + num_rows)]
    vertical_lines <- list(
      list(type = "line", x0 = time_start, x1 = time_start, xref = "x", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", dash = "dash")),
      list(type = "line", x0 = time_end, x1 = time_end, xref = "x", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", dash = "dash"))
    )
    p <- p %>%
      layout(
        shapes = c(vertical_lines, list(
          list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 1, y1 = 1, yref = "paper", line = list(color = "black", width = 1)),
          list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0, yref = "paper", line = list(color = "black", width = 1)),
          list(type = "line", x0 = 0, x1 = 0, xref = "paper", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1)),
          list(type = "line", x0 = 1, x1 = 1, xref = "paper", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1))
        )),
        annotations = list(
          list(x = (time_start + time_end) / 2, y = 0.8, text = '1s time window', showarrow = FALSE, xref = 'x', yref = 'paper', font = list(color = 'black', size = 12), xanchor = 'center', yanchor = 'middle')
        )
      )
  }
  
  if (stage == 4) {
    # ROI 2
    if ("max_accmag_roi2" %in% names(sensor_summary) && "max_accmag_time_roi2" %in% names(sensor_summary)) {
      p <- p %>%
        add_markers(
          x = sensor_summary[["max_accmag_time_roi2"]],
          y = sensor_summary[["max_accmag_roi2"]],
          name = "Max Acc nadir region",
          yaxis = 'y2',
          marker = list(color = 'orange', size = 10, line = list(color = 'black', width = 1))
        )
    }
    
    # ROI 3 and 4
    for (roi in 3:4) {
      max_accmag_col <- paste0("max_accmag_roi", roi)
      event_time_col <- paste0("t_event_roi", roi)
      if (max_accmag_col %in% names(sensor_summary) && event_time_col %in% names(sensor_summary)) {
        p <- p %>%
          add_markers(
            x = sensor_summary[[event_time_col]],
            y = sensor_summary[[max_accmag_col]],
            name = paste0("Max Acc ROI ", roi),
            yaxis = 'y2',
            marker = list(color = 'orange', size = 10, line = list(color = 'black', width = 1))
          )
      }
    }
  }
  
  # Add markers for pressure nadir and max 1s nadir
  if (stage == 1) {
    p <- p %>%
      add_markers(
        x = sensor_summary[["pres_min[time]"]],
        y = sensor_summary[["pres_min[mbar]"]],
        name = 'Pressure nadir',
        yaxis = 'y1',
        marker = list(color = 'red', size = 10, line = list(color = 'black', width = 1))
      ) %>%
      add_annotations(
        x = sensor_summary[["pres_min[time]"]],
        y = sensor_summary[["pres_min[mbar]"]],
        text = 'Pressure nadir',
        showarrow = TRUE,
        arrowhead = 2,
        ax = 150,
        ay = -50,
        yanchor = "top",
        font = list(color = 'red')
      )
  } else if (stage >= 2 && stage <= 4) {
    p <- p %>%
      add_markers(
        x = sensor_summary[["t_nadir"]],
        y = sensor_summary[["nadir"]],
        name = 'Pressure nadir',
        yaxis = 'y1',
        marker = list(color = 'red', size = 10, line = list(color = 'black', width = 1))
      ) %>%
      add_annotations(
        x = sensor_summary[["t_nadir"]],
        y = sensor_summary[["nadir"]],
        text = 'Pressure nadir',
        showarrow = TRUE,
        arrowhead = 2,
        ax = 150,
        ay = -50,
        yanchor = "top",
        font = list(color = 'red')
      ) %>%
      add_markers(
        x = sensor_summary[["t_max_p_1s_nadir"]],
        y = sensor_summary[["max_p_1s_nadir"]],
        name = 'Max 1s nadir',
        yaxis = 'y1',
        marker = list(color = 'blue', size = 10, line = list(color = 'black', width = 1))
      ) %>%
      add_annotations(
        x = sensor_summary[["t_max_p_1s_nadir"]],
        y = sensor_summary[["max_p_1s_nadir"]],
        text = 'Max 1s nadir',
        showarrow = TRUE,
        arrowhead = 2,
        ax = 150,
        ay = -50,
        yanchor = "top",
        font = list(color = 'blue')
      )
  }
  
  # Layout
  p <- p %>%
    layout(
      xaxis = list(
        title = if(stage == 5) "Normalized Time" else "Time [s]",
        zeroline = FALSE,
        showline = FALSE,
        showgrid = FALSE,
        ticks = 'outside',
        tickcolor = 'black'
      ),
      yaxis = list(
        title = "Pressure [mbar]",
        side = 'left',
        zeroline = FALSE,
        showgrid = FALSE,
        showline = FALSE,
        ticks = 'outside',
        tickcolor = 'black',
        color = 'red'
      ),
      yaxis2 = list(
        title = "Acceleration magnitude [m/s2]",
        overlaying = 'y',
        side = 'right',
        showgrid = FALSE,
        zeroline = FALSE,
        showline = FALSE,
        ticks = 'outside',
        tickcolor = 'black',
        color = 'blue',
        automargin = TRUE 
      ),
      plot_bgcolor = 'white',
      paper_bgcolor = 'white',
      showlegend = FALSE,
      shapes = list(
        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 1, y1 = 1, yref = "paper", line = list(color = "black", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0, yref = "paper", line = list(color = "black", width = 1)),
        list(type = "line", x0 = 0, x1 = 0, xref = "paper", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1)),
        list(type = "line", x0 = 1, x1 = 1, xref = "paper", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1))
      )
    )
  
  return(p)
}

check_sensor <- function(batch_summary, selected_sensor) {
  is_sens_good <- readline(prompt = "Is sensor data appropriate for analysis? (Y/N): ")
  
  if (toupper(is_sens_good) == "Y") {
    cat("Sensor data good\n")
    return(TRUE)
  } else {
    cat("Sensor data not appropriate for analysis\n")
    batch_summary <<- batch_summary %>%
      mutate(badsens = ifelse(file == selected_sensor, "Y", badsens))
    return(FALSE)
  }
}

process_nadir_value <- function(batch_summary, selected_sensor) {
  is_nadir_correct <- readline(prompt = "Is nadir value correct? (Y/N): ")
  
  if (toupper(is_nadir_correct) == "Y") {
    batch_summary <- batch_summary %>%
      mutate(
        nadir = ifelse(file == selected_sensor, `pres_min[mbar]`, nadir),
        t_nadir = ifelse(file == selected_sensor, `pres_min[time]`, t_nadir)
      )
    cat("Nadir values recorded in batch summary\n")
  } else if (toupper(is_nadir_correct) == "N") {
    new_nadir_value <- as.numeric(readline(prompt = "Enter new nadir value: "))
    new_nadir_time <- as.numeric(readline(prompt = "Enter new nadir time: "))
    if (is.na(new_nadir_value) || is.na(new_nadir_time)) {
      stop("Invalid nadir value or time")
    }
    
    batch_summary <- batch_summary %>%
      mutate(
        nadir = ifelse(file == selected_sensor, new_nadir_value, nadir),
        t_nadir = ifelse(file == selected_sensor, new_nadir_time, t_nadir)
      )
    cat("New nadir values recorded in batch summary\n")
  } else {
    cat("Invalid input. Please enter 'Y' or 'N'.\n")
  }
  
  batch_summary
}

calculate_max_pressure_before_nadir <- function(data, batch_summary, selected_sensor, num_rows) {
  t_nadir_value <- batch_summary %>% filter(file == selected_sensor) %>% pull(t_nadir)
  nadir_index <- which(data$time == t_nadir_value & data$long_id == selected_sensor)
  
  if (length(nadir_index) == 0) {
    cat("Error: t_nadir value not found in the data for the selected sensor.\n")
    return(NULL)
  } else if (length(nadir_index) > 1) {
    cat("Error: Multiple t_nadir values found in the data for the selected sensor. Using the first occurrence.\n")
    nadir_index <- nadir_index[1]
  }
  
  start_index <- max(nadir_index - num_rows, 1)
  end_index <- nadir_index - 1
  
  if (start_index <= end_index) {
    max_p_1s_nadir_value <- max(data$pres[start_index:end_index], na.rm = TRUE)
    max_p_1s_nadir_time <- data$time[which.max(data$pres[start_index:end_index]) + start_index - 1]
  } else {
    max_p_1s_nadir_value <- NA
    max_p_1s_nadir_time <- NA
  }
  
  batch_summary <- batch_summary %>%
    mutate(
      max_p_1s_nadir = ifelse(file == selected_sensor, max_p_1s_nadir_value, max_p_1s_nadir),
      t_max_p_1s_nadir = ifelse(file == selected_sensor, max_p_1s_nadir_time, t_max_p_1s_nadir)
    )
  cat("Max pressure within 1s prior to nadir calculated and recorded in batch summary\n")
  
  batch_summary
}

process_max_p_1s_nadir <- function(batch_summary, selected_sensor) {
  # Retrieve the current max_p_1s_nadir and t_max_p_1s_nadir values
  sensor_summary <- batch_summary %>% filter(file == selected_sensor)
  current_max_p_1s_nadir <- sensor_summary$max_p_1s_nadir
  current_t_max_p_1s_nadir <- sensor_summary$t_max_p_1s_nadir
  
  # Display current values and prompt user for validation
  cat("Current max pressure within 1s prior to nadir value:", current_max_p_1s_nadir, "\n")
  cat("Current time of max pressure within 1s prior to nadir:", current_t_max_p_1s_nadir, "\n")
  is_correct <- readline(prompt = "Is max pressure within 1s prior to nadir value correct? (Y/N): ")
  
  if (toupper(is_correct) == "N") {
    new_max_p_1s_nadir <- as.numeric(readline(prompt = "Enter new max pressure within 1s prior to nadir value: "))
    new_t_max_p_1s_nadir <- as.numeric(readline(prompt = "Enter new time of max pressure within 1s prior to nadir: "))
    
    # Update batch_summary with new values
    batch_summary <- batch_summary %>%
      mutate(
        max_p_1s_nadir = ifelse(file == selected_sensor, new_max_p_1s_nadir, max_p_1s_nadir),
        t_max_p_1s_nadir = ifelse(file == selected_sensor, new_t_max_p_1s_nadir, t_max_p_1s_nadir)
      )
    cat("New max pressure within 1s prior to nadir and time recorded in batch summary\n")
  } else if (toupper(is_correct) == "Y") {
    cat("Max pressure within 1s prior to nadir value and time remain unchanged\n")
  } else {
    cat("Invalid input. Please enter 'Y' or 'N'.\n")
  }
  
  return(batch_summary)
}

calculate_rate_pressure_change <- function(batch_summary, selected_sensor) {
  rate_pc_value <- batch_summary %>% filter(file == selected_sensor) %>% 
    mutate(rate_pc = max_p_1s_nadir - nadir) %>% pull(rate_pc)
  
  batch_summary <- batch_summary %>%
    mutate(
      rate_pc = ifelse(file == selected_sensor, rate_pc_value, rate_pc)
    )
  
  cat("Rate pressure change =", rate_pc_value, "\n")
  cat("Rate pressure change calculated and recorded in batch summary\n")
  
  batch_summary
}

set_max_acclimation_pressure <- function(default_pressure = 1000) {
  cat("Default maximum acclimation pressure is set to", default_pressure, "mbar.\n")
  change_max_acclim <- readline(prompt = "Do you want to change the maximum acclimation pressure? (Y/N): ")
  
  if (toupper(change_max_acclim) == "Y") {
    max_acclim <- as.numeric(readline(prompt = "Enter new maximum acclimation pressure: "))
    while (is.na(max_acclim)) {
      cat("Invalid input. Please enter a numeric value.\n")
      max_acclim <- as.numeric(readline(prompt = "Enter new maximum acclimation pressure: "))
    }
  } else {
    max_acclim <- default_pressure
  }
  
  cat("Maximum acclimation pressure set to", max_acclim, "mbar.\n")
  return(max_acclim)
}

calculate_ratio_pressure_change <- function(batch_summary, selected_sensor, max_acclim) {
  # Calculate ratio pressure change
  ratio_pc_value <- batch_summary %>% filter(file == selected_sensor) %>% 
    mutate(ratio_pc = max_acclim / nadir) %>% pull(ratio_pc)
  
  # Calculate log ratio pressure change
  log_RPC_value <- log(ratio_pc_value)
  
  # Update batch summary with the calculated values
  batch_summary <- batch_summary %>%
    mutate(
      ratio_pc = ifelse(file == selected_sensor, ratio_pc_value, ratio_pc),
      log_RPC = ifelse(file == selected_sensor, log_RPC_value, log_RPC)
    )
  
  # Print calculated values
  cat("Ratio pressure change =", ratio_pc_value, "\n")
  cat("Log ratio pressure change =", log_RPC_value, "\n")
  cat("Ratio and log ratio pressure change calculated and recorded in batch summary\n")
  
  batch_summary
}

prompt_for_injection_tailwater_times <- function(data, batch_summary, selected_sensor, num_rows) {
  sensor_data <- data %>% filter(long_id == selected_sensor)
  t_nadir_val <- batch_summary %>% filter(file == selected_sensor) %>% pull(t_nadir)
  nadir_index <- which(sensor_data$time == t_nadir_val)
  
  if (length(nadir_index) == 0) {
    stop("Nadir time not found in the data")
  }
  
  # Create ROI columns if they don't exist
  roi_columns <- c(
    paste0("t_start_roi", 1:4),
    paste0("t_end_roi", 1:4),
    "t_event_roi3", "t_event_roi4"
  )
  
  for (col in roi_columns) {
    if (!(col %in% colnames(batch_summary))) {
      batch_summary[[col]] <- NA_real_
    }
  }
  
  for (roi in 1:4) {
    cat(paste("\nROI", roi, ":"))
    
    if (roi == 1) {
      cat(" (Main overall sensor passage ROI)\n")
      row_start <- max(1, nadir_index - num_rows * 3)
      row_end <- min(nrow(sensor_data), nadir_index + num_rows * 3)
      t_start_val <- sensor_data$time[row_start]
      t_end_val <- sensor_data$time[row_end]
      cat("Default start time (nadir - 3s):", t_start_val, "\n")
      cat("Default end time (nadir + 3s):", t_end_val, "\n")
      prompt <- "Use these default times for the main passage ROI? (Y/N): "
    } else if (roi == 2) {
      cat(" (Nadir event ROI)\n")
      row_start <- max(1, nadir_index - 20)
      row_end <- min(nrow(sensor_data), nadir_index + 20)
      t_start_val <- sensor_data$time[row_start]
      t_end_val <- sensor_data$time[row_end]
      cat("Default nadir event start time (nadir - 0.2s):", t_start_val, "\n")
      cat("Default nadir event end time (nadir + 0.2s):", t_end_val, "\n")
      prompt <- "Use these default times for nadir event? (Y/N): "
    } else {
      event_label <- if(roi == 3) "pre-nadir" else "post-nadir"
      cat("\nEnter event time series value for", event_label, "event\n")
      event_time <- as.numeric(readline(prompt = "Enter event time value: "))
      event_index <- which.min(abs(sensor_data$time - event_time))
      row_start <- max(1, event_index - 10)
      row_end <- min(nrow(sensor_data), event_index + 10)
      t_start_val <- sensor_data$time[row_start]
      t_end_val <- sensor_data$time[row_end]
      cat("Default", event_label, "start time (event time - 0.1s):", t_start_val, "\n")
      cat("Default", event_label, "end time (event time + 0.1s):", t_end_val, "\n")
      prompt <- paste("Use these default times for", event_label, "event? (Y/N): ")
      
      # Record event time in batch_summary
      batch_summary <- batch_summary %>%
        mutate(
          !!paste0("t_event_roi", roi) := ifelse(file == selected_sensor, event_time, !!sym(paste0("t_event_roi", roi)))
        )
    }
    
    use_default <- readline(prompt = prompt)
    
    if (toupper(use_default) != "Y") {
      cat(paste("Enter custom start and end times for ROI", roi, "\n"))
      t_start_val <- as.numeric(readline(prompt = "Enter start time value: "))
      t_end_val <- as.numeric(readline(prompt = "Enter end time value: "))
      
      while (is.na(t_start_val) || is.na(t_end_val) || t_start_val >= t_end_val) {
        cat("Invalid input. Please ensure start time is less than end time and both are numeric.\n")
        t_start_val <- as.numeric(readline(prompt = "Enter start time value: "))
        t_end_val <- as.numeric(readline(prompt = "Enter end time value: "))
      }
    }
    
    batch_summary <- batch_summary %>%
      mutate(
        !!paste0("t_start_roi", roi) := ifelse(file == selected_sensor, t_start_val, !!sym(paste0("t_start_roi", roi))),
        !!paste0("t_end_roi", roi) := ifelse(file == selected_sensor, t_end_val, !!sym(paste0("t_end_roi", roi)))
      )
  }
  
  cat("Start and end times for all ROIs recorded in batch summary\n")
  
  batch_summary
}

assign_passage_points <- function(data, batch_summary, selected_sensor, sample_rate) {
  cat("Assigning passage points for sensor:", selected_sensor, "\n")
  
  # Filter the data for the selected sensor
  data_sensor <- data %>% filter(long_id == selected_sensor)
  
  # ROI 1
  t_start_val <- batch_summary %>% 
    filter(file == selected_sensor) %>% 
    pull(t_start_roi1)
  t_nadir_val <- batch_summary %>% 
    filter(file == selected_sensor) %>% 
    pull(t_nadir)
  t_end_val <- batch_summary %>% 
    filter(file == selected_sensor) %>% 
    pull(t_end_roi1)
  
  data_sensor <- data_sensor %>%
    mutate(passage_point = case_when(
      time >= t_start_val & time < t_nadir_val ~ "pre_nadir",
      time >= t_nadir_val & time <= t_end_val ~ "post_nadir",
      TRUE ~ NA_character_
    ))
  
  # ROI 2, 3, and 4
  for (roi in 2:4) {
    t_start_val <- batch_summary %>% 
      filter(file == selected_sensor) %>% 
      pull(!!paste0("t_start_roi", roi))
    t_end_val <- batch_summary %>% 
      filter(file == selected_sensor) %>% 
      pull(!!paste0("t_end_roi", roi))
    
    data_sensor <- data_sensor %>%
      mutate(!!paste0("roi", roi) := case_when(
        time >= t_start_val & time <= t_end_val ~ 
          case_when(
            roi == 2 ~ "nadir_event",
            roi == 3 ~ "pre_nadir_event",
            roi == 4 ~ "post_nadir_event",
            TRUE ~ NA_character_
          ),
        TRUE ~ NA_character_
      ))
  }
  
  cat("Passage points assigned to", selected_sensor, "for all ROIs\n")
  
  # Update the data for the selected sensor
  new_columns <- setdiff(names(data_sensor), names(data))
  if (length(new_columns) > 0) {
    for (col in new_columns) {
      data[[col]] <- NA
    }
  }
  
  rows_to_update <- which(data$long_id == selected_sensor)
  data[rows_to_update, names(data_sensor)] <- data_sensor
  
  return(data)
}


max_acceleration_extract <- function(data, batch_summary, selected_sensor) {
  # Filter data for the selected sensor
  sensor_data <- data %>% filter(long_id == selected_sensor)
  
  # Create columns if they don't exist
  new_columns <- c(
    paste0("max_accmag_roi", 1:4),
    paste0("max_accmag_time_roi", 1:4)
  )
  for (col in new_columns) {
    if (!(col %in% colnames(batch_summary))) {
      batch_summary[[col]] <- NA_real_
    }
  }
  
  # Extract max acceleration for ROI 1 and 2
  for (roi in 1:2) {
    t_start <- batch_summary %>% 
      filter(file == selected_sensor) %>% 
      pull(!!paste0("t_start_roi", roi))
    t_end <- batch_summary %>% 
      filter(file == selected_sensor) %>% 
      pull(!!paste0("t_end_roi", roi))
    
    roi_data <- sensor_data %>% 
      filter(time >= t_start & time <= t_end)
    
    max_accmag <- max(roi_data$accmag, na.rm = TRUE)
    max_accmag_time <- roi_data$time[which.max(roi_data$accmag)]
    
    # Update batch_summary
    batch_summary <- batch_summary %>%
      mutate(
        !!paste0("max_accmag_roi", roi) := ifelse(file == selected_sensor, max_accmag, !!sym(paste0("max_accmag_roi", roi))),
        !!paste0("max_accmag_time_roi", roi) := ifelse(file == selected_sensor, max_accmag_time, !!sym(paste0("max_accmag_time_roi", roi)))
      )
  }
  
  # Extract accmag values for ROI 3 and 4
  for (roi in 3:4) {
    t_event <- batch_summary %>% 
      filter(file == selected_sensor) %>% 
      pull(!!paste0("t_event_roi", roi))
    
    event_index <- which.min(abs(sensor_data$time - t_event))
    event_accmag <- sensor_data$accmag[event_index]
    
    # Update batch_summary
    batch_summary <- batch_summary %>%
      mutate(
        !!paste0("max_accmag_roi", roi) := ifelse(file == selected_sensor, event_accmag, !!sym(paste0("max_accmag_roi", roi))),
        !!paste0("max_accmag_time_roi", roi) := ifelse(file == selected_sensor, t_event, !!sym(paste0("max_accmag_time_roi", roi)))
      )
  }
  
  # Print results to console
  cat("\nMaximum acceleration values:\n")
  cat("ROI 1 (Main overall sensor passage) max acceleration:", batch_summary$max_accmag_roi1[batch_summary$file == selected_sensor], "\n")
  cat("ROI 2 (nadir event) max acceleration:", batch_summary$max_accmag_roi2[batch_summary$file == selected_sensor], "\n")
  cat("ROI 3 (pre-nadir event) acceleration:", batch_summary$max_accmag_roi3[batch_summary$file == selected_sensor], "\n")
  cat("ROI 4 (post-nadir event) acceleration:", batch_summary$max_accmag_roi4[batch_summary$file == selected_sensor], "\n")
  
  # Remove max_accmag_time_roi columns for ROI 1, 3, and 4, but keep ROI 2
  batch_summary <- batch_summary %>%
    select(-matches("max_accmag_time_roi[134]"))
  
  return(batch_summary)
}


time_normalization <- function(data, selected_sensor) {
  
  # Ensure time_norm column exists and is initialized as 0
  if (!"time_norm" %in% colnames(data)) {
    data$time_norm <- 0
  }
  
  # Filter data for the selected sensor
  data_sensor <- filter(data, long_id == selected_sensor)
  
  # Identify the start and end times for normalization
  start_time <- min(data_sensor$time[data_sensor$passage_point == 'pre_nadir'], na.rm = TRUE)
  mid_time <- min(data_sensor$time[data_sensor$passage_point == 'post_nadir'], na.rm = TRUE)
  end_time <- max(data_sensor$time[data_sensor$passage_point == 'post_nadir'], na.rm = TRUE)
  
  
  # Calculate normalized time values
  data_sensor <- data_sensor %>%
    mutate(time_norm_new = case_when(
      time <= start_time ~ 0,
      time >= end_time ~ 1,
      time > start_time & time < mid_time ~ (time - start_time) / (mid_time - start_time) * 0.5,
      time >= mid_time & time <= end_time ~ 0.5 + (time - mid_time) / (end_time - mid_time) * 0.5
    ))
  
  # Merge the normalized time back into the original data
  data <- data %>%
    left_join(select(data_sensor, time, long_id, time_norm_new), by = c("time", "long_id")) %>%
    mutate(time_norm = coalesce(time_norm_new, time_norm)) %>%
    select(-time_norm_new)
  
  cat("Injection to tailwater passage time normalization complete\n")
  return(data)
}

update_global_environment <- function(data, batch_summary, sample_rate) {
  if (sample_rate == 250) {
    data_250hz <- data
    assign("data_250hz", data_250hz, envir = .GlobalEnv)
  } else if (sample_rate == 100) {
    data_100hz <- data
    assign("data_100hz", data_100hz, envir = .GlobalEnv)
  } else if (sample_rate == "100_imp") {
    data_100_imp <- data
    assign("data_100_imp", data_100_imp, envir = .GlobalEnv)
  } else if (sample_rate == "2000_hig") {
    data_2000_hig <- data
    assign("data_2000_hig", data_2000_hig, envir = .GlobalEnv)
  }
  assign("batch_summary", batch_summary, envir = .GlobalEnv)
}


#BDS analysis tool prototype
BDSAnalysisTool <- function(batch_summary, data_250hz, data_100hz, data_100_imp, data_2000_hig, sample_rate) {
  
  filter_result <- filter_batch_summary(batch_summary, sample_rate)
  if (!is.null(filter_result$message)) {
    cat(filter_result$message)
    cat("BDS Analysis tool ended.\n")
    return(invisible())
  }
  
  # Update the data assignment logic
  data <- if (sample_rate == 250) {
    data_250hz
  } else if (sample_rate == 100) {
    data_100hz
  } else if (sample_rate == "100_imp") {
    data_100_imp
  } else if (sample_rate == "2000_hig") {
    data_2000_hig
  }
  
  # Create list of sensors available for processing.
  # Add num_rows for calculating max pressure <1s nadir
  sensors <- filter_result$sensors
  num_rows <- filter_result$num_rows
  # Prompt user to select a sensor
  selected_sensor <- prompt_user_to_select_sensor(sensors)
  # Create sensor summary
  sensor_summary <- batch_summary %>% filter(file == selected_sensor)
  time_range <- max(data$time) - min(data$time)
  
  # Draw plot for nadir value check
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 1)
  print(plotly_plot)
  
  # Check the sensor for good data
  result <- check_sensor(batch_summary, selected_sensor)
  
  #rest of function is within a loop determined by result of check_sensor function
  
  if (!result) {
    cat("Stopping BDS analysis...\n")
  } else {
    cat("Proceeding with BDS analysis...\n")
    
    # Process nadir value. User can input new value if needed
    batch_summary <- process_nadir_value(batch_summary, selected_sensor)
    
    # Calculate max pressure within 1s before the nadir
    batch_summary <- calculate_max_pressure_before_nadir(data, batch_summary, selected_sensor, num_rows)
    
    # Update plot with new nadir value and 1s max pressure
    sensor_summary <- batch_summary %>% filter(file == selected_sensor)
    plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 2, num_rows)
    print(plotly_plot)
    
    # Process 1s max nadir value. User can input new value if needed
    batch_summary <- process_max_p_1s_nadir(batch_summary, selected_sensor)
    
    # Calculate rate pressure change
    batch_summary <- calculate_rate_pressure_change(batch_summary, selected_sensor)
    
    #Use default acclimaiton pressure or set new
    max_acclim <- set_max_acclimation_pressure()
    
    # Calculate ratio pressure change and log ratio pressure change
    batch_summary <- calculate_ratio_pressure_change(batch_summary, selected_sensor, max_acclim)
    
    # Prompt user to continue with ROI selection
    continue_injection <- readline(prompt = "Enter Y when ready to continue with ROI selection: ")
    if (toupper(continue_injection) != "Y") {
      cat("BDS analysis tool ended.\n")
      return(invisible())
    }
    
    # Draw plot for ROI selection
    plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 3)
    print(plotly_plot)
    
    # Prompt user to enter ROI
    batch_summary <- prompt_for_injection_tailwater_times(data, batch_summary, selected_sensor, num_rows)
    
    # Assign passage points
    data <- assign_passage_points(data, batch_summary, selected_sensor, sample_rate)
    
    #Find max acceleration values within each ROI. ROI 1 and 2 use a range, RO1 3 and 4 use event time
    batch_summary <- max_acceleration_extract(data, batch_summary, selected_sensor)
    
    #Update sensor_summary and ensure NA columns are dropped from other sensors. 
    #May be redundant, but plot fucntion relises on it at the moment
    sensor_summary <- batch_summary %>%
      filter(file == selected_sensor) %>%
      select_if(~ any(!is.na(.)))
    
    # Draw plot using ROI 1 (either 6s, or user input)
    plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 4)
    print(plotly_plot)
    
    #Prompt user to continue with time series normalization
    continue_normal <- readline(prompt = "Enter Y when ready to continue with passage time series normalization: ")
    if (toupper(continue_normal) != "Y") {
      cat("BDS analysis tool ended.\n")
      return(invisible())
    }
    
    # Perform time series normalization on ROI 1 and update data
    data <- time_normalization(data, selected_sensor)
    
    # Draw plot with normalized ROI 1
    plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 5)
    print(plotly_plot)
    
    # Update batch summary with processed status
    batch_summary <- batch_summary %>%
      mutate(pres_processed = ifelse(file == selected_sensor, "Y", pres_processed),
             acc_processed = ifelse(file == selected_sensor, "Y", acc_processed))
  }
  
  #if the result was not true (e.g., user said Yes to bad sensor), then update batch_summary to say bad sensor
  
  if (!result) {
    batch_summary <- batch_summary %>%
      mutate(badsens = ifelse(file == selected_sensor, "Y", badsens))
  }
  
  # Update global environment to reassign batch_summary and data.
  update_global_environment(data, batch_summary, sample_rate)
  
  # BDS analysis complete. Prompt user to continue or end
  cat("BDS analysis complete for", selected_sensor, "\n")
  continue <- readline(prompt = "Continue with BDS analysis? (Y/N): ")
  if (toupper(continue) == "Y") {
    BDSAnalysisTool(batch_summary, data_250hz, data_100hz, data_100_imp, data_2000_hig, sample_rate)
  } else {
    cat("BDS Analysis tool ended\n")
  }
  # End of BDS analysis   
}



##4.1 Misc####

# find_acceleration_peaks <- function(data, batch_summary, selected_sensor, threshold = 49.03, min_gap = 12) {
#   # Set threshold to desired magnitude and set min_gap to time series (0.1s = 12)
#   data <- filter(data, long_id == selected_sensor)
#   peaks <- findpeaks(data$accmag, minpeakheight = threshold, minpeakdistance = min_gap)
#   
#   if (is.null(peaks)) {
#     cat("No acceleration peaks >= 49.03 m/s2 found\n")
#     return(batch_summary)
#   }
#   
#   peak_indices <- peaks[, 2]
#   peak_values <- data$accmag[peak_indices]
#   peak_times <- data$time[peak_indices]
#   peak_pres <- data$pres[peak_indices]
#   num_peaks <- length(peak_values)
#   
#   # Update batch_summary with peak values dynamically
#   for (i in seq_along(peak_values)) {
#     accmag_col <- paste0("accmag_", i)
#     time_col <- paste0("accmag_", i, "_t")
#     pres_col <- paste0("accmag_", i, "_p")
#     
#     if (!(accmag_col %in% colnames(batch_summary))) {
#       batch_summary[[accmag_col]] <- NA_real_
#     }
#     if (!(time_col %in% colnames(batch_summary))) {
#       batch_summary[[time_col]] <- NA_real_
#     }
#     if (!(pres_col %in% colnames(batch_summary))) {
#       batch_summary[[pres_col]] <- NA_real_
#     }
#     
#     batch_summary <- batch_summary %>%
#       mutate(
#         !!accmag_col := ifelse(file == selected_sensor, peak_values[i], .data[[accmag_col]]),
#         !!time_col := ifelse(file == selected_sensor, peak_times[i], .data[[time_col]]),
#         !!pres_col := ifelse(file == selected_sensor, peak_pres[i], .data[[pres_col]])
#       )
#   }
#   
#   cat(num_peaks, "acceleration peaks >= 49.03 m/s2 identified and recorded in batch summary\n")
#   
#   return(batch_summary)
# }



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

subset_by_long_id <- function(data) {
  unique_ids <- unique(data$long_id)
  
  for (id in unique_ids) {
    subset_data <- data %>% filter(long_id == id)
    subset_name <- as.character(id)
    assign(subset_name, subset_data, envir = .GlobalEnv)
    cat(paste("Created subset:", subset_name, "\n"))
  }
}

save_data <- function(batch_summary, data_100hz, data_250hz, data_100_imp, data_2000_hig) {
  # Check relevant data class (100, 250, or 100_imp)
  cat("Enter the sample rate (100, 250, 100_imp or 2000_hig\n")
  sample_rate <- readline(prompt = "Sample rate: ")
  if (!sample_rate %in% c("100", "250", "100_imp", "2000_imp")) {
    stop("Invalid sample rate. Please enter 100, 250, 100_imp or 2000_hig.")
  }
  
  # Check if pres_processed and acc_processed all = Y
  if (all(batch_summary$pres_processed == "Y") && all(batch_summary$acc_processed == "Y")) {
    message <- "Batch processing complete, proceed with saving? (Y/N): "
  } else {
    message <- "Caution: Batch processing incomplete, proceed with saving? (Y/N): "
  }
  
  proceed <- readline(prompt = message)
  if (toupper(proceed) != "Y") {
    cat("Saving process aborted.\n")
    return(invisible())
  }
  
  # Prompt user to enter batch name
  batch_name <- readline(prompt = "Enter batch name: ")
  
  cat("Saving batch summary and data file...\n")
  
  # Save batch_summary and data_ as .csv files
  output_dir <- "./R_output_files"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  batch_summary_file <- file.path(output_dir, paste0(batch_name, "_batch_summary.csv"))
  write.csv(batch_summary, batch_summary_file, row.names = FALSE)
  
  # Filter data_ by passage_point to remove any NA values
  if (sample_rate == "100") {
    data_to_save <- data_100hz %>% filter(!is.na(passage_point))
  } else if (sample_rate == "250") {
    data_to_save <- data_250hz %>% filter(!is.na(passage_point))
  } else if (sample_rate == "100_imp") {
    data_to_save <- data_100_imp %>% filter(!is.na(passage_point))
  } else if (sample_rate == "2000_hig") {
    data_to_save <- data_2000_hig %>% filter(!is.na(passage_point))
  } else {
    stop("Invalid sample rate.")
  }
  
  data_file <- file.path(output_dir, paste0(batch_name, "_data_", sample_rate, ".csv"))
  write.csv(data_to_save, data_file, row.names = FALSE)
  
  cat("Files saved successfully.\n")
}



assign_passage_points_manual <- function(sample_rate, batch_summary) {
  cat("Assigning passage points for sample rate:", sample_rate, "\n")
  
  # Determine which dataframe to use based on sample_rate
  data_name <- switch(sample_rate,
                      "100" = "data_100hz",
                      "250" = "data_250hz",
                      "100_imp" = "data_100_imp",
                      "2000_hig" = "data_2000_hig",
                      stop("Invalid sample_rate. Must be '100', '250', '100_imp', or '2000_hig'"))
  
  # Get the data from the global environment
  data <- get(data_name, envir = .GlobalEnv)
  
  # Ensure passage_point column exists and is initialized as NA
  if (!"passage_point" %in% colnames(data)) {
    data$passage_point <- NA_character_
  }
  
  # Process each file in batch_summary
  for (file in batch_summary$file) {
    # Extract injection and tailwater times from batch_summary
    file_summary <- batch_summary[batch_summary$file == file, ]
    t_injection_val <- file_summary$t_injection
    t_nadir_val <- file_summary$t_nadir
    t_tailwater_val <- file_summary$t_tailwater
    
    # Only process if both t_injection and t_tailwater are not 0.00
    if (t_injection_val != 0.00 && t_tailwater_val != 0.00) {
      cat("Processing file:", file, "\n")
      cat("Injection time:", t_injection_val, "\n")
      cat("Nadir time:", t_nadir_val, "\n")
      cat("Tailwater time:", t_tailwater_val, "\n")
      
      # Assign passage points based on time ranges
      data <- data %>%
        mutate(passage_point = case_when(
          long_id == file & time >= t_injection_val & time < t_nadir_val ~ "pre_nadir",
          long_id == file & time >= t_nadir_val & time <= t_tailwater_val ~ "post_nadir",
          TRUE ~ passage_point
        ))
      
      cat("Passage points assigned to", file, "\n")
    } else {
      cat("Skipping file:", file, "(t_injection or t_tailwater is 0.00)\n")
    }
  }
  
  # Convert passage_point to factor
  data$passage_point <- factor(data$passage_point, levels = c("pre_nadir", "post_nadir"))
  
  # Assign the modified data back to the global environment
  assign(data_name, data, envir = .GlobalEnv)
  
  cat("Passage points assignment completed and data updated in the global environment.\n")
}



time_normalization_manual <- function(sample_rate, batch_summary) {
  cat("Performing time normalization for sample rate:", sample_rate, "\n")
  
  # Determine which dataframe to use based on sample_rate
  data_name <- switch(sample_rate,
                      "100" = "data_100hz",
                      "250" = "data_250hz",
                      "100_imp" = "data_100_imp",
                      "2000_hig" = "data_2000_hig",
                      stop("Invalid sample_rate. Must be '100', '250', '100_imp', or '2000_hig'"))
  
  # Get the data from the global environment
  data <- get(data_name, envir = .GlobalEnv)
  
  # Ensure time_norm column exists and is initialized as 0
  if (!"time_norm" %in% colnames(data)) {
    data$time_norm <- 0
  }
  
  # Process each file in batch_summary
  for (file in batch_summary$file) {
    # Filter data for the current file
    data_file <- data[data$long_id == file, ]
    
    # Check if there are any pre_nadir and post_nadir points
    if (any(data_file$passage_point == 'pre_nadir', na.rm = TRUE) && 
        any(data_file$passage_point == 'post_nadir', na.rm = TRUE)) {
      
      cat("Processing file:", file, "\n")
      
      # Identify the start and end times for normalization
      start_time <- min(data_file$time[data_file$passage_point == 'pre_nadir'], na.rm = TRUE)
      mid_time <- min(data_file$time[data_file$passage_point == 'post_nadir'], na.rm = TRUE)
      end_time <- max(data_file$time[data_file$passage_point == 'post_nadir'], na.rm = TRUE)
      
      # Calculate normalized time values
      data <- data %>%
        mutate(time_norm = case_when(
          long_id == file & time <= start_time ~ 0,
          long_id == file & time >= end_time ~ 1,
          long_id == file & time > start_time & time < mid_time ~ 
            (time - start_time) / (mid_time - start_time) * 0.5,
          long_id == file & time >= mid_time & time <= end_time ~ 
            0.5 + (time - mid_time) / (end_time - mid_time) * 0.5,
          TRUE ~ time_norm
        ))
      
      cat("Time normalization completed for", file, "\n")
    } else {
      cat("Skipping file:", file, "(missing pre_nadir or post_nadir points)\n")
    }
  }
  
  # Assign the modified data back to the global environment
  assign(data_name, data, envir = .GlobalEnv)
  
  cat("Time normalization completed and data updated in the global environment.\n")
}
