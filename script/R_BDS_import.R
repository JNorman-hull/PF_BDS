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
    "acc_max[m/s2]", "acc_max[time]", "max_impact[g-force]")
  
  for (col in cols_to_numeric) {
    if (col %in% colnames(batch_summary)) {
      batch_summary[[col]] <- as.numeric(batch_summary[[col]])
    } else {
      batch_summary[[col]] <- NA_real_  # Initialize with NA if the column does not exist
    }
  }
  
  additional_cols <- c(
    "t_injection", "t_max_p_1s_nadir","t_nadir", "t_tailwater", "nadir",
    "max_p_1s_nadir", "rate_pc", "ratio_pc", "log_RPC")
  for (col in additional_cols) {
    batch_summary[[col]] <- as.numeric(0)
  }
  
  # Add 'processed' column as a factor with default value "N"
  batch_summary$pres_processed <- as.character("N")
  batch_summary$acc_processed <- as.character("N")
  
  cat("Batch summary loaded successfully.\n")
  return(batch_summary)
}

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

subset_by_long_id <- function(data) {
  unique_ids <- unique(data$long_id)
  
  for (id in unique_ids) {
    subset_data <- data %>% filter(long_id == id)
    subset_name <- as.character(id)
    assign(subset_name, subset_data, envir = .GlobalEnv)
    cat(paste("Created subset:", subset_name, "\n"))
  }
}

save_data <- function(batch_summary, data_100hz, data_250hz) {
  # Check relevant data class (100 or 250)
  cat("Enter the sample rate (100 or 250)\n")
  sample_rate <- readline(prompt = "Sample rate: ")
  if (!sample_rate %in% c("100", "250")) {
    stop("Invalid sample rate. Please enter either 100 or 250.")
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
  } else {
    stop("Invalid sample rate.")
  }
  
  data_file <- file.path(output_dir, paste0(batch_name, "_data_", sample_rate, "hz.csv"))
  write.csv(data_to_save, data_file, row.names = FALSE)
  
  cat("Files saved successfully.\n")
}

filter_batch_summary <- function(batch_summary, sample_rate) {
  if (sample_rate == 250) {
    data_class <- "250hz"
    num_rows <- 96
  } else if (sample_rate == 100) {
    data_class <- "100hz"
    num_rows <- 96
  } else {
    stop("Invalid sample rate")
  }
  
  #num_rows is 96 for both sensors as 96 rows = 1s. Only included for future redundancy 
  
  #modify this to allow for either to be Y but not both
  sensors <- batch_summary %>%
    filter(class == data_class & (pres_processed == "N" | acc_processed == "N")) %>%
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
  if (stage == 5 || stage == 6) {
    data <- filter(data, long_id == selected_sensor, !is.na(passage_point))
  } else {
    data <- filter(data, long_id == selected_sensor)
  }
  
  if (stage == 2) {
    nadir_time <- sensor_summary$t_nadir
    nadir_index <- which(data$time == nadir_time)
    row_start <- max(1, nadir_index - num_rows * 2)
    row_end <- min(nrow(data), nadir_index + num_rows * 2)
    data <- data[row_start:row_end, ]
  }
  
  # Base plot
  p <- plot_ly(data)
  
  if (stage != 3) {
    if (stage == 6) {
      p <- p %>%
        add_lines(x = ~time_norm, y = ~pres, name = 'Pressure', yaxis = 'y1', line = list(color = 'red', width = 1))
    } else {
      p <- p %>%
        add_lines(x = ~time, y = ~pres, name = 'Pressure', yaxis = 'y1', line = list(color = 'red', width = 1))
    }
  }
  
  if (stage >= 3) {
    if (stage == 6) {
      p <- p %>%
        add_lines(x = ~time_norm, y = ~accmag, name = 'Acceleration', yaxis = 'y2', line = list(color = 'blue', width = 1))
    } else {
      p <- p %>%
        add_lines(x = ~time, y = ~accmag, name = 'Acceleration', yaxis = 'y2', line = list(color = 'blue', width = 1))
    }
    
    # Add horizontal lines for stage 3
    if (stage == 3) {
      p <- p %>%
        layout(
          shapes = list(
            list(type = "line", 
                 x0 = min(data$time), x1 = max(data$time), xref = "x",
                 y0 = 49.03, y1 = 49.03, yref = "y2",
                 line = list(color = "black", dash = "dash")),
            list(type = "line", 
                 x0 = min(data$time), x1 = max(data$time), xref = "x",
                 y0 = 98.07, y1 = 98.07, yref = "y2",
                 line = list(color = "red", dash = "dash")),
            list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 1, y1 = 1, yref = "paper", line = list(color = "black", width = 1)),
            list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0, yref = "paper", line = list(color = "black", width = 1)),
            list(type = "line", x0 = 0, x1 = 0, xref = "paper", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1)),
            list(type = "line", x0 = 1, x1 = 1, xref = "paper", y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1))
          )
        ) %>%
        add_annotations(
          x = max(data$time) - (max(data$time) - min(data$time)) * 0.05,
          y = 49.03,
          text = 'G-force >= 5g',
          yref = "y2",
          showarrow = FALSE,
          xanchor = "right",
          font = list(color = 'black')
        ) %>%
        add_annotations(
          x = max(data$time) - (max(data$time) - min(data$time)) * 0.05,
          y = 98.07,
          text = 'G-force >= 10g',
          yref = "y2",
          showarrow = FALSE,
          xanchor = "right",
          font = list(color = 'red')
        )
    }
  }
  
  if (stage == 2) {
    nadir_index <- which(data$time == sensor_summary$t_nadir)
    time_start <- data$time[max(1, nadir_index - num_rows)]
    time_end <- data$time[min(nrow(data), nadir_index + num_rows)]
    vertical_lines <- list(
      list(type = "line", 
           x0 = time_start, x1 = time_start, xref = "x",
           y0 = 0, y1 = 1, yref = "paper",
           line = list(color = "black", dash = "dash")),
      list(type = "line", 
           x0 = time_end, x1 = time_end, xref = "x",
           y0 = 0, y1 = 1, yref = "paper",
           line = list(color = "black", dash = "dash"))
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
          list(
            x = (time_start + time_end) / 2,
            y = 0.8,
            text = '1s time window',
            showarrow = FALSE,
            xref = 'x',
            yref = 'paper',
            font = list(color = 'black', size = 12),
            xanchor = 'center',
            yanchor = 'middle'
          )
        )
      )
  }
  
  p <- p %>%
    layout(
      xaxis = list(
        title = "Time [s]",
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
  
  # Add dynamic points for pressure nadir based on the stage
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
  }
  
  if (stage == 2 || stage == 4 || stage == 5) {
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
  
  if (stage == 3 || stage == 4 || stage == 5) {
    peak_cols <- grep("^accmag_\\d+$", names(sensor_summary), value = TRUE)
    for (col in peak_cols) {
      peak_times <- sensor_summary[[paste0(col, "_t")]]
      peak_values <- sensor_summary[[col]]
      p <- p %>%
        add_markers(
          x = peak_times,
          y = peak_values,
          name = paste0(col),
          yaxis = 'y2',
          marker = list(color = 'orange', size = 10, line = list(color = 'black', width = 1))
        )
    }
  }
  
  p
}

find_acceleration_peaks <- function(data, batch_summary, selected_sensor, threshold = 49.03, min_gap = 12) {
  # Set threshold to desired magnitude and set min_gap to time series (0.1s = 12)
  data <- filter(data, long_id == selected_sensor)
  peaks <- findpeaks(data$accmag, minpeakheight = threshold, minpeakdistance = min_gap)
  
  if (is.null(peaks)) {
    cat("No acceleration peaks >= 49.03 m/s2 found\n")
    return(batch_summary)
  }
  
  peak_indices <- peaks[, 2]
  peak_values <- data$accmag[peak_indices]
  peak_times <- data$time[peak_indices]
  peak_pres <- data$pres[peak_indices]
  num_peaks <- length(peak_values)
  
  # Update batch_summary with peak values dynamically
  for (i in seq_along(peak_values)) {
    accmag_col <- paste0("accmag_", i)
    time_col <- paste0("accmag_", i, "_t")
    pres_col <- paste0("accmag_", i, "_p")
    
    if (!(accmag_col %in% colnames(batch_summary))) {
      batch_summary[[accmag_col]] <- NA_real_
    }
    if (!(time_col %in% colnames(batch_summary))) {
      batch_summary[[time_col]] <- NA_real_
    }
    if (!(pres_col %in% colnames(batch_summary))) {
      batch_summary[[pres_col]] <- NA_real_
    }
    
    batch_summary <- batch_summary %>%
      mutate(
        !!accmag_col := ifelse(file == selected_sensor, peak_values[i], .data[[accmag_col]]),
        !!time_col := ifelse(file == selected_sensor, peak_times[i], .data[[time_col]]),
        !!pres_col := ifelse(file == selected_sensor, peak_pres[i], .data[[pres_col]])
      )
  }
  
  cat(num_peaks, "acceleration peaks >= 49.03 m/s2 identified and recorded in batch summary\n")
  
  return(batch_summary)
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

prompt_for_injection_tailwater_times <- function(batch_summary, selected_sensor) {
  cat("Visually inspect the plot to determine injection and tailwater points. Enter injection and tailwater time\n")
  t_injection_val <- as.numeric(readline(prompt = "Enter injection time value: "))
  t_tailwater_val <- as.numeric(readline(prompt = "Enter tailwater time value: "))
  if (is.na(t_injection_val) || is.na(t_tailwater_val)) {
    stop("Invalid injection or tailwater time")
  }
  
  batch_summary <- batch_summary %>%
    mutate(
      t_injection = ifelse(file == selected_sensor, t_injection_val, t_injection),
      t_tailwater = ifelse(file == selected_sensor, t_tailwater_val, t_tailwater)
    )
  
  cat("Injection and tailwater times recorded in batch summary\n")
  
  batch_summary
}

assign_passage_points <- function(data, batch_summary, selected_sensor, sample_rate) {
  cat("Assigning passage points for sensor:", selected_sensor, "\n")
  
  # Extract injection, nadir, and tailwater times from batch_summary
  t_injection_val <- batch_summary %>% filter(file == selected_sensor) %>% pull(t_injection)
  t_nadir_val <- batch_summary %>% filter(file == selected_sensor) %>% pull(t_nadir)
  t_tailwater_val <- batch_summary %>% filter(file == selected_sensor) %>% pull(t_tailwater)
  
  cat("Injection time:", t_injection_val, "\n")
  cat("Nadir time:", t_nadir_val, "\n")
  cat("Tailwater time:", t_tailwater_val, "\n")
  

  # Ensure passage_point column exists and is initialized as NA
  if (!"passage_point" %in% colnames(data)) {
    data$passage_point <- NA_character_
  }
  
  # Filter the data for the selected sensor
  data_sensor <- data %>% filter(long_id == selected_sensor)
  
  # Assign passage points based on time ranges
  data_sensor <- data_sensor %>%
    mutate(passage_point_new = case_when(
      time >= t_injection_val & time < t_nadir_val ~ "pre_nadir",
      time >= t_nadir_val & time <= t_tailwater_val ~ "post_nadir",
      TRUE ~ NA_character_
    ))
  
  # Convert passage_point_new to factor
  data_sensor <- data_sensor %>%
    mutate(passage_point_new = factor(passage_point_new, levels = c("pre_nadir", "post_nadir")))
  
  cat("Passage points assigned to", selected_sensor, "\n")
  
  # Merge the passage points back into the original data
  data <- data %>%
    left_join(select(data_sensor, time, long_id, passage_point_new), by = c("time", "long_id")) %>%
    mutate(passage_point = coalesce(passage_point_new, passage_point)) %>%
    select(-passage_point_new)
  
  return(data)
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


#BDS analysis tool prototype
BDSAnalysisTool <- function(batch_summary, data_250hz, data_100hz, sample_rate) {

  filter_result <- filter_batch_summary(batch_summary, sample_rate)
  if (!is.null(filter_result$message)) {
    cat(filter_result$message)
    cat("BDS Analysis tool ended.\n")
    return(invisible())
  }
  
  data <- if (sample_rate == 250) data_250hz else data_100hz
  
  #Create list of sensors available for processing. 
  #Add num_rows for calculating max pressure <1s nadir
  sensors <- filter_result$sensors
  num_rows <- filter_result$num_rows
  #Prompt user to select a sensor
  selected_sensor <- prompt_user_to_select_sensor(sensors)
  #Create sensor summary
  sensor_summary <- batch_summary %>% filter(file == selected_sensor)
  time_range <- max(data$time) - min(data$time)

  #Draw plot for nadir value check
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 1)
  print(plotly_plot)
  
  #Process nadir value. User can input new value if needed
  batch_summary <- process_nadir_value(batch_summary, selected_sensor)
  
  #Calculate max pressure within 1s before the nadir
  batch_summary <- calculate_max_pressure_before_nadir(data, batch_summary, selected_sensor, num_rows)

  #Update plot with new nadir value and 1s max pressure
  sensor_summary <- batch_summary %>% filter(file == selected_sensor)
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 2, num_rows)
  print(plotly_plot)
  
  #Process 1s max nadir value. User can input new value if needed
  batch_summary <- process_max_p_1s_nadir(batch_summary, selected_sensor)
  
  #Calculate rate pressure change
  batch_summary <- calculate_rate_pressure_change(batch_summary, selected_sensor)
  
  cat("Enter maximum acclimation pressure (e.g., atmospheric pressure or static pressure within system)\n")
  max_acclim <- as.numeric(readline(prompt = "Maximum acclimation pressure: "))
  if (is.na(max_acclim)) {
    stop("Invalid maximum acclimation pressure")
  }
  
  #Calculate ratio pressure change and log ratio pressure change
  batch_summary <- calculate_ratio_pressure_change(batch_summary, selected_sensor, max_acclim)
  
  #Calculate acceleration peaks
  batch_summary <- find_acceleration_peaks(data, batch_summary, selected_sensor)

  #Draw plot for acceleration magnitude 
  #First update sensor_summary and ensure accmag_ NA columns are dropped from other sensors
  sensor_summary <- batch_summary %>%
    filter(file == selected_sensor) %>%
    select_if(~ any(!is.na(.)))
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 3)
  print(plotly_plot)
  
  #Prompt user to continue with injection and tailwater selection
  continue_injection <- readline(prompt = "Enter Y when ready to continue with injection and tailwater selection: ")
  if (toupper(continue_injection) != "Y") {
    cat("BDS analysis tool ended.\n")
    return(invisible())
  }

  #Draw plot for ROI selection
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 4)
  print(plotly_plot)
  
  #Prompt user to enter injection and tailwater time
  batch_summary <- prompt_for_injection_tailwater_times(batch_summary, selected_sensor)
  
  #Assign passage points
  data <- assign_passage_points(data, batch_summary, selected_sensor, sample_rate)
  
  #Draw plot with updated passage points
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 5)
  print(plotly_plot)
  
  #Prompt user to continue with time series normalization
  continue_normal<- readline(prompt = "Enter Y when ready to continue with passage time series normalization: ")
  if (toupper(continue_normal) != "Y") {
    cat("BDS analysis tool ended.\n")
    return(invisible())
  }

  #Perform time series normalization and update data
  data <- time_normalization(data, selected_sensor)

  #Draw plot with normalized passage
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 6)
  print(plotly_plot)

  #Update batch summary with processed status
  batch_summary <- batch_summary %>%
    mutate(pres_processed = ifelse(file == selected_sensor, "Y", pres_processed),
           acc_processed = ifelse(file == selected_sensor, "Y", acc_processed))

  #Update global environment
  if (sample_rate == 250) {
    data_250hz <- data
    assign("data_250hz", data_250hz, envir = .GlobalEnv)
  } else if (sample_rate == 100) {
    data_100hz <- data
    assign("data_100hz", data_100hz, envir = .GlobalEnv)
  }
  assign("batch_summary", batch_summary, envir = .GlobalEnv)
  
  #BDS analysis complete. Prompt user to continue or end
  cat("BDS analysis complete for", selected_sensor, "\n")
    continue <- readline(prompt = "Continue with BDS analysis? (Y/N): ")
  if (toupper(continue) == "Y") {
    plot_selected_sensor(batch_summary, data_250hz, data_100hz, sample_rate)
  } else {
    cat("BDS Analysis tool ended\n")
  }
  #End of BDS analysis   
}
