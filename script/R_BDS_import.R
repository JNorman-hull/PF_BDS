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

py_install("pillow")
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

add_video_metrics <- function(batch_meta_dir) {
  # Load video metrics CSV
  video_metrics_path <- file.path(batch_meta_dir, "video_metrics.csv")
  if (!file.exists(video_metrics_path)) {
    stop("video_metrics.csv not found in batch_meta_dir")
  }
  
  # Define the columns to extract
  cols_to_extract <- c(
    "treatment", "inj_pos", "frame_in", "frame_out", "pl_frame", "c_1_f", "c_2_f",
    "sens_wc_pos_ent", "sens_wc_pos_blade","sens_wc_pos_exit",
    "clear_passage", "rotate_hub", "pres_force", "passage_loc",
    "centre_hub_contact", "blade_contact", "n_contact",
    "c_1_bl", "c_1_sl", "c_1_bv", "c_1_type",
    "c_2_bl", "c_2_sl", "c_2_bv", "c_2_type"
  )
  
  #Removed cols not needed for 400_1, 500_1
  # Kept frame info incase needed
  # removed orientation as not likely to analyse  "sens_ec_orien", 
  #sens_velocity", "pre_swirl", "c_3_bl", "c_3_sl", "c_3_bv", "c_3_f"
  # Read video metrics and select columns
  video_metrics <- read_csv(video_metrics_path, show_col_types = FALSE) %>%
    select(sens_file, all_of(cols_to_extract))
  
  # Update batch_summary in global environment
  batch_summary <<- batch_summary %>%
    mutate(join_key = gsub("_imp$", "", file)) %>%  # Remove _imp for joining
    left_join(
      video_metrics,
      by = c("join_key" = "sens_file")
    ) %>%
    select(-join_key)  # Remove temporary join key
  
  cat("Video metrics added to batch_summary\n")
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
  sensor_labels <- sensors %>%
    left_join(batch_summary %>% select(file, treatment, clear_passage, blade_contact, centre_hub_contact), 
              by = c("file")) %>%
    mutate(display_name = paste0(
      file,
      " Treatment = ", coalesce(treatment, "Unknown"),
      case_when(
        clear_passage == "Y" ~ " Clear sensor passage",
        TRUE ~ paste0(
          case_when(
            blade_contact == "Y" ~ " Blade strike",
            blade_contact == "N" ~ " No blade strike",
            TRUE ~ ""
          ),
          case_when(
            centre_hub_contact == "Y" ~ " Centre hub contact",
            centre_hub_contact == "N" ~ " No centre hub contact",
            TRUE ~ ""
          )
        )
      )
    ))
  
  # Rest of the function remains the same
  cat("Select a sensor dataset to begin BDS analysis:\n")
  for (i in seq_along(sensor_labels$display_name)) {
    cat(i, ": ", sensor_labels$display_name[i], "\n", sep = "")
  }
  
  selected_sensor_index <- as.integer(readline(prompt = "Enter the number corresponding to the sensor: "))
  if (is.na(selected_sensor_index) || selected_sensor_index < 1 || selected_sensor_index > nrow(sensor_labels)) {
    stop("Invalid sensor selection")
  }
  
  sensor_labels$file[selected_sensor_index]
}

create_combined_plot <- function(data, sensor_summary, selected_sensor, stage, num_rows) {
  # Filter data for the selected sensor
  data <- filter(data, long_id == selected_sensor)
  
  # Get sensor information for title
  sensor_title <- batch_summary %>%
    filter(file == selected_sensor) %>%
    mutate(display_title = paste0(
      file,
      " Treatment = ", coalesce(treatment, "Unknown"),
      case_when(
        clear_passage == "Y" ~ " Clear sensor passage",
        TRUE ~ paste0(
          case_when(
            blade_contact == "Y" ~ " Blade strike",
            blade_contact == "N" ~ " No blade strike",
            TRUE ~ ""
          ),
          case_when(
            centre_hub_contact == "Y" ~ " Centre hub contact",
            centre_hub_contact == "N" ~ " No centre hub contact",
            TRUE ~ ""
          )
        )
      )
    )) %>%
    pull(display_title)
  
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
  
  # Base plot with title
  p <- plot_ly(data) %>%
    layout(
      title = list(
        text = sensor_title,
        x = 0.5,
        xanchor = 'center',
        yanchor = 'top',
        y = 0.95
      )
    )
  
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
    # Define colors for each ROI
    roi_colors <- c(
      "inj_to_pre_nadir" = "lightblue",
      "pre_nadir_event" = "lightgreen",
      "nadir_event" = "pink",
      "post_nadir_event" = "lightgreen",
      "post_nadir_to_rec" = "lightblue"
    )
    
    # Add shaded regions for each ROI
    for (roi_name in names(roi_colors)) {
      roi_data <- data %>% filter(roi == roi_name)
      if (nrow(roi_data) > 0) {
        p <- p %>% add_ribbons(
          data = roi_data,
          x = ~time,
          ymin = min(data$pres),
          ymax = max(data$pres),
          fillcolor = roi_colors[roi_name],
          opacity = 0.3,
          line = list(width = 0),
          name = roi_name,
          showlegend = TRUE
        )
      }
    }
    
    p <- p %>%
      add_lines(x = ~time, y = ~pres, name = 'Pressure', yaxis = 'y1', line = list(color = 'red', width = 1)) %>%
      add_lines(x = ~time, y = ~accmag, name = 'Acceleration', yaxis = 'y2', line = list(color = 'blue', width = 1))
    
    #acceleration markers
    accmag_cols <- names(sensor_summary)[grep("^accmag_\\d+$", names(sensor_summary))]
    for (accmag_col in accmag_cols) {
      i <- as.numeric(sub("accmag_", "", accmag_col))
      time_col <- paste0(accmag_col, "_t")
      
      if (time_col %in% names(sensor_summary) &&
          !is.na(sensor_summary[[accmag_col]]) && !is.na(sensor_summary[[time_col]])) {
        p <- p %>%
          add_markers(
            x = sensor_summary[[time_col]],
            y = sensor_summary[[accmag_col]],
            name = paste0("Acc Peak ", i),
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

process_nadir_value <- function(batch_summary, selected_sensor, data) {
  is_nadir_correct <- readline(prompt = "Is nadir value correct? (Y/N): ")
  
  if (toupper(is_nadir_correct) == "Y") {
    batch_summary <- batch_summary %>%
      mutate(
        nadir = ifelse(file == selected_sensor, `pres_min[mbar]`, nadir),
        t_nadir = ifelse(file == selected_sensor, `pres_min[time]`, t_nadir)
      )
    cat("Nadir values recorded in batch summary\n")
  } else if (toupper(is_nadir_correct) == "N") {
    new_nadir_time <- as.numeric(readline(prompt = "Enter new nadir time: "))
    if (is.na(new_nadir_time)) {
      stop("Invalid time value")
    }
    
    # Get corresponding pressure value from data
    new_nadir_value <- data %>%
      filter(long_id == selected_sensor, time == new_nadir_time) %>%
      pull(pres)
    
    if (length(new_nadir_value) == 0) {
      stop("No pressure value found for the specified time")
    }
    
    batch_summary <- batch_summary %>%
      mutate(
        nadir = ifelse(file == selected_sensor, new_nadir_value, nadir),
        t_nadir = ifelse(file == selected_sensor, new_nadir_time, t_nadir)
      )
    cat(sprintf("New nadir time: %.2f, corresponding pressure: %.1f mbar\n", 
                new_nadir_time, new_nadir_value))
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

# process_max_p_1s_nadir <- function(batch_summary, selected_sensor) {
#   # Retrieve the current max_p_1s_nadir and t_max_p_1s_nadir values
#   sensor_summary <- batch_summary %>% filter(file == selected_sensor)
#   current_max_p_1s_nadir <- sensor_summary$max_p_1s_nadir
#   current_t_max_p_1s_nadir <- sensor_summary$t_max_p_1s_nadir
#   
#   # Display current values and prompt user for validation
#   cat("Current max pressure within 1s prior to nadir value:", current_max_p_1s_nadir, "\n")
#   cat("Current time of max pressure within 1s prior to nadir:", current_t_max_p_1s_nadir, "\n")
#   is_correct <- readline(prompt = "Is max pressure within 1s prior to nadir value correct? (Y/N): ")
#   
#   if (toupper(is_correct) == "N") {
#     new_max_p_1s_nadir <- as.numeric(readline(prompt = "Enter new max pressure within 1s prior to nadir value: "))
#     new_t_max_p_1s_nadir <- as.numeric(readline(prompt = "Enter new time of max pressure within 1s prior to nadir: "))
#     
#     # Update batch_summary with new values
#     batch_summary <- batch_summary %>%
#       mutate(
#         max_p_1s_nadir = ifelse(file == selected_sensor, new_max_p_1s_nadir, max_p_1s_nadir),
#         t_max_p_1s_nadir = ifelse(file == selected_sensor, new_t_max_p_1s_nadir, t_max_p_1s_nadir)
#       )
#     cat("New max pressure within 1s prior to nadir and time recorded in batch summary\n")
#   } else if (toupper(is_correct) == "Y") {
#     cat("Max pressure within 1s prior to nadir value and time remain unchanged\n")
#   } else {
#     cat("Invalid input. Please enter 'Y' or 'N'.\n")
#   }
#   
#   return(batch_summary)
# }

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

# set_max_acclimation_pressure <- function(default_pressure = 1000) {
#   cat("Default maximum acclimation pressure is set to", default_pressure, "mbar.\n")
#   change_max_acclim <- readline(prompt = "Do you want to change the maximum acclimation pressure? (Y/N): ")
#   
#   if (toupper(change_max_acclim) == "Y") {
#     max_acclim <- as.numeric(readline(prompt = "Enter new maximum acclimation pressure: "))
#     while (is.na(max_acclim)) {
#       cat("Invalid input. Please enter a numeric value.\n")
#       max_acclim <- as.numeric(readline(prompt = "Enter new maximum acclimation pressure: "))
#     }
#   } else {
#     max_acclim <- default_pressure
#   }
#   
#   cat("Maximum acclimation pressure set to", max_acclim, "mbar.\n")
#   return(max_acclim)
# }

set_max_acclimation_pressure <- function(default_pressure = 1000) {
  cat("Default maximum acclimation pressure is set to", default_pressure, "mbar.\n")
  max_acclim <- default_pressure
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

#Modify times to reflect increased passage duration for 50% data

prompt_for_injection_tailwater_times <- function(data, batch_summary, selected_sensor, num_rows) {
  sensor_data <- data %>% filter(long_id == selected_sensor)
  t_nadir_val <- batch_summary %>% filter(file == selected_sensor) %>% pull(t_nadir)
  nadir_index <- which(sensor_data$time == t_nadir_val)
  
  if (length(nadir_index) == 0) {
    stop("Nadir time not found in the data")
  }
  
  # Create ROI columns if they don't exist
  roi_columns <- paste0(c("t_start_roi", "t_end_roi"), rep(1:6, each = 2))
  
  for (col in roi_columns) {
    if (!(col %in% colnames(batch_summary))) {
      batch_summary[[col]] <- NA_real_
    }
  }
  
  for (roi in 1:6) {
    cat(paste("\nROI", roi, ":"))
    
    if (roi == 1) {
      cat(" (Main overall sensor passage ROI)\n")
      row_start <- max(1, nadir_index - num_rows * 2)
      row_end <- min(nrow(sensor_data), nadir_index + num_rows * 2)
      t_start_val <- sensor_data$time[row_start]
      t_end_val <- sensor_data$time[row_end]
      cat("Default start time for PF (nadir - 2s):", t_start_val, "\n")
      cat("Default end time for PF (nadir + 2s):", t_end_val, "\n")
      
    } else if (roi == 2) {
      cat(" (Nadir event ROI)\n")
      row_start <- max(1, nadir_index - 20)
      row_end <- min(nrow(sensor_data), nadir_index + 20)
      t_start_val <- sensor_data$time[row_start]
      t_end_val <- sensor_data$time[row_end]
      cat("Default nadir event start time for PF (nadir - 0.2s):", t_start_val, "\n")
      cat("Default nadir event end time for PF (nadir + 0.2s):", t_end_val, "\n")
      
    } else if (roi == 3) {
      cat(" (Pre-nadir ROI)\n")
      t_start_val <- sensor_data$time[max(1, which(sensor_data$time == batch_summary$t_start_roi2[batch_summary$file == selected_sensor]) - num_rows / 1.92)]
      t_end_val <- sensor_data$time[which(sensor_data$time == batch_summary$t_start_roi2[batch_summary$file == selected_sensor])]
      cat("Default start time for PF (nadir ROI start - 0.5s):", t_start_val, "\n")
      cat("End time:", t_end_val, "\n")
      
    } else if (roi == 4) {
      cat(" (Post-nadir ROI)\n")
      t_start_val <- sensor_data$time[which(sensor_data$time == batch_summary$t_end_roi2[batch_summary$file == selected_sensor])]
      t_end_val <- sensor_data$time[min(nrow(sensor_data), which(sensor_data$time == batch_summary$t_end_roi2[batch_summary$file == selected_sensor]) + num_rows / 1.92)]
      cat("Start time:", t_start_val, "\n")
      cat("Default end time for PF (nadir ROI end + 0.5s):", t_end_val, "\n")
      
    } else if (roi == 5) {
      t_start_val <- batch_summary$t_start_roi1[batch_summary$file == selected_sensor]
      t_end_val <- batch_summary$t_start_roi3[batch_summary$file == selected_sensor]
      cat(" (inj_to_pre_nadir)\n")
      cat("Start time:", t_start_val, "\n")
      cat("End time:", t_end_val, "\n")
    } else if (roi == 6) {
      t_start_val <- batch_summary$t_end_roi4[batch_summary$file == selected_sensor]
      t_end_val <- batch_summary$t_end_roi1[batch_summary$file == selected_sensor]
      cat(" (post_nadir_to_rec)\n")
      cat("Start time:", t_start_val, "\n")
      cat("End time:", t_end_val, "\n")
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
  cat("\nAssigning passage points for sensor:", selected_sensor, "\n")
  
  # First ensure the columns exist
  if(!"passage_point" %in% names(data)) data$passage_point <- NA_character_
  if(!"roi" %in% names(data)) data$roi <- NA_character_
  
  # Filter the data for the selected sensor
  data_sensor <- data %>% filter(long_id == selected_sensor)
  
  # Get all ROI times
  roi_times <- batch_summary %>%
    filter(file == selected_sensor) %>%
    select(starts_with("t_start_roi"), starts_with("t_end_roi"), t_nadir)
  
  data_sensor <- data_sensor %>%
    mutate(
      passage_point = case_when(
        time >= roi_times$t_start_roi1 & time < roi_times$t_nadir ~ "pre_nadir",
        time >= roi_times$t_nadir & time <= roi_times$t_end_roi1 ~ "post_nadir",
        TRUE ~ NA_character_
      ),
      roi = case_when(
        time >= roi_times$t_start_roi5 & time < roi_times$t_end_roi5 ~ "inj_to_pre_nadir",
        time >= roi_times$t_start_roi3 & time < roi_times$t_end_roi3 ~ "pre_nadir_event",
        time >= roi_times$t_start_roi2 & time < roi_times$t_end_roi2 ~ "nadir_event",
        time >= roi_times$t_start_roi4 & time < roi_times$t_end_roi4 ~ "post_nadir_event",
        time >= roi_times$t_start_roi6 & time <= roi_times$t_end_roi6 ~ "post_nadir_to_rec",
        TRUE ~ NA_character_
      )
    )
  
  cat("Passage points and ROIs assigned to", selected_sensor, "\n")
  
  # Update the data for the selected sensor
  data <- data %>%
    rows_update(data_sensor, by = c("time", "long_id"))
  
  return(data)
}


max_acceleration_extract <- function(data, batch_summary, selected_sensor) {
  # extract the max acceleration value in each ROI
  # Filter data for the selected sensor
  sensor_data <- data %>% filter(long_id == selected_sensor)
  
  # Create columns if they don't exist
  new_columns <- c(
    paste0("max_accmag_roi", 1:6),
    paste0("max_accmag_time_roi", 1:6),
    "nadir_acceleration"  # Add new column
  )
  for (col in new_columns) {
    if (!(col %in% colnames(batch_summary))) {
      batch_summary[[col]] <- NA_real_
    }
  }
  
  # Extract max acceleration
  for (roi in 1:6) {
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
  
  # Extract acceleration at nadir
  t_nadir <- batch_summary %>% 
    filter(file == selected_sensor) %>% 
    pull(t_nadir)
  
  nadir_acceleration <- sensor_data %>%
    filter(time == t_nadir) %>%
    pull(accmag)
  
  cat("Debug: t_nadir:", t_nadir, "\n")
  cat("Debug: nadir_acceleration:", nadir_acceleration, "\n")
  
  # Update batch_summary with nadir acceleration 
  batch_summary[batch_summary$file == selected_sensor, "nadir_acceleration"] <- nadir_acceleration
  
  # Print results to console
  cat("\nMaximum acceleration values:\n")
  cat("ROI 1 (Main overall sensor passage) max acceleration:", batch_summary$max_accmag_roi1[batch_summary$file == selected_sensor], "\n")
  cat("ROI 2 (nadir event) max acceleration:", batch_summary$max_accmag_roi2[batch_summary$file == selected_sensor], "\n")
  cat("ROI 3 (pre-nadir ROI) acceleration:", batch_summary$max_accmag_roi3[batch_summary$file == selected_sensor], "\n")
  cat("ROI 4 (post-nadir ROI) acceleration:", batch_summary$max_accmag_roi4[batch_summary$file == selected_sensor], "\n")
  cat("ROI 5 (inection to pre-nadir ROI) acceleration:", batch_summary$max_accmag_roi5[batch_summary$file == selected_sensor], "\n")
  cat("ROI 6 (post-nadir to recovery ROI) acceleration:", batch_summary$max_accmag_roi6[batch_summary$file == selected_sensor], "\n")
  cat("Acceleration at nadir:", nadir_acceleration, "\n")  # Add print statement for nadir acceleration
  
  return(batch_summary)
}

max_hig_extract <- function(data_2000_hig, batch_summary, selected_sensor) {
  # Convert selected_sensor ID from _imp to _hig
  hig_sensor <- gsub("_imp$", "_hig", selected_sensor)
  
  # Create columns if they don't exist
  new_columns <- c(
    paste0("max_hig_roi", 1:6),
    paste0("max_hig_time_roi", 1:6)
  )
  for (col in new_columns) {
    if (!(col %in% colnames(batch_summary))) {
      batch_summary[[col]] <- NA_real_
    }
  }
  
  # Check if HIG data exists for this sensor
  if (!hig_sensor %in% unique(data_2000_hig$long_id)) {
    cat("\nNo HIG data available for this sensor\n")
    batch_summary <- batch_summary %>%
      mutate(
        across(starts_with("max_hig_roi"), ~ifelse(file == selected_sensor, NA_real_, .)),
        across(starts_with("max_hig_time_roi"), ~ifelse(file == selected_sensor, NA_real_, .))
      )
    return(batch_summary)
  }
  
  # Extract max HIG for each ROI
  for (roi in 1:6) {
    t_start <- batch_summary %>% 
      filter(file == selected_sensor) %>% 
      pull(!!paste0("t_start_roi", roi))
    t_end <- batch_summary %>% 
      filter(file == selected_sensor) %>% 
      pull(!!paste0("t_end_roi", roi))
    
    roi_data <- data_2000_hig %>% 
      filter(long_id == hig_sensor,
             time >= t_start & time <= t_end)
    
    max_hig <- if(nrow(roi_data) > 0) max(roi_data$HIGaccmag, na.rm = TRUE) else NA_real_
    max_hig_time <- if(nrow(roi_data) > 0) roi_data$time[which.max(roi_data$HIGaccmag)] else NA_real_
    
    batch_summary <- batch_summary %>%
      mutate(
        !!paste0("max_hig_roi", roi) := ifelse(file == selected_sensor, max_hig, !!sym(paste0("max_hig_roi", roi))),
        !!paste0("max_hig_time_roi", roi) := ifelse(file == selected_sensor, max_hig_time, !!sym(paste0("max_hig_time_roi", roi)))
      )
  }
  
  cat("\nMaximum HIG acceleration values:\n")
  for(roi in 1:6) {
    val <- batch_summary[[paste0("max_hig_roi", roi)]][batch_summary$file == selected_sensor]
    cat(sprintf("ROI %d max HIG: %s\n", 
                roi, 
                ifelse(is.na(val), "NA", sprintf("%f", val))))
  }
  
  return(batch_summary)
}

nadir_acceleration_distance <- function(data, batch_summary, selected_sensor) {
  # Extract relevant values from batch_summary
  sensor_summary <- batch_summary[batch_summary$file == selected_sensor, ]
  t_nadir <- sensor_summary$t_nadir
  max_accmag_time_roi2 <- sensor_summary$max_accmag_time_roi2
  
  # Find the row indices for nadir and max acceleration
  nadir_index <- which.min(abs(data$time - t_nadir))
  max_acc_index <- which.min(abs(data$time - max_accmag_time_roi2))
  
  # Calculate the row difference
  row_diff <- max_acc_index - nadir_index
  
  # Calculate the time difference using fixed time step
  fixed_time_step <- 0.01  # 0.01 seconds per row
  max_nadir_acc_dist <- row_diff * fixed_time_step
  
  # Determine the position relative to nadir
  if (row_diff < 0) {
    max_nadir_acc_position <- "before_nadir"
  } else if (row_diff > 0) {
    max_nadir_acc_position <- "after_nadir"
  } else {
    max_nadir_acc_position <- "on_nadir"
  }
  
  # Update batch_summary with new values
  batch_summary[batch_summary$file == selected_sensor, "max_nadir_acc_dist"] <- max_nadir_acc_dist
  batch_summary[batch_summary$file == selected_sensor, "max_nadir_acc_position"] <- max_nadir_acc_position
  
  # Print message to console
  cat("\nnadir ROI: Max acceleration time differnce from nadir")
  cat(sprintf("\nMaximum acceleration in nadir ROI occurred %.3f seconds (%d rows) %s.\n", 
              max_nadir_acc_dist, abs(row_diff), max_nadir_acc_position))
  
  return(batch_summary)
}

acceleration_rms <- function(data, batch_summary, selected_sensor) {
  # Get ROI times for nadir event
  roi_start <- batch_summary$t_start_roi2[batch_summary$file == selected_sensor]
  roi_end <- batch_summary$t_end_roi2[batch_summary$file == selected_sensor]
  
  # Filter data and calculate RMS
  roi_data <- data %>%
    filter(long_id == selected_sensor,
           time >= roi_start,
           time <= roi_end)
  
  # Calculate RMS acceleration
  rms_acc <- roi_data %>%
    summarise(acc_rms = sqrt(mean(accmag^2))) %>%
    pull(acc_rms)
  
  # Add RMS column if it doesn't exist
  if(!("acc_rms" %in% names(batch_summary))) {
    batch_summary$acc_rms <- NA_real_
  }
  
  # Update batch_summary
  batch_summary$acc_rms[batch_summary$file == selected_sensor] <- rms_acc

  cat("RMS acceleration for nadir event ROI:", round(rms_acc, 2), "m/s²\n")
  
  return(batch_summary)
}

find_acceleration_peaks <- function(data, batch_summary, selected_sensor, 
                                    peak = 98.1, peak_gap = 5, drop = 35, drop_gap = 1,
                                    group_window_multiplier = 3) {
  # Filter data within t_start_roi1 to t_end_roi1 (overall passage)
  roi_start <- batch_summary$t_start_roi1[batch_summary$file == selected_sensor]
  roi_end <- batch_summary$t_end_roi1[batch_summary$file == selected_sensor]
  data <- filter(data, long_id == selected_sensor, time >= roi_start, time <= roi_end)
  
  #cat("\nDebug: Processing sensor:", selected_sensor, "\n")
  #cat("Debug: ROI start:", roi_start, "ROI end:", roi_end, "\n\n")
  
  # Find initial peaks
  peaks <- findpeaks(data$accmag, minpeakheight = peak, minpeakdistance = peak_gap)
  
  if (is.null(peaks) || nrow(peaks) == 0) {
    cat("No acceleration peaks meeting the criteria were found\n")
    return(batch_summary)
  }
  
  peak_indices <- peaks[, 2]
  peak_values <- data$accmag[peak_indices]
  peak_times <- data$time[peak_indices]
  
  #cat("Debug: Initial peaks found:", length(peak_indices), "\n")
  #cat("Debug: Peak times:", paste(peak_times, collapse = ", "), "\n")
  #cat("Debug: Peak values:", paste(peak_values, collapse = ", "), "\n\n")
  
  # Sort peaks by time
  sorted_order <- order(peak_times)
  peak_indices <- peak_indices[sorted_order]
  peak_values <- peak_values[sorted_order]
  peak_times <- peak_times[sorted_order]
  
  #cat("Debug: Sorted peak times:", paste(peak_times, collapse = ", "), "\n")
  #cat("Debug: Sorted peak values:", paste(peak_values, collapse = ", "), "\n\n")
  
  valid_peaks <- numeric()
  
  if (length(peak_indices) > 1) {
    group_window <- peak_gap * group_window_multiplier
    #cat("Debug: New calculated group window:", group_window, "\n")
    
    groups <- cumsum(c(1, diff(peak_indices) > group_window))
    #cat("Debug: Group assignments:", paste(groups, collapse = ", "), "\n\n")
    
    for (g in unique(groups)) {
      group_indices <- peak_indices[groups == g]
      group_values <- peak_values[groups == g]
      group_times <- peak_times[groups == g]
      
      #cat("\nDebug: Processing group", g, "\n")
      #cat("Debug: Group indices:", paste(group_indices, collapse = ", "), "\n")
      #cat("Debug: Group values:", paste(group_values, collapse = ", "), "\n")
      #cat("Debug: Group times:", paste(group_times, collapse = ", "), "\n")
      
      if (length(group_indices) == 1) {
        valid_peaks <- c(valid_peaks, group_indices)
        #cat("Debug: Single peak in group, keeping it:", group_values, "\n")
      } else {
        # Check for drops between consecutive peaks
        drops <- sapply(1:(length(group_indices)-1), function(i) {
          between_peaks <- data$accmag[(group_indices[i]+1):(group_indices[i+1]-1)]
          drop_runs <- rle(between_peaks <= drop)
          any(drop_runs$lengths[drop_runs$values] >= drop_gap)
        })
        
        #cat("Debug: Drops detected at positions:", paste(which(drops), collapse = ", "), "\n")
        #cat("Debug: Drop values:", paste(drops, collapse = ", "), "\n")
        
        if (any(drops)) {
          # Find highest peak before first drop
          drop_points <- which(drops)
          first_section_end <- drop_points[1]
          first_section_peak_idx <- which.max(group_values[1:first_section_end])
          valid_peaks_in_group <- c(group_indices[first_section_peak_idx])
          
          #cat("\nDebug: Keeping highest peak before drop:", 
          #    group_values[first_section_peak_idx], 
          #    "at time", group_times[first_section_peak_idx], "\n")
          
          # For each drop, find the highest peak between this drop and the next
          for (i in seq_along(drop_points)) {
            start_idx <- drop_points[i] + 1
            end_idx <- if (i == length(drop_points)) length(group_indices) else drop_points[i + 1]
            
            #cat("\nDebug: Processing drop point", i, "\n")
            #cat("Debug: Looking at peaks from index", start_idx, "to", end_idx, "\n")
            #cat("Debug: Values in range:", paste(group_values[start_idx:end_idx], collapse = ", "), "\n")
            
            peak_values_in_range <- group_values[start_idx:end_idx]
            max_idx_in_range <- which.max(peak_values_in_range) + start_idx - 1
            valid_peaks_in_group <- c(valid_peaks_in_group, group_indices[max_idx_in_range])
            
            #cat("Debug: Selected highest peak in range:", group_values[max_idx_in_range], 
            #    "at time", group_times[max_idx_in_range], "\n")
          }
          
          valid_peaks <- c(valid_peaks, valid_peaks_in_group)
          #cat("\nDebug: Final peaks kept in group after drops:", 
          #    paste(data$accmag[valid_peaks_in_group], collapse = ", "), "\n")
          #cat("Debug: At times:", paste(data$time[valid_peaks_in_group], collapse = ", "), "\n")
        } else {
          highest_peak <- group_indices[which.max(group_values)]
          valid_peaks <- c(valid_peaks, highest_peak)
          #cat("Debug: No drops, keeping only the highest peak:", 
          #    data$accmag[highest_peak], "at time", data$time[highest_peak], "\n")
        }
      }
    }
  } else {
    valid_peaks <- peak_indices
    #cat("\nDebug: Only one peak found, keeping it\n")
  }
  
  # Update batch_summary with peak values
  peak_values <- data$accmag[valid_peaks]
  peak_times <- data$time[valid_peaks]
  peak_pres <- data$pres[valid_peaks]
  
  #cat("\nDebug: Final Results:\n")
  #cat("Debug: Total peaks kept:", length(valid_peaks), "\n")
  #cat("Debug: Final peak values:", paste(peak_values, collapse = ", "), "\n")
  #cat("Debug: Final peak times:", paste(peak_times, collapse = ", "), "\n")
  
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
    
    if (!("guide_vane_contact" %in% colnames(batch_summary))) {
      batch_summary$guide_vane_contact <- 0
    }
    
    batch_summary <- batch_summary %>%
      mutate(
        !!accmag_col := ifelse(file == selected_sensor, peak_values[i], .data[[accmag_col]]),
        !!time_col := ifelse(file == selected_sensor, peak_times[i], .data[[time_col]]),
        !!pres_col := ifelse(file == selected_sensor, peak_pres[i], .data[[pres_col]])
      )
  }
  
  # Count peaks in each ROI
  roi_categories <- unique(data$roi)
  for (category in roi_categories) {
    if (!is.na(category)) {  # Skip NA categories
      count_col <- paste0(category, "_accmag_count")
      if (!(count_col %in% colnames(batch_summary))) {
        batch_summary[[count_col]] <- 0
      }
      
      category_peaks <- sum(data$roi[valid_peaks] == category, na.rm = TRUE)
      
      batch_summary <- batch_summary %>%
        mutate(!!count_col := ifelse(file == selected_sensor, category_peaks, .data[[count_col]]))
      
      cat("\nROI Category:", category, "\n")
      cat("Number of peaks:", category_peaks, "\n")
      if (category_peaks > 0) {
        cat("Peak times:", paste(peak_times[data$roi[valid_peaks] == category], collapse = ", "), "\n")
        cat("Peak values:", paste(peak_values[data$roi[valid_peaks] == category], collapse = ", "), "\n")
      }
    }
  }
  
  if ("post_nadir_event_accmag_count" %in% colnames(batch_summary)) {
    batch_summary <- batch_summary %>%
      mutate(guide_vane_contact = ifelse(
        file == selected_sensor & post_nadir_event_accmag_count > 0,
        1,
        guide_vane_contact
      ))
  }
  
  return(batch_summary)
}


time_normalization <- function(data, selected_sensor) {
  # normalize time series for overall passage (ROI 1) by 0 - 1 with nadir at 0.5
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
    rows_update(
      select(data_sensor, time, long_id, time_norm_new) %>%
        rename(time_norm = time_norm_new),
      by = c("time", "long_id")
    )
  
  cat("Injection to tailwater passage time normalization complete\n")
  return(data)
}

update_global_environment <- function(data, batch_summary, sample_rate, selected_sensor) {
  # Get current global data
  current_data <- if (sample_rate == 250) {
    get("data_250hz", envir = .GlobalEnv)
  } else if (sample_rate == 100) {
    get("data_100hz", envir = .GlobalEnv)
  } else if (sample_rate == "100_imp") {
    get("data_100_imp", envir = .GlobalEnv)
  } else if (sample_rate == "2000_hig") {
    get("data_2000_hig", envir = .GlobalEnv)
  }
  
  # Initialize new columns in current_data if they don't exist
  if(!"passage_point" %in% names(current_data)) {
    current_data$passage_point <- factor(NA, levels = c("post_nadir", "pre_nadir"))
  }
  if(!"roi" %in% names(current_data)) {
    current_data$roi <- factor(NA, levels = c("inj_to_pre_nadir", "pre_nadir_event", 
                                              "nadir_event", "post_nadir_event", "post_nadir_to_rec"))
  }
  if(!"time_norm" %in% names(current_data)) {
    current_data$time_norm <- 0
  }
  
  # Handle specific factor columns in update data
  data$long_id <- factor(data$long_id, levels = levels(current_data$long_id))
  if(!is.null(data$passage_point)) {
    data$passage_point <- factor(data$passage_point, levels = c("post_nadir", "pre_nadir"))
  }
  if(!is.null(data$roi)) {
    data$roi <- factor(data$roi, levels = c("inj_to_pre_nadir", "pre_nadir_event", 
                                            "nadir_event", "post_nadir_event", "post_nadir_to_rec"))
  }
  
  # Update only rows for current sensor
  current_data <- current_data %>%
    rows_update(
      data %>% filter(long_id == selected_sensor),
      by = c("time", "long_id")
    )
  
  # Assign back to global environment
  if (sample_rate == 250) {
    assign("data_250hz", current_data, envir = .GlobalEnv)
  } else if (sample_rate == 100) {
    assign("data_100hz", current_data, envir = .GlobalEnv)
  } else if (sample_rate == "100_imp") {
    assign("data_100_imp", current_data, envir = .GlobalEnv)
  } else if (sample_rate == "2000_hig") {
    assign("data_2000_hig", current_data, envir = .GlobalEnv)
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
    batch_summary <- process_nadir_value(batch_summary, selected_sensor, data)
    
    # Calculate max pressure within 1s before the nadir
    batch_summary <- calculate_max_pressure_before_nadir(data, batch_summary, selected_sensor, num_rows)
    
    # Update plot with new nadir value and 1s max pressure
    sensor_summary <- batch_summary %>% filter(file == selected_sensor)
    plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 2, num_rows)
    print(plotly_plot)
    
    # # Process 1s max nadir value. User can input new value if needed
    # batch_summary <- process_max_p_1s_nadir(batch_summary, selected_sensor)
    
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
    # plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 3)
    # print(plotly_plot)
    
    # Prompt user to enter ROI
    batch_summary <- prompt_for_injection_tailwater_times(data, batch_summary, selected_sensor, num_rows)
    
    # Assign passage points
    data <- assign_passage_points(data, batch_summary, selected_sensor, sample_rate)
    
    #Find max acceleration values within each ROI
    batch_summary <- max_acceleration_extract(data, batch_summary, selected_sensor)
    
    # Add HIG extraction for 100_imp sensors
    if (sample_rate == "100_imp") {
      batch_summary <- max_hig_extract(data_2000_hig, batch_summary, selected_sensor)
    }
    
    #determine if max acceleration in nadir roi occured before, on or after nadir, and time difference
    batch_summary <- nadir_acceleration_distance(data, batch_summary, selected_sensor)
    
    #Calculate RMS acceleration for nadir event
    batch_summary <- acceleration_rms(data, batch_summary, selected_sensor)
    
    #Find acceleration peaks with given criteria
    batch_summary <- find_acceleration_peaks(data, batch_summary, selected_sensor)
    
    #Update sensor_summary and ensure NA columns are dropped from other sensors. 
    #May be redundant, but plot fucntion relises on it at the moment
    sensor_summary <- batch_summary %>%
      filter(file == selected_sensor) %>%
      select_if(~ any(!is.na(.)))
    
    # Draw plot using ROI 1 
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
  update_global_environment(data, batch_summary, sample_rate, selected_sensor)
  
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
