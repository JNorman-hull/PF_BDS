#development script

#User choice ROI selection

prompt_for_injection_tailwater_times <- function(data, batch_summary, selected_sensor, num_rows) {
  #ROI selection tool
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
      row_start <- max(1, nadir_index - num_rows * 1.5)
      row_end <- min(nrow(sensor_data), nadir_index + num_rows * 1.5)
      t_start_val <- sensor_data$time[row_start]
      t_end_val <- sensor_data$time[row_end]
      cat("Default start time for PF (nadir - 1.5s):", t_start_val, "\n")
      cat("Default end time for PF (nadir + 1.5s):", t_end_val, "\n")
      prompt <- "Use these default times for the main passage ROI? (Y/N): "
      
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
      
    } else if (roi == 2) {
      cat(" (Nadir event ROI)\n")
      row_start <- max(1, nadir_index - 10)
      row_end <- min(nrow(sensor_data), nadir_index + 10)
      t_start_val <- sensor_data$time[row_start]
      t_end_val <- sensor_data$time[row_end]
      cat("Default nadir event start time for PF (nadir - 0.1s):", t_start_val, "\n")
      cat("Default nadir event end time for PF (nadir + 0.1s):", t_end_val, "\n")
      prompt <- "Use these default times for nadir event? (Y/N): "
      
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
      
    } else if (roi == 3) {
      cat(" (Pre-nadir ROI)\n")
      t_start_val <- sensor_data$time[max(1, which(sensor_data$time == batch_summary$t_start_roi2[batch_summary$file == selected_sensor]) - num_rows / 3.2)]
      t_end_val <- sensor_data$time[which(sensor_data$time == batch_summary$t_start_roi2[batch_summary$file == selected_sensor])]
      cat("Default start time for PF (nadir ROI start - 0.3s):", t_start_val, "\n")
      cat("End time:", t_end_val, "\n")
      prompt <- "Use default start time for pre-nadir ROI? (Y/N): "
      
      use_default <- readline(prompt = prompt)
      
      if (toupper(use_default) != "Y") {
        cat("Enter custom start time for pre-nadir ROI\n")
        t_start_val <- as.numeric(readline(prompt = "Enter start time value: "))
        
        while (is.na(t_start_val) || t_start_val >= t_end_val) {
          cat("Invalid input. Please ensure start time is less than end time and is numeric.\n")
          t_start_val <- as.numeric(readline(prompt = "Enter start time value: "))
        }
      }
      
    } else if (roi == 4) {
      cat(" (Post-nadir ROI)\n")
      t_start_val <- sensor_data$time[which(sensor_data$time == batch_summary$t_end_roi2[batch_summary$file == selected_sensor])]
      t_end_val <- sensor_data$time[min(nrow(sensor_data), which(sensor_data$time == batch_summary$t_end_roi2[batch_summary$file == selected_sensor]) + num_rows /3.2)]
      cat("Start time:", t_start_val, "\n")
      cat("Default end time for PF (nadir ROI end + 0.3s):", t_end_val, "\n")
      prompt <- "Use default end time for post-nadir ROI? (Y/N): "
      
      use_default <- readline(prompt = prompt)
      
      if (toupper(use_default) != "Y") {
        cat("Enter custom end time for post-nadir ROI\n")
        t_end_val <- as.numeric(readline(prompt = "Enter end time value: "))
        
        while (is.na(t_end_val) || t_end_val <= t_start_val) {
          cat("Invalid input. Please ensure end time is greater than start time and is numeric.\n")
          t_end_val <- as.numeric(readline(prompt = "Enter end time value: "))
        }
      }
      
    } else if (roi == 5) {
      t_start_val <- batch_summary$t_start_roi1[batch_summary$file == selected_sensor]
      t_end_val <- batch_summary$t_start_roi3[batch_summary$file == selected_sensor]
      cat("ROI 5 automatically labelled (inj_to_pre_nadir)\n")
      cat("Start time:", t_start_val, "\n")
      cat("End time:", t_end_val, "\n")
    } else if (roi == 6) {
      t_start_val <- batch_summary$t_end_roi4[batch_summary$file == selected_sensor]
      t_end_val <- batch_summary$t_end_roi1[batch_summary$file == selected_sensor]
      cat("ROI 6 automatically labelled (post_nadir_to_rec)\n")
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




#####
#removed user choice for acceleration from BDSANlysisTool
#Find acceleration peaks with given criteria and show plot in a loop until user is satisfied
repeat {
  # Clear previous peak data from batch_summary
  if (exists("peak")) {  # Only clear if not first run
    # Clear accmag_n columns
    batch_summary <- batch_summary %>%
      select(-matches("^accmag_\\d+$"), 
             -matches("^accmag_\\d+_t$"),
             -matches("^accmag_\\d+_p$"),
             -matches("_accmag_count$"))
  }
  
  # Apply peak detection
  if (!exists("peak")) {
    batch_summary <- find_acceleration_peaks(data, batch_summary, selected_sensor)
  } else {
    batch_summary <- find_acceleration_peaks(data, batch_summary, selected_sensor,
                                             peak = peak, peak_gap = peak_gap, 
                                             drop = drop, drop_gap = drop_gap,
                                             group_window_multiplier = group_window_multiplier)
  }
  
  #Update sensor_summary and ensure NA columns are dropped from other sensors
  sensor_summary <- batch_summary %>%
    filter(file == selected_sensor) %>%
    select_if(~ any(!is.na(.)))
  
  # Draw plot using ROI 1
  plotly_plot <- create_combined_plot(data, sensor_summary, selected_sensor, stage = 4)
  print(plotly_plot)
  
  # Check if peaks are satisfactory
  peaks_ok <- readline(prompt = "Are the acceleration peaks satisfactory? (Y/N): ")
  if (toupper(peaks_ok) == "Y") {
    break
  } else {
    # Get new parameters with validation
    repeat {
      peak <- suppressWarnings(as.numeric(readline(prompt = "Enter new peak threshold (default 49.03): ")))
      if (!is.na(peak)) break
      cat("Please enter a valid number\n")
    }
    repeat {
      peak_gap <- suppressWarnings(as.numeric(readline(prompt = "Enter new peak gap (default 5): ")))
      if (!is.na(peak_gap)) break
      cat("Please enter a valid number\n")
    }
    repeat {
      drop <- suppressWarnings(as.numeric(readline(prompt = "Enter new drop threshold (default 5): ")))
      if (!is.na(drop)) break
      cat("Please enter a valid number\n")
    }
    repeat {
      drop_gap <- suppressWarnings(as.numeric(readline(prompt = "Enter new drop gap (default 1): ")))
      if (!is.na(drop_gap)) break
      cat("Please enter a valid number\n")
    }
    repeat {
      group_window_multiplier <- suppressWarnings(as.numeric(readline(prompt = "Enter new group window multiplier (default 3): ")))
      if (!is.na(group_window_multiplier)) break
      cat("Please enter a valid number\n")
    }
  }
}

