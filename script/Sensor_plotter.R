#Simple plotting tool for batch processing sensors. Lots of bugs.

#Development reasons only, uncomment if needed
#batch_summary <-batch_summary %>%
#  mutate(plotted = "N")

plot_sensors <- function(data_100hz, data_250hz, data_100_imp, sample_rate) {
  batch_summary <- get("batch_summary", envir = .GlobalEnv)
  
  # Initialize 'plotted' if it doesn't exist
  if (!"plotted" %in% colnames(batch_summary)) {
    batch_summary$plotted <- "N"
  }
  
  repeat {
    # Identify sensors based on sample rate
    sensors <- batch_summary %>%
      filter(plotted == "N" & (pres_processed == "Y" | acc_processed == "Y") & class == sample_rate) %>%
      pull(file) %>%
      unique()
    
    if (length(sensors) == 0) {
      cat("All sensors plotted for sample rate", sample_rate, "\n")
      cat("Plotting processes ended\n")
      break
    }
    
    # User selection with numbered list
    cat("Select a sensor dataset to plot pressure and acceleration profiles:\n")
    for (i in seq_along(sensors)) {
      cat(i, ": ", sensors[i], "\n", sep = "")
    }
    cat(length(sensors) + 1, ": All\n", sep = "")
    
    sensor_choice_index <- as.integer(readline(prompt = "Enter the number corresponding to the sensor: "))
    if (is.na(sensor_choice_index) || sensor_choice_index < 1 || sensor_choice_index > length(sensors) + 1) {
      cat("Invalid sensor selection. Please try again.\n")
      next
    }
    
    if (sensor_choice_index == length(sensors) + 1) {
      sensors_to_plot <- sensors
    } else {
      sensors_to_plot <- sensors[sensor_choice_index]
    }
    
    for (sensor in sensors_to_plot) {
      cat("Plotting pressure and acceleration for sensor", sensor, "\n")
      
      # Determine the dataset based on sample rate
      if (sample_rate == "250hz") {
        data <- data_250hz
      } else if (sample_rate == "100hz") {
        data <- data_100hz
      } else if (sample_rate == "100_imp") {
        data <- data_100_imp
      } else {
        cat("Invalid sample rate for sensor", sensor, "\n")
        next
      }
      
      # Filter data for the selected sensor
      filtered_data <- data %>%
        filter(!is.na(passage_point) & long_id == sensor)
      
      if (nrow(filtered_data) == 0) {
        cat("No data available for sensor", sensor, "\n")
        next
      }
      # Filter data for the selected sensor
      summary_stats <- filtered_data %>%
        summarize(
          min_pres = floor(min(pres) / 50) * 50,
          max_pres = ceiling(max(pres) / 50) * 50,
          min_accmag = floor(min(accmag, na.rm = TRUE) / 10) * 10,
          max_accmag = ceiling(max(accmag, na.rm = TRUE) / 10) * 10
        )
      
      primary_breaks <- seq(summary_stats$min_pres, summary_stats$max_pres, by = 50)
      secondary_breaks <- seq(summary_stats$min_accmag, summary_stats$max_accmag, by = 10)
      
      filtered_data <- filtered_data %>%
        mutate(time_seconds = (row_number() - 1) * 0.01) %>%
        filter(time_seconds <= floor(max(time_seconds))) %>%
        left_join(batch_summary %>% select(file, nadir), by = c("long_id" = "file")) %>%
        filter(!is.na(accmag)) %>%
        mutate(
          min_time = min(time_seconds),
          max_time = max(time_seconds),
          y_center = mean(pres, na.rm = TRUE),
          min_time_seconds = floor(min(time_seconds)),
          max_time_seconds = ceiling(max(time_seconds)),
          nadir_time = time_seconds[which.min(abs(pres - nadir))],
          rounded_nadir_time = round(nadir_time),
          start_time = round(nadir_time) - 2,
          end_time = round(nadir_time) + 2,
          pres_range = summary_stats$max_pres - summary_stats$min_pres,
          accmag_range = summary_stats$max_accmag - summary_stats$min_accmag,
          normalized_accmag = (accmag - summary_stats$min_accmag) / accmag_range * pres_range + summary_stats$min_pres
        )
      
      min_time_seconds <- filtered_data %>% pull(min_time_seconds) %>% unique()
      max_time_seconds <- filtered_data %>% pull(max_time_seconds) %>% unique()
      nadir_time <- filtered_data$nadir_time[1]
      nadir <- filtered_data$nadir[1]
      
      primary_to_secondary <- function(x) {
        (x - summary_stats$min_pres) / (summary_stats$max_pres - summary_stats$min_pres) * (summary_stats$max_accmag - summary_stats$min_accmag) + summary_stats$min_accmag
      }
      
      plot_1 <- ggplot(data = filtered_data, aes(x = time_seconds, y = pres)) +
        geom_line() +
        geom_vline(xintercept = filtered_data$min_time[1], linetype = "dashed") +
        geom_vline(xintercept = filtered_data$max_time[1], linetype = "dashed") +
        geom_point(aes(x = nadir_time, y = nadir), shape = 21, fill = 'orange', size = 2) +
        annotate("text", x = filtered_data$min_time[1], y = filtered_data$y_center[1], label = "Injection", angle = 90, vjust = 1.5) +
        annotate("text", x = filtered_data$max_time[1], y = filtered_data$y_center[1], label = "Recovery", angle = 90, vjust = -1.5) + 
        annotate("text", x = filtered_data$nadir_time[1], y = filtered_data$nadir[1], label = paste("Pressure nadir:", round(filtered_data$nadir[1], 2)), hjust = -0.1) +
        scale_y_continuous(
          breaks = primary_breaks,
          limits = c(summary_stats$min_pres, summary_stats$max_pres)
        ) +
        scale_x_continuous(
          breaks = seq(min_time_seconds, max_time_seconds, by = 1),
          limits = c(min_time_seconds, max_time_seconds)
        ) +
        labs(x = "Passage time (seconds)", y = "Pressure (mbar)") +
        coord_cartesian(clip = "off") +
        theme_JN() +
        theme(
          plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"),
          panel.grid.minor.y = element_line(color = 'lightgrey'),
          panel.grid.major.y = element_line(color = 'grey')
        )
      
      plot_2 <- ggplot(data = filtered_data %>% filter(time_seconds >= start_time & time_seconds <= end_time), aes(x = time_seconds, y = pres)) +
        geom_line() +
        geom_point(aes(x = nadir_time, y = nadir), shape = 21, fill = 'orange', size = 2) +
        annotate("text", x = nadir_time, y = nadir, label = paste("Pressure nadir:", round(nadir, 2)), hjust = -0.1) +
        scale_y_continuous(
          breaks = primary_breaks,
          limits = c(summary_stats$min_pres, summary_stats$max_pres)
        ) +
        scale_x_continuous(
          breaks = seq(filtered_data$start_time[1], filtered_data$end_time[1], by = 0.5),
          limits = c(filtered_data$start_time[1], filtered_data$end_time[1]), expand = c(0, 0)
        ) +
        labs(x = "Passage time (seconds)", y = "Pressure (mbar)") +
        coord_cartesian(clip = "off") +
        theme_JN() +
        theme(
          plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"),
          panel.grid.minor.y = element_line(color = 'lightgrey'),
          panel.grid.major.y = element_line(color = 'grey')
        )
      
      plot_3 <- ggplot(data = filtered_data, aes(x = time_seconds)) +
        geom_line(aes(y = pres), color = "black") +
        geom_line(aes(y = normalized_accmag), color = "blue") +
        geom_vline(xintercept = filtered_data$min_time[1], linetype = "dashed") +
        geom_vline(xintercept = filtered_data$max_time[1], linetype = "dashed") +
        geom_point(aes(x = nadir_time, y = nadir), shape = 21, fill = 'orange', size = 2) +
        annotate("text", x = filtered_data$min_time[1], y = filtered_data$y_center[1], label = "Injection", angle = 90, vjust = 1.5) +
        annotate("text", x = filtered_data$max_time[1], y = filtered_data$y_center[1], label = "Recovery", angle = 90, vjust = -1.5) +
        annotate("text", x = filtered_data$nadir_time[1], y = filtered_data$nadir[1], label = paste("Pressure nadir:", round(filtered_data$nadir[1], 2)), hjust = -0.1) +
        scale_y_continuous(
          name = "Pressure (mbar)",
          breaks = primary_breaks,
          limits = c(summary_stats$min_pres, summary_stats$max_pres),
          sec.axis = sec_axis(
            trans = ~ primary_to_secondary(.),
            name = expression("Acceleration Magnitude (m/s"^2*")"),
            breaks = secondary_breaks
          )
        ) +
        scale_x_continuous(
          breaks = seq(min_time_seconds, max_time_seconds, by = 1),
          limits = c(min_time_seconds, max_time_seconds)
        ) +
        labs(x = "Passage time (seconds)") +
        coord_cartesian(clip = "off") +
        theme_JN() +
        theme(
          plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"),
          panel.grid.minor.y = element_line(color = 'lightgrey'),
          panel.grid.major.y = element_line(color = 'grey')
        )
      
      plot_4 <- ggplot(data = filtered_data %>% filter(time_seconds >= start_time & time_seconds <= end_time), aes(x = time_seconds)) +
        geom_line(aes(y = pres), color = "black") +
        geom_line(aes(y = normalized_accmag), color = "blue") +
        geom_point(aes(x = nadir_time, y = nadir), shape = 21, fill = 'orange', size = 2) +
        annotate("text", x = filtered_data$nadir_time[1], y = filtered_data$nadir[1], label = paste("Pressure nadir:", round(filtered_data$nadir[1], 2)), hjust = -0.1) +
        scale_y_continuous(
          name = "Pressure (mbar)",
          breaks = primary_breaks,
          limits = c(summary_stats$min_pres, summary_stats$max_pres),
          sec.axis = sec_axis(
            trans = ~ primary_to_secondary(.),
            name = expression("Acceleration Magnitude (m/s"^2*")"),
            breaks = secondary_breaks
          )
        ) +
        scale_x_continuous(
          breaks = seq(filtered_data$start_time[1], filtered_data$end_time[1], by = 0.5),
          limits = c(filtered_data$start_time[1], filtered_data$end_time[1]), expand = c(0, 0)
        ) +
        labs(x = "Passage time (seconds)") +
        coord_cartesian(clip = "off") +
        theme_JN() +
        theme(
          plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"),
          panel.grid.minor.y = element_line(color = 'lightgrey'),
          panel.grid.major.y = element_line(color = 'grey')
        )
      
      pressure_plot <- plot_grid(plot_1, plot_2,
                                 ncol = 1, nrow = 2, align = "h")
      
      pres_acc_plot <- plot_grid(plot_3, plot_4,
                                 ncol = 1, nrow = 2, align = "h")
      
      ggsave(filename=paste0("./R_output_files/", sensor, "_pressure_plot.png"), plot=pressure_plot, device = "png", units="cm", width=18, height=20)
      ggsave(filename=paste0("./R_output_files/", sensor, "_pres_acc_plot.png"), plot=pres_acc_plot, device = "png", units="cm", width=18, height=20)
      
     # Update 'plotted' status
      batch_summary$plotted[batch_summary$file == sensor] <- "Y"
      
      cat("Plotting complete for sensor", sensor, "\n")
    }
    
    # Save the updated batch_summary to the global environment
    assign("batch_summary", batch_summary, envir = .GlobalEnv)
    
    # Check if all sensors are plotted
    if (all(batch_summary$plotted == "Y")) {
      cat("All sensors plotted for sample rate", sample_rate, "\n")
      cat("Plotting processes ended\n")
      break
    } else {
      continue_choice <- readline("Do you want to plot another sensor? (Y/N): ")
      if (toupper(continue_choice) != "Y") {
        cat("Plotting processes ended\n")
        break
      }
    }
  }
}

plot_sensors(data_100hz, data_250hz, data_100_imp, "100hz")

#Input options
# "250hz"
# "100hz"
# "100_imp"