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

BDSAnalysisTool(batch_summary, data_250hz, data_100hz, 100)

#After analysis, add treatment variables
add_treatment(data_250hz, "data_250hz")

#When batch complete, crop data to passage points, save and export
#If continuing with analysis, remove NA passage_points to reduce memory consumption from data files
save_data(batch_summary, data_100hz, data_250hz)

batch_summary <- batch_summary %>%
  filter(pres_processed != 'N'| acc_processed != 'N')




primary_to_secondary <- function(x) {
  (x - summary_stats$min_pres) / (summary_stats$max_pres - summary_stats$min_pres) * (summary_stats$max_accmag - summary_stats$min_accmag) + summary_stats$min_accmag
}

summary_stats <- data_100hz %>%
  filter(!is.na(passage_point) & long_id == 'C330522135009') %>%
  summarize(
    min_pres = floor(min(pres) / 50) * 50,
    max_pres = ceiling(max(pres) / 50) * 50,
    min_accmag = floor(min(accmag, na.rm = TRUE) / 10) * 10,
    max_accmag = ceiling(max(accmag, na.rm = TRUE) / 10) * 10
  )

primary_breaks <- seq(summary_stats$min_pres, summary_stats$max_pres, by = 50)
secondary_breaks <- seq(summary_stats$min_accmag, summary_stats$max_accmag, by = 10)

# Calculate other necessary statistics and merge
filtered_data <- data_100hz %>%
  filter(!is.na(passage_point) & long_id == 'C330522135009') %>%
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
      name = expression("Acceleration Magnitude (accmag) (m/s"^2*")"),
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
      name = expression("Acceleration Magnitude (accmag) (m/s"^2*")"),
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


ggsave(filename="pres_acc_plot.png", plot=pres_acc_plot, device = "png",units="cm",width=18,height=20)





data_100hz=read_csv("./R_output_files/230524_test/test_data_100hz.csv")
batch_summary=read_csv("./R_output_files/230524_test/test_batch_summary.csv")
  

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