#R 4.3.3
#RStudio 2024.04.2+764 "Chocolate Cosmos" Release for Windows

#27/08/24 Interim executive summary blade strike analysis

###NOTE###
#These analyese are based upon incomplete data and satisfy only the raw observed strike rates of RAPID
###NOTE###

#Required libraries####

library(tidyverse)
library(funModeling)
library(gt)
#library(plotly)
#library(cowplot)
#library(pracma)

#
BS_data=read_csv("./PROCESSED_data/UoH_JN_BSONLY_270824.csv")

#Filter data to remove unwanted columns, sensors with no video files and to RAPID only
BS_data_1 <- BS_data %>% 
  select(-c(sens_file,
            video_file,
            misc_comments,
            frame_in,
            frame_out,
            sens_wc_pos_ent,
            sens_wc_pos_blade,
            sens_ec_orien_flow,
            passage_loc,
            blade_contact_loc,
            sens_contact_loc,
            sens_velocity,
            blade_velocity,
            sens_rot,
            sens_rotation_dir,
            exit_pos,
            passage_comments
            )) %>%
  filter(video_data == 'G') %>%
  filter(!sensor %in% c('C41', 'C43', 'M28'))

#Convert logical coulmns to factors. Ensure inj_pos is assigned levels, so we know 'low' is the reference level. 
#For factors with N/Y, N is reference level by default 
#Convert n_contacts to integer as it represents count data (numerical could introduce fractional counts)

BS_data_1 <- BS_data_1 %>%
  mutate(
    treatment = as.factor(treatment),
    inj_pos = factor(inj_pos, levels = c("low", "high")),
    sensor = as.factor(sensor),
    pre_swirl = as.factor(pre_swirl),
    clear_passage = as.factor(clear_passage),
    centre_hub_contact = as.factor(centre_hub_contact),
    blade_contact = as.factor(blade_contact),
    n_contact = as.integer(n_contact)
  )

#Check structure of resulting dataframe
glimpse(BS_data_1)

#No central tendancies to report at this stage, only a single interger value (n_contacts)
#No consideration for outliers or normality need at this stage as they are 
#1) categorical data, 2) count data with no normality assumption, 3) analysed with non-parametric methods

#Summary table of descriptive statistics

calculate_proportion_ci <- function(n, x, confidence_level = 0.95) {
  # Calculate observed proportion
  p <- x / n
  
  # Calculate z-score for the confidence level
  z <- qnorm((1 + confidence_level) / 2)
  
  # Calculate margin of error
  margin_of_error <- z * sqrt((p * (1 - p)) / n)
  
  # Calculate the confidence interval
  lower_bound <- max(0, p - margin_of_error)
  upper_bound <- min(1, p + margin_of_error)
  
  # Return lower and upper bounds
  return(list(lower_bound = lower_bound, upper_bound = upper_bound))
}

summary_stats <- function(data) {
  data %>%
    summarise(
      total_rapid = n(),
      total_no_collision = sum(clear_passage == "Y"),
      total_collision = sum(clear_passage == "N"),
      total_blade_strike = sum(blade_contact == "Y"),
      single_contact = sum(blade_contact == "Y" & n_contact == 1),
      multiple_contact = sum(blade_contact == "Y" & n_contact > 1),
      total_contacts = sum(n_contact),
      impeller_collision_prob = total_collision / total_rapid,
      blade_strike_prob = total_blade_strike / total_rapid
    ) %>%
    mutate(
      impeller_collision_ci = map2(total_rapid, total_collision, ~calculate_proportion_ci(..1, ..2)),
      blade_strike_ci = map2(total_rapid, total_blade_strike, ~calculate_proportion_ci(..1, ..2))
    )
}

summary_long <- BS_data_1 %>%
  group_by(treatment) %>%
  summary_stats() %>%
  ungroup() %>%
  mutate(
    impeller_collision_lower = map_dbl(impeller_collision_ci, "lower_bound"),
    impeller_collision_upper = map_dbl(impeller_collision_ci, "upper_bound"),
    blade_strike_lower = map_dbl(blade_strike_ci, "lower_bound"),
    blade_strike_upper = map_dbl(blade_strike_ci, "upper_bound")
  ) %>%
  select(-impeller_collision_ci, -blade_strike_ci)




formatted_st <- summary_long %>%
  mutate(
    `Pump operation` = case_when(
      treatment == "500_1" ~ "500 RPM (100% BEP)",
      treatment == "400_1" ~ "400 RPM (100% BEP)",
      treatment == "400_07" ~ "400 RPM (70% BEP)",
      treatment == "400_05" ~ "400 RPM (50% BEP)"
    ),
    `Flow rate m³/s` = case_when(
      treatment == "500_1" ~ 1.32,
      treatment == "400_1" ~ 1.05,
      treatment == "400_07" ~ 0.73,
      treatment == "400_05" ~ 0.52
    ),
    impeller_collision_prob = sprintf("%.1f%% [%.1f, %.1f]", 
                                      impeller_collision_prob * 100, 
                                      impeller_collision_lower * 100, 
                                      impeller_collision_upper * 100),
    blade_strike_prob = sprintf("%.1f%% [%.1f, %.1f]", 
                                blade_strike_prob * 100, 
                                blade_strike_lower * 100, 
                                blade_strike_upper * 100)
  ) %>%
  arrange(desc(treatment)) %>%
  select(`Pump operation`, `Flow rate m³/s`, total_rapid, total_no_collision, total_collision, impeller_collision_prob, 
         total_blade_strike, single_contact, multiple_contact, total_contacts, blade_strike_prob) %>%
  gt() %>%
  cols_label(
    `Flow rate m³/s` = "Flow rate (m³/s)",
    total_rapid = "Total sensor observation",
    total_no_collision = "Total impeller passage",
    total_collision = "Total impeller collision",
    impeller_collision_prob = "Impeller collision probability [95% CI]",
    total_blade_strike = "Total sensor blade strike",
    single_contact = md("*n single contact*"),
    multiple_contact = md("*n multiple contact*"),
    total_contacts = md("*total contacts*"),
    blade_strike_prob = "Blade strike probability* [95% CI]"
  ) %>%
  fmt_number(
    columns = c(total_rapid, total_no_collision, total_collision, total_blade_strike, 
                single_contact, multiple_contact, total_contacts),
    decimals = 0
  ) %>%
  tab_style(
    style = cell_fill(color = "#E8F5E9"), 
    locations = cells_body(columns = total_no_collision)
  ) %>%
  tab_style(
    style = cell_fill(color = "#FFE0B2"),  
    locations = cells_body(columns = c(total_collision, impeller_collision_prob))
  ) %>%
  tab_style(
    style = cell_fill(color = "#FFCC80"),  
    locations = cells_body(
      columns = c(total_blade_strike, single_contact, multiple_contact, total_contacts, blade_strike_prob)
    )
  ) %>%
  text_transform(
    locations = cells_body(columns = c(impeller_collision_prob, blade_strike_prob)),
    fn = function(x) paste0("<u>", x, "</u>")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_options(
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.top.color = "black",
    table.border.bottom.style = "solid",
    table.border.bottom.color = "black",
    table.border.left.style = "none",
    table.border.right.style = "none",
    column_labels.border.top.style = "solid",
    column_labels.border.top.color = "black",
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.style = "solid",
    table_body.border.bottom.color = "black",
    table_body.hlines.style = "none"
  )

formatted_st

gtsave(formatted_st, filename = "./R_output_files/blade_strike_table.html")








max_contacts <- 120

blade_strike_plot <- summary_long %>%
  mutate(treatment = factor(treatment, levels = c("500_1", "400_1", "400_07", "400_05"),
                            labels = c("500 RPM\n(100% BEP)", "400 RPM\n(100% BEP)", "400 RPM\n(70% BEP)", "400 RPM\n(50% BEP)"))) %>%
  ggplot(aes(x = treatment)) +
  geom_point(aes(y = blade_strike_prob), size = 3) +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = 0.1) + 
  geom_text(aes(y = blade_strike_prob, label = sprintf("%.1f%%", blade_strike_prob)), 
            vjust = 2, 
            hjust = -0.2,
            fontface = "bold") +
  geom_line(aes(y = blade_strike_prob, group = 1), linetype = "solid") +
  geom_text(aes(y = 0, label = sprintf("n = %d", total_rapid)), 
            vjust = 0, 
            fontface = "italic", 
            size = 8 / .pt)+
  geom_line(aes(y = total_contacts / max_contacts  * 100, group = 1), linetype = "dashed", , position = position_nudge(x=0.2)) +
  geom_point(aes(y = total_contacts / max_contacts  * 100), size = 3, shape = 2, position = position_nudge(x=0.2)) +
  geom_text(aes(y = total_contacts / max_contacts  * 100, label = total_contacts),
            hjust = -1.5,
            vjust = -1,
            size = 8 / .pt) +
  annotate("text", x = 0.5, y = 100, label = "CI = 95%", hjust = 0, vjust = 1, fontface = "bold") +
    scale_y_continuous(
    name = "Blade strike probability",
    labels = scales::percent_format(scale = 1), 
    limits = c(0, 100),
    sec.axis = sec_axis(~ . * max_contacts / 100, 
                        name = "Total blade contacts",
                        breaks = seq(0, max_contacts, by = 20))
  ) +
  labs(x = "Pump operation") +
  coord_cartesian(clip = "off") +
  theme_JN() +
  theme(axis.text.x = element_text(lineheight = 0.8))


ggsave(filename="./R_output_files/blade_strike_plot.svg", plot=blade_strike_plot,device = "svg",units="cm",width=14,height=10)

chisq_result <- chisq.test(table(BS_data_1$treatment, BS_data_1$blade_contact))
print(chisq_result)

blade_strike_table <- table(BS_data_1$treatment, BS_data_1$blade_contact)
pairwise_prop <- pairwise.prop.test(blade_strike_table[,"Y"], 
                                    rowSums(blade_strike_table), 
                                    p.adjust.method = "bonferroni")
print(pairwise_prop)


pairwise_prop_bh <- pairwise.prop.test(blade_strike_table[,"Y"], 
                                       rowSums(blade_strike_table), 
                                       p.adjust.method = "BH")
print(pairwise_prop_bh)




rpm_400_data <- BS_data_1 %>%
  filter(treatment %in% c("400_1", "400_07", "400_05")) %>%
  mutate(efficiency = case_when(
    treatment == "400_1" ~ 100,
    treatment == "400_07" ~ 70,
    treatment == "400_05" ~ 50
  ),
  blade_strike = as.numeric(blade_contact == "Y"))

cor_test_individual <- cor.test(rpm_400_data$efficiency, rpm_400_data$blade_strike, 
                                method = "spearman", exact = FALSE)

print(cor_test_individual)

cor_test_contacts <- cor.test(rpm_400_data$efficiency, rpm_400_data$n_contact, 
                              method = "spearman", exact = FALSE)

print(cor_test_contacts)


# Filter the data for the two treatments we want to compare
contact_data <- BS_data_1 %>%
  filter(treatment %in% c("500_1", "400_1"))

# Perform the Mann-Whitney U test
wilcox_test <- wilcox.test(n_contact ~ treatment, data = contact_data)

print(wilcox_test)

