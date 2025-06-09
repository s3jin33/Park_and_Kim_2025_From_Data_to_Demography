library(dplyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(readr)
library(readxl)
library(stringr)
library(tidyr)
library(writexl)
library(ggsci)
library(patchwork)
install.packages(ggtext)


getwd()
setwd("/Volumes/SAMSUNG/23-2학회/BHDC/Revision")

# 
# #Load Data
# Data<-read.csv("After_Combine_250602.csv")
# 
# # Calibration
# x <- calibrate(x = Data$BP, errors = Data$error, calCurves = 'intcal20')
# 
# # Summary
# summaryx <- summary(x)
# 
# # MedianBP 
# medians <- as.numeric(summaryx$MedianBP)
# 
# # Data_clean에 붙이기 (행 수 일치하는지 반드시 확인)
# Data$median <- medians
# 
# write_xlsx(Data, "After_Combine_250604.xlsx")

Data<-read_csv("After_Combine_250604.csv")

Data_clean <- Data %>% filter(is.na(Mismatch_Flag))


Site_Minmax <- Data_clean %>% 
  group_by(Site_Code, Sub_Region, Total_Houses, Period_New) %>% 
  summarise(
    median_min = min(median, na.rm = TRUE),
    median_max = max(median, na.rm = TRUE),
    .groups = "drop"
  )

# Counting Total Houses per sub region
Regional_Total_Houses <- Site_Minmax %>%
  group_by(Sub_Region) %>%
  summarise(
    total_houses = sum(Total_Houses, na.rm = TRUE),
    site_count = n_distinct(Site_Code),
    .groups = "drop"
  )


Regional_Total_Houses_Period <- Site_Minmax %>%
  group_by(Sub_Region,Period_New) %>%
  summarise(
    total_houses = sum(Total_Houses, na.rm = TRUE),
    site_count = n_distinct(Site_Code),
    .groups = "drop"
  )


# Create a sequence of time intervals
intervals <- seq(3500, 1100, by = -200)

# Function to calculate phases and phases string for each site
calculate_phases <- function(site_data) {
  presence_phases <- 0
  phases_string <- character(0)
  
  for (i in (length(intervals) - 1):1) {
    interval_start <- intervals[i]
    interval_end <- intervals[i + 1]
    
    min_val <- as.numeric(site_data[["median_min"]])
    max_val <- as.numeric(site_data[["median_max"]])
    
    if (is.na(min_val) || is.na(max_val)) next
    
    if (max_val >= interval_end && min_val <= interval_start) {
      presence_phases <- presence_phases + 1
      phases_string <- c(phases_string, sprintf("%d-%d", interval_start, interval_end))
    }
  }
  
  phases_string <- rev(phases_string)
  return(c(presence_phases, paste(phases_string, collapse = ",")))
}

results <- t(apply(Site_Minmax, 1, calculate_phases))

# Add results to the data frame
Site_Minmax$presence_phases <- as.integer(results[, 1])
Site_Minmax$phases_string <- results[, 2]

# Calculate values based on number of presence phases
Site_Minmax$phase_weight <- round(1 / Site_Minmax$presence_phases, 2)
Site_Minmax$adjusted_house_count <- Site_Minmax$phase_weight * Site_Minmax$Total_Houses

# Remove sites with zero presence phases
Site_Minmax <- Site_Minmax %>% filter(presence_phases != 0)

# Expand rows so that each phase has a separate row
Site_Minmax <- separate_rows(Site_Minmax, phases_string, sep = ",")

# List of all unique phases present
unique(Site_Minmax$phases_string)

# Define phase orders and labels
phase_order <- c("3500-3300", "3300-3100", "3100-2900", "2900-2700",
                 "2700-2500", "2500-2300", "2300-2100", "2100-1900",
                 "1900-1700", "1700-1500", "1500-1300", "1300-1100")

phase_labels <- c("35-33", "33-31", "31-29", "29-27",
                  "27-25", "25-23", "23-21", "21-19",
                  "19-17", "17-15", "15-13", "13-11")

Site_Minmax$Sub_Region <- factor(Site_Minmax$Sub_Region, levels = c('Northern Gyeonggi',
                                                                    'Southern Gyeonggi',
                                                                    'North Han River Basin',
                                                                    'South Han River Basin',
                                                                    'West Coast'))



Site_Minmax_total <- Site_Minmax %>%
  group_by(phases_string) %>%
  summarise(
    phase_weight = sum(phase_weight, na.rm = TRUE),
    adjusted_house_count = sum(adjusted_house_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Sub_Region = "Han River Basin(Total)")

# Total plot1
plot1_total <- Site_Minmax_total %>%
  ggplot(aes(x = phases_string, y = phase_weight)) +
  geom_bar(stat = "identity", fill = "#80796BFF") + 
  facet_wrap(~ Sub_Region) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 40, face = "bold"),
    text = element_text(size = 40, face = "bold"),
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),   
    strip.text = element_text(size = 50, face = "bold")
  ) +
  scale_x_discrete(limits = phase_order, labels = phase_labels)

# Total plot2
plot2_total <- Site_Minmax_total %>%
  ggplot(aes(x = phases_string, y = adjusted_house_count)) +
  geom_bar(stat = "identity", fill = "#80796BFF") +
  facet_wrap(~ Sub_Region) +
  scale_y_continuous(position = "right") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 40, face = "bold"),
    text = element_text(size = 40, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),     # y축 제목 제거
    strip.text = element_text(size = 50, face = "bold")
  ) +
  scale_x_discrete(limits = phase_order, labels = phase_labels) 

# Sub_Region plot1 


plot1_regions <- Site_Minmax %>%
  ggplot(aes(x = phases_string, y = phase_weight, fill = Sub_Region)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sub_Region, nrow = 5) +
  scale_fill_jama() +
  coord_cartesian(ylim = c(0, 30)) +  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 40, face = "bold"),
    axis.text.y = element_text(size = 40, face = "bold"),
    text = element_text(size = 40, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 40), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 5), face = "bold"),
    strip.text = element_text(size = 50, face = "bold")
  ) +
  scale_x_discrete(limits = phase_order, labels = phase_labels) +
  labs(x = "kyr cal BP", y = "Weighted Site Counts")


# Sub_Region plot2 
plot2_regions <- Site_Minmax %>%
  ggplot(aes(x = phases_string, y = adjusted_house_count, fill = Sub_Region)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sub_Region, nrow = 5) +
  scale_fill_jama() +
  coord_cartesian(ylim = c(0, 600)) +  
  scale_y_continuous(position = "right", limits = c(0, 600)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 40, face = "bold"),
    axis.text.y = element_text(size = 40, face = "bold"),
    text = element_text(size = 40, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 40), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 5), face = "bold"),
    strip.text = element_text(size = 50, face = "bold")
  ) +
  scale_x_discrete(limits = phase_order, labels = phase_labels) +
  labs(x = "kyr cal BP", y = "Weighted House Counts")


plot1_total <- plot1_total + 
  ggtitle("(a) Weighted Values of Site Counts") +
  theme(plot.title = element_text(size = 43, face = "bold", hjust = 0.8))

plot2_total <- plot2_total + 
  ggtitle("(b) Weighted Values of Total Houses") +
  theme(plot.title = element_text(size = 43, face = "bold", hjust = 0.8))

final_plot <- (plot1_total | plot2_total) / (plot1_regions | plot2_regions) +
  plot_layout(heights = c(1, 5))

ggsave("combined_graphs_differential_y_JAMA_250604.png", 
       plot = final_plot,
       width = 24, 
       height = 32, 
       units = "in", 
       dpi = 500)
