install.packages("tidyverse")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("knitr") # For kable
install.packages("kableExtra") # For kable styling

library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(tidyr) # For pivot_longer

# Load the human and mouse data 
base_path <- "/Users/addisonyam/Documents/KLF4/"
human_data <- read.csv(file.path(base_path, "Human_motif.csv"))
mouse_data <- read.csv(file.path(base_path, "Mouse_motif.csv"))

human_cell_lines <- c("BJ", "H9_ESC", "iPSC", "LIS49_hESC", "MCF7", "U87", "HN_SCC")
mouse_cell_lines <- c("bone_marrow", "Mbd3f", "MEF", "SCC", "V6_5_ESC", "CD19", "ESC", "pre_iPSC", "Prostate_Stem")

# Function to process data for a given species
process_motif_data <- function(data, cell_lines, species_name) {
  results <- data.frame(
    Cell_Line = character(),
    Total_KLF4_Bindings = numeric(),
    KLF4_Motif_Count = numeric(),
    KLF4_Motif_Percentage = numeric(),
    Average_Klf4_Motif_Score = numeric(),
    SD_Klf4_Motif_Score = numeric(),
    Total_Number_of_Klf4_Motifs = numeric(),
    Average_Total_Motif_Score = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col_name in cell_lines) {
    # Filter rows where the current cell line column is 1 (indicating KLF4 binding)
    cell_line_bound_data <- data %>% filter(!!sym(col_name) == 1)
    
    total_bindings <- nrow(cell_line_bound_data)
    print(cell_line_bound_data)
    klf4_motif_count <- sum(cell_line_bound_data$`Klf4.motif` == 1, na.rm = TRUE)
    
    klf4_motif_percentage <- ifelse(total_bindings > 0, (klf4_motif_count / total_bindings) * 100, 0)
    
    avg_klf4_motif_score <- ifelse(all(is.na(cell_line_bound_data$`Klf4.motif.score`)), NA, mean(cell_line_bound_data$`Klf4.motif.score`, na.rm = TRUE))
    sd_klf4_motif_score <- ifelse(all(is.na(cell_line_bound_data$`Klf4.motif.score`)), NA, sd(cell_line_bound_data$`Klf4.motif.score`, na.rm = TRUE))
    
    total_num_klf4_motifs <- sum(cell_line_bound_data$`Number.of.Klf4.motifs`, na.rm = TRUE)
    avg_total_motif_score <- ifelse(all(is.na(cell_line_bound_data$`Total.motif.score`)), NA, mean(cell_line_bound_data$`Total.motif.score`, na.rm = TRUE))
    
    
    results <- rbind(results, data.frame(
      Cell_Line = col_name,
      Total_KLF4_Bindings = total_bindings,
      KLF4_Motif_Count = klf4_motif_count,
      KLF4_Motif_Percentage = klf4_motif_percentage,
      Average_Klf4_Motif_Score = avg_klf4_motif_score,
      SD_Klf4_Motif_Score = sd_klf4_motif_score,
      Total_Number_of_Klf4_Motifs = total_num_klf4_motifs,
      Average_Total_Motif_Score = avg_total_motif_score
    ))
  }
  return(results)
}

# Process Data for Human and Mouse
human_results <- process_motif_data(human_data, human_cell_lines, "Human")
mouse_results <- process_motif_data(mouse_data, mouse_cell_lines, "Mouse")

# Bar Graphs for KLF4 Motif Percentage 

# Human Bar Graph
human_bar_plot <- ggplot(human_results, aes(x = Cell_Line, y = KLF4_Motif_Percentage, fill = Cell_Line)) +
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = "Percentage of KLF4 Motif in Human Cell Lines",
    x = "Cell Line",
    y = "Percentage of KLF4 Motif (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  geom_text(aes(label = sprintf("%.2f%%", KLF4_Motif_Percentage)), vjust = -0.5, size = 3) # Percentage labels

print(human_bar_plot)

# Mouse Bar Graph
mouse_bar_plot <- ggplot(mouse_results, aes(x = Cell_Line, y = KLF4_Motif_Percentage, fill = Cell_Line)) +
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = "Percentage of KLF4 Motif in Mouse Cell Lines",
    x = "Cell Line",
    y = "Percentage of KLF4 Motif (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  geom_text(aes(label = sprintf("%.2f%%", KLF4_Motif_Percentage)), vjust = -0.5, size = 3) # Percentage labels

print(mouse_bar_plot)

# Generate Detailed Tables (Charts)

cat("\n--- Human KLF4 Motif Analysis Table ---\n")
human_table <- human_results %>%
  mutate(
    KLF4_Motif_Percentage = sprintf("%.2f%%", KLF4_Motif_Percentage),
    Average_Klf4_Motif_Score = sprintf("%.2f", Average_Klf4_Motif_Score),
    SD_Klf4_Motif_Score = sprintf("%.2f", SD_Klf4_Motif_Score),
    Average_Total_Motif_Score = sprintf("%.2f", Average_Total_Motif_Score)
  )

# Use kable for a nicely formatted table
print(kable(human_table, caption = "Human KLF4 Motif Analysis") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F))


cat("\n--- Mouse KLF4 Motif Analysis Table ---\n")
mouse_table <- mouse_results %>%
  mutate(
    KLF4_Motif_Percentage = sprintf("%.2f%%", KLF4_Motif_Percentage),
    Average_Klf4_Motif_Score = sprintf("%.2f", Average_Klf4_Motif_Score),
    SD_Klf4_Motif_Score = sprintf("%.2f", SD_Klf4_Motif_Score),
    Average_Total_Motif_Score = sprintf("%.2f", Average_Total_Motif_Score)
  )

# Use kable for a nicely formatted table
print(kable(mouse_table, caption = "Mouse KLF4 Motif Analysis") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F))
