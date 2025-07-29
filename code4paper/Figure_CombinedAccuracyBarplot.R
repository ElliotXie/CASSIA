# ==============================================================================
# Figure: Combined Accuracy Bar Plot
# Publication-ready R code for comparing all methods across all datasets
# ==============================================================================

# Load required libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(scales)

# ==============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ==============================================================================

# Specify the path to your data file
file_path <- "Supplementary_Data_4.xlsx"

# Read the first sheet of the xlsx file
data <- read_excel(file_path, sheet = 1)

# Function to reorder and group datasets
reorder_and_group_datasets <- function(data) {
  # Define the dataset ordering and grouping
  dataset_mapping <- c(
    # Core datasets (first five)
    "GTEx" = "1_GTEX",
    "HCL" = "2_HCL",
    "TS" = "3_TS",
    "MCA" = "4_MCA",
    "Azimuth" = "5_Azimuth",
    
    # Cancer datasets
    "Colon Cancer" = "6_Cancer",
    "Lympho node met" = "6_Cancer",
    "Brain met2" = "6_Cancer",
    "CVS" = "6_Cancer",
    "non small cell lung" = "6_Cancer",
    
    # Immune datasets
    "PBMC68k" = "7_Immune",
    "ProjectTILs" = "7_Immune",
    
    # Rare datasets
    "Shark" = "8_Rare",
    "Non-model mammal" = "8_Rare"
  )
  
  # First, group and reorder datasets
  data <- data %>%
    mutate(
      dataset_grouped = factor(
        dataset_mapping[dataset],
        levels = unique(dataset_mapping),
        labels = c("GTEx", "HCL", "TS", "MCA", "Azimuth",
                  "Cancer", "Immune", "Rare")
      ),
      dataset_original = dataset,
      dataset = dataset_grouped
    )
  
  return(data)
}

# Apply reordering function
result <- reorder_and_group_datasets(data)

# Function to clean NA values and convert types
clean_result_df <- function(df) {
  # First, convert all string "NA" to real NA
  df[] <- lapply(df, function(x) {
    if(is.character(x)) {
      # Replace variations of "NA" strings with real NA
      x[x %in% c("NA", "na", "N/A", "n/a")] <- NA
    }
    return(x)
  })
  
  # Convert correctness columns to numeric
  correctness_cols <- c(
    "CASSIA_correctness",
    "gptcelltype_4_correctness",
    "gptcelltype_4o_correctness",
    "sctype_correctness",
    "singleR_correctness",
    "celltypist_correctness",
    "sccatch_correctness"
  )
  
  # Convert correctness columns to numeric
  df[correctness_cols] <- lapply(df[correctness_cols], function(x) {
    as.numeric(x)
  })
  
  # Ensure dataset remains as factor
  df$dataset <- as.factor(df$dataset)
  
  return(df)
}

# Apply the cleaning function
result_cleaned <- clean_result_df(result)
df <- result_cleaned

# ==============================================================================
# 2. COMBINED ACCURACY BAR PLOT
# ==============================================================================

# Calculate overall accuracy for each method
overall_accuracy <- df %>%
  summarise(
    CASSIA = mean(as.numeric(CASSIA_correctness), na.rm = TRUE),
    GPTcelltype4 = mean(as.numeric(gptcelltype_4_correctness), na.rm = TRUE),
    GPTcelltype4o = mean(as.numeric(gptcelltype_4o_correctness), na.rm = TRUE),
    ScType = mean(as.numeric(sctype_correctness), na.rm = TRUE),
    SingleR = mean(as.numeric(singleR_correctness), na.rm = TRUE),
    CellTypist = mean(as.numeric(celltypist_correctness), na.rm = TRUE),
    scCATCH = mean(as.numeric(sccatch_correctness), na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Method",
    values_to = "Accuracy"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("CASSIA", "GPTcelltype4", "GPTcelltype4o", 
                                      "ScType", "SingleR", "CellTypist", "scCATCH"))
  )

# Create color palette
method_colors <- c(
  "CASSIA" = "#E64B35",
  "GPTcelltype4" = "#4DBBD5",
  "GPTcelltype4o" = "#00A087",
  "ScType" = "#3C5488",
  "SingleR" = "#F39B7F",
  "CellTypist" = "#8491B4",
  "scCATCH" = "#91D1C2"
)

# Create the combined accuracy bar plot
p_combined <- ggplot(overall_accuracy, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Accuracy)), 
            vjust = -0.5, 
            size = 4,
            fontface = "bold") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Overall Annotation Accuracy Across All Datasets",
    x = "Method",
    y = "Average Accuracy"
  )

# Save the plot
ggsave("combined_accuracy_barplot.jpg", p_combined, width = 10, height = 8, dpi = 600)
ggsave("combined_accuracy_barplot.pdf", p_combined, width = 10, height = 8)
ggsave("combined_accuracy_barplot.png", p_combined, width = 10, height = 8, dpi = 300)

# ==============================================================================
# 3. DATASET-SPECIFIC ACCURACY BAR PLOT
# ==============================================================================

# Calculate accuracy by dataset
dataset_accuracy <- df %>%
  group_by(dataset) %>%
  summarise(
    CASSIA = mean(as.numeric(CASSIA_correctness), na.rm = TRUE),
    GPTcelltype4 = mean(as.numeric(gptcelltype_4_correctness), na.rm = TRUE),
    GPTcelltype4o = mean(as.numeric(gptcelltype_4o_correctness), na.rm = TRUE),
    ScType = mean(as.numeric(sctype_correctness), na.rm = TRUE),
    SingleR = mean(as.numeric(singleR_correctness), na.rm = TRUE),
    CellTypist = mean(as.numeric(celltypist_correctness), na.rm = TRUE),
    scCATCH = mean(as.numeric(sccatch_correctness), na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -dataset,
    names_to = "Method",
    values_to = "Accuracy"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("CASSIA", "GPTcelltype4", "GPTcelltype4o", 
                                      "ScType", "SingleR", "CellTypist", "scCATCH"))
  )

# Create grouped bar plot by dataset
p_dataset <- ggplot(dataset_accuracy, aes(x = dataset, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Annotation Accuracy by Dataset",
    x = "Dataset",
    y = "Average Accuracy"
  )

# Save the grouped plot
ggsave("dataset_accuracy_barplot.jpg", p_dataset, width = 12, height = 8, dpi = 600)
ggsave("dataset_accuracy_barplot.pdf", p_dataset, width = 12, height = 8)
ggsave("dataset_accuracy_barplot.png", p_dataset, width = 12, height = 8, dpi = 300)

# ==============================================================================
# 4. SOURCE DATA EXPORT
# ==============================================================================

# Export overall accuracy data
write.csv(overall_accuracy, "SourceData_Fig_CombinedAccuracy.csv", row.names = FALSE)

# Export dataset-specific accuracy data
write.csv(dataset_accuracy, "SourceData_Fig_DatasetAccuracy.csv", row.names = FALSE)

# Optional Excel export
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "Overall_Accuracy" = overall_accuracy,
    "Dataset_Accuracy" = dataset_accuracy
  ),
  file = "SourceData_Fig_AccuracyBarplots.xlsx",
  asTable = TRUE
)

print("Combined accuracy bar plots generated successfully!")