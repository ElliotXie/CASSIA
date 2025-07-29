# ==============================================================================
# Figure: Fully Accurate Annotation Bar Plot
# Publication-ready R code for showing percentage of fully accurate annotations
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
# 2. CALCULATE FULLY ACCURATE PERCENTAGES
# ==============================================================================

# Calculate percentage of fully accurate (score = 1) annotations for each method
fully_accurate_overall <- df %>%
  summarise(
    CASSIA = mean(CASSIA_correctness == 1, na.rm = TRUE),
    GPTcelltype4 = mean(gptcelltype_4_correctness == 1, na.rm = TRUE),
    GPTcelltype4o = mean(gptcelltype_4o_correctness == 1, na.rm = TRUE),
    ScType = mean(sctype_correctness == 1, na.rm = TRUE),
    SingleR = mean(singleR_correctness == 1, na.rm = TRUE),
    CellTypist = mean(celltypist_correctness == 1, na.rm = TRUE),
    scCATCH = mean(sccatch_correctness == 1, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Method",
    values_to = "FullyAccurateRate"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("CASSIA", "GPTcelltype4", "GPTcelltype4o", 
                                      "ScType", "SingleR", "CellTypist", "scCATCH"))
  )

# Calculate by dataset
fully_accurate_dataset <- df %>%
  group_by(dataset) %>%
  summarise(
    CASSIA = mean(CASSIA_correctness == 1, na.rm = TRUE),
    GPTcelltype4 = mean(gptcelltype_4_correctness == 1, na.rm = TRUE),
    GPTcelltype4o = mean(gptcelltype_4o_correctness == 1, na.rm = TRUE),
    ScType = mean(sctype_correctness == 1, na.rm = TRUE),
    SingleR = mean(singleR_correctness == 1, na.rm = TRUE),
    CellTypist = mean(celltypist_correctness == 1, na.rm = TRUE),
    scCATCH = mean(sccatch_correctness == 1, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -dataset,
    names_to = "Method",
    values_to = "FullyAccurateRate"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("CASSIA", "GPTcelltype4", "GPTcelltype4o", 
                                      "ScType", "SingleR", "CellTypist", "scCATCH"))
  )

# ==============================================================================
# 3. CREATE PLOTS
# ==============================================================================

# Define color palette
method_colors <- c(
  "CASSIA" = "#E64B35",
  "GPTcelltype4" = "#4DBBD5",
  "GPTcelltype4o" = "#00A087",
  "ScType" = "#3C5488",
  "SingleR" = "#F39B7F",
  "CellTypist" = "#8491B4",
  "scCATCH" = "#91D1C2"
)

# Overall fully accurate rate plot
p_fully_overall <- ggplot(fully_accurate_overall, aes(x = Method, y = FullyAccurateRate, fill = Method)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", FullyAccurateRate * 100)), 
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
    title = "Percentage of Fully Accurate Annotations (Score = 1)",
    x = "Method",
    y = "Fully Accurate Rate"
  )

# Save overall plot
ggsave("fully_accurate_overall.jpg", p_fully_overall, width = 10, height = 8, dpi = 600)
ggsave("fully_accurate_overall.pdf", p_fully_overall, width = 10, height = 8)
ggsave("fully_accurate_overall.png", p_fully_overall, width = 10, height = 8, dpi = 300)

# Dataset-specific fully accurate rate plot
p_fully_dataset <- ggplot(fully_accurate_dataset, aes(x = dataset, y = FullyAccurateRate, fill = Method)) +
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
    title = "Fully Accurate Annotations by Dataset",
    x = "Dataset",
    y = "Fully Accurate Rate"
  )

# Save dataset plot
ggsave("fully_accurate_dataset.jpg", p_fully_dataset, width = 12, height = 8, dpi = 600)
ggsave("fully_accurate_dataset.pdf", p_fully_dataset, width = 12, height = 8)
ggsave("fully_accurate_dataset.png", p_fully_dataset, width = 12, height = 8, dpi = 300)

# ==============================================================================
# 4. CREATE STACKED BAR PLOT FOR ACCURACY LEVELS
# ==============================================================================

# Calculate distribution of accuracy levels (0, 0.5, 1)
accuracy_distribution <- df %>%
  pivot_longer(
    cols = ends_with("_correctness"),
    names_to = "Method",
    values_to = "Score"
  ) %>%
  mutate(
    Method = str_remove(Method, "_correctness"),
    Method = case_when(
      Method == "CASSIA" ~ "CASSIA",
      Method == "gptcelltype_4" ~ "GPTcelltype4",
      Method == "gptcelltype_4o" ~ "GPTcelltype4o",
      Method == "sctype" ~ "ScType",
      Method == "singleR" ~ "SingleR",
      Method == "celltypist" ~ "CellTypist",
      Method == "sccatch" ~ "scCATCH",
      TRUE ~ Method
    ),
    Method = factor(Method, levels = c("CASSIA", "GPTcelltype4", "GPTcelltype4o", 
                                      "ScType", "SingleR", "CellTypist", "scCATCH")),
    ScoreCategory = case_when(
      Score == 0 ~ "Incorrect (0)",
      Score == 0.5 ~ "Partially Correct (0.5)",
      Score == 1 ~ "Fully Correct (1)",
      TRUE ~ "Missing"
    ),
    ScoreCategory = factor(ScoreCategory, 
                          levels = c("Fully Correct (1)", "Partially Correct (0.5)", 
                                    "Incorrect (0)", "Missing"))
  ) %>%
  filter(!is.na(Score)) %>%
  group_by(Method, ScoreCategory) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(Percentage = Count / sum(Count))

# Create stacked bar plot
p_stacked <- ggplot(accuracy_distribution, aes(x = Method, y = Percentage, fill = ScoreCategory)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(
    "Fully Correct (1)" = "#2E8B57",      # Sea green
    "Partially Correct (0.5)" = "#FFA500", # Orange
    "Incorrect (0)" = "#DC143C",          # Crimson
    "Missing" = "#808080"                 # Gray
  )) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Distribution of Annotation Accuracy Scores",
    x = "Method",
    y = "Percentage",
    fill = "Accuracy Level"
  )

# Save stacked plot
ggsave("accuracy_distribution_stacked.jpg", p_stacked, width = 10, height = 8, dpi = 600)
ggsave("accuracy_distribution_stacked.pdf", p_stacked, width = 10, height = 8)
ggsave("accuracy_distribution_stacked.png", p_stacked, width = 10, height = 8, dpi = 300)

# ==============================================================================
# 5. SOURCE DATA EXPORT
# ==============================================================================

# Export all data
write.csv(fully_accurate_overall, "SourceData_Fig_FullyAccurateOverall.csv", row.names = FALSE)
write.csv(fully_accurate_dataset, "SourceData_Fig_FullyAccurateDataset.csv", row.names = FALSE)
write.csv(accuracy_distribution, "SourceData_Fig_AccuracyDistribution.csv", row.names = FALSE)

# Optional Excel export
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "Fully_Accurate_Overall" = fully_accurate_overall,
    "Fully_Accurate_Dataset" = fully_accurate_dataset,
    "Accuracy_Distribution" = accuracy_distribution
  ),
  file = "SourceData_Fig_FullyAccurateAnalysis.xlsx",
  asTable = TRUE
)

print("Fully accurate annotation plots generated successfully!")