# ==============================================================================
# Figure: Average Annotation Accuracy Heatmap
# Publication-ready R code for generating heatmap visualization
# ==============================================================================

# Load required libraries
library(tidyverse)
library(readxl)

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
  
  # Define the desired column order
  method_cols <- c(
    "CASSIA_correctness",
    "gptcelltype_4_correctness",
    "gptcelltype_4o_correctness",
    "sctype_correctness",
    "singleR_correctness",
    "celltypist_correctness",
    "sccatch_correctness"
  )
  prediction_cols <- c(
    "Predicted.Main.Cell.Type",  # CASSIA's predictions
    "gptcelltype_4",
    "gptcelltype_4o",
    "sctype",
    "singleR",
    "celltypist",
    "sccatch"
  )
  
  # Reorder columns
  other_cols <- setdiff(names(data), 
                       c(method_cols, prediction_cols, 
                         "dataset", "dataset_original", "dataset_grouped"))
  
  data <- data %>%
    select(
      dataset, dataset_original,
      all_of(other_cols),
      all_of(prediction_cols),
      all_of(method_cols)
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
  
  # Convert Marker.Number and Iterations to numeric
  df$Marker.Number <- as.numeric(df$Marker.Number)
  df$Iterations <- as.numeric(df$Iterations)
  
  # Ensure dataset remains as factor
  df$dataset <- as.factor(df$dataset)
  
  # Convert logical columns properly
  logical_cols <- names(df)[sapply(df, is.logical)]
  df[logical_cols] <- lapply(df[logical_cols], as.logical)
  
  return(df)
}

# Apply the cleaning function
result_cleaned <- clean_result_df(result)
df <- result_cleaned

# ==============================================================================
# 2. HEATMAP GENERATION
# ==============================================================================

# Heatmap color gradient based on Nature's red accent
nature_gradient <- colorRampPalette(c("#FFFFFF", "#FFC4C4", "#CB1B45"))(100)

# Prepare heatmap data
heatmap_data <- df %>%
  mutate(across(ends_with("correctness"), as.numeric)) %>%
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
    names_to = "Methods",
    values_to = "Average_Score"
  ) %>%
  mutate(Methods = factor(Methods, 
                         levels = c("CASSIA", "GPTcelltype4", "GPTcelltype4o", "ScType", 
                                  "SingleR", "CellTypist", "scCATCH"))) %>%
  mutate(dataset = factor(dataset, levels = rev(unique(dataset))))

# Create the heatmap plot
p1 <- ggplot(heatmap_data, aes(x = Methods, y = dataset, fill = Average_Score)) +
  geom_tile(color = "#E5E5E5", linewidth = 0.3) +
  geom_text(data = subset(heatmap_data, !is.na(Average_Score)),
            aes(label = sprintf("%.2f", Average_Score)), 
            color = "black",
            size = 7,
            family = "Times New Roman") +
  scale_fill_gradientn(
    colors = nature_gradient,
    limits = c(0, 1),
    na.value = "white",
    name = "Value",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.5,
      barwidth = 1,
      barheight = 10,
      title.position = "top",    # Move title to top
      title.hjust = 0.5         # Center the title
    )
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(family = "Times New Roman", size = 14, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    plot.subtitle = element_text(family = "Times New Roman", size = 11, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_blank(),
    axis.text.x = element_text(family = "Times New Roman", angle = 45, hjust = 1, size = 20, face = "bold", color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 20, face = "bold", color = "black"),
    legend.position = "right",
    legend.title = element_text(family = "Times New Roman", size = 20, face = "bold", margin = margin(b = 25)),  # Add bottom margin to title
    legend.text = element_text(family = "Times New Roman", size = 20),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "Average Annotation Accuracy"
  ) +
  coord_fixed()

# Save the heatmap
ggsave("heatmap_annotation_accuracy.png", p1, width = 12, height = 8)  # in inches

# ==============================================================================
# 3. SOURCE DATA EXPORT
# ==============================================================================

# Create a clean version of heatmap_data for export
source_heatmap <- heatmap_data %>%
  arrange(dataset, Methods) %>%
  rename(
    method = Methods,
    average_score = Average_Score
  )

# Export to CSV
write.csv(source_heatmap, "SourceData_Fig_Heatmap.csv", row.names = FALSE)

# Optional: Save to Excel format (recommended for Nature submissions)
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "AverageAnnotationAccuracy" = source_heatmap
  ),
  file = "SourceData_Fig_Heatmap.xlsx",
  asTable = TRUE
)

print("Heatmap figure generated successfully!")