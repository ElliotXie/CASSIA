# ==============================================================================
# Figure: Cancer Detection Accuracy Heatmap
# Publication-ready R code for generating cancer detection heatmap
# ==============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(extrafont)

# Load fonts
# font_import(prompt = FALSE)  # Run once if needed
loadfonts(device = "win")

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

# Define data
df <- tibble(
  dataset = c("Brain_met1", "Brain_met2", "Brain_met3", "Lung_pri", "Colon_pri"),
  CASSIA_correctness = c(1, 8/12, 12/17, 4/7, 2/4),
  Cassia_enhanced_correctness = c(1, 15/17, 1, 1, 1),
  gptcelltype_4_correctness = c(2/5, 3/12, 3/17, 1/7, 0),
  gptcelltype_4o_correctness = c(1/5, 0, 0, 1/7, 0),
  sctype_correctness = c(0, 1/12, 2/17, 1/7, 0)
)

# Method order and labels
method_order <- c("CASSIA", "Cassia_enhanced", "gptcelltype_4", "gptcelltype_4o", "sctype")
method_labels <- c("CASSIA", "CASSIA Enhanced", "GPTcelltype4", "GPTcelltype4o", "ScType")

# Reshape data
heatmap_data <- df %>%
  pivot_longer(
    cols = ends_with("_correctness"),
    names_to = "method",
    values_to = "fully_correct_rate"
  ) %>%
  mutate(
    method = str_remove(method, "_correctness"),
    method = factor(method, levels = method_order, labels = method_labels),
    fully_correct_rate = replace_na(fully_correct_rate, 0),
    label = scales::percent(fully_correct_rate, accuracy = 1),
    dataset = factor(dataset, levels = c("Brain_met1", "Brain_met2", "Brain_met3", "Lung_pri", "Colon_pri"))
  )

# ==============================================================================
# 2. HEATMAP GENERATION
# ==============================================================================

# Define Nature-style gradient
nature_gradient <- colorRampPalette(c("#FFFFFF", "#FFC4C4", "#CB1B45"))(100)

# Plot heatmap
p <- ggplot(heatmap_data, aes(x = method, y = dataset, fill = fully_correct_rate)) +
  geom_tile(color = NA, linewidth = 0) +  # No gridlines between tiles
  geom_text(aes(label = label), size = 5, family = "Times New Roman", fontface = "bold", color = "black") +
  scale_fill_gradientn(colors = nature_gradient, name = "Accuracy") +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title = element_blank(),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

# Save to multiple formats
ggsave("cancer_detection_heatmap.jpg", plot = p, width = 8, height = 5, dpi = 600)
ggsave("cancer_detection_heatmap.pdf", plot = p, width = 8, height = 5)
ggsave("cancer_detection_heatmap.png", plot = p, width = 8, height = 5, dpi = 300)

# ==============================================================================
# 3. SOURCE DATA EXPORT
# ==============================================================================

# Prepare source data
source_cancer <- heatmap_data %>%
  select(dataset, method, fully_correct_rate) %>%
  mutate(
    accuracy_percent = round(fully_correct_rate * 100, 1)
  ) %>%
  arrange(dataset, method)

# Export to CSV
write.csv(source_cancer, "SourceData_Fig_CancerDetection.csv", row.names = FALSE)

# Optional Excel export
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "Cancer_Detection" = source_cancer
  ),
  file = "SourceData_Fig_CancerDetection.xlsx",
  asTable = TRUE
)

print("Cancer detection heatmap generated successfully!")