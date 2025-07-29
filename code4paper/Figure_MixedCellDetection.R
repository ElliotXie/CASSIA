# ==============================================================================
# Figure: Mixed Cell Detection Accuracy
# Publication-ready R code for generating bar plots by marker number
# ==============================================================================

# Load required libraries
library(tidyverse)
library(scales)
library(tidyr)
library(showtext)
library(Cairo)

# Add Times New Roman font
showtext_auto()
font_add("Times New Roman", regular = "times.ttf", bold = "timesbd.ttf")  # Path may need adjustment

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

# Define method order and colors
method_order <- c("CASSIA (default+scoring)", "CASSIA (default)")
ratio_labels <- c("50_50" = "50:50 Mix", "80_20" = "80:20 Mix", 
                  "40_60" = "40:60 Mix", "30_70" = "30:70 Mix")
ratio_colors <- c("50:50 Mix" = "#3C5488",
                  "80:20 Mix" = "#E64B35",
                  "40:60 Mix" = "#00A087",
                  "30:70 Mix" = "#4DBBD5")

# Create dataset
mix_dat <- tibble::tribble(
  ~ratio,   ~marker_number, ~`CASSIA (default)`, ~`CASSIA (default+scoring)`,
  "30_70",  30,             0.65,               0.8,
  "40_60",  30,             0.65,               1.0,
  "50_50",  30,             0.8,                0.9,
  "80_20",  30,             0.5,                0.75,
  "30_70",  50,             0.55,               0.85,
  "50_50",  50,             0.6,                0.9,
  "40_60",  50,             0.7,                1.0,
  "80_20",  50,             0.55,               0.8,
  "30_70",  70,             0.55,               0.85,
  "40_60",  70,             0.65,               0.95,
  "50_50",  70,             0.8,                0.8,
  "80_20",  70,             0.65,               0.85
)

# Transform data to long format
long_dat <- mix_dat %>%
  pivot_longer(cols = starts_with("CASSIA"),
               names_to = "Method", values_to = "Value") %>%
  mutate(
    Method = factor(Method, levels = method_order),
    Scenario = recode(ratio, !!!ratio_labels)
  )

# Complete all combinations
all_complete <- long_dat %>%
  complete(
    marker_number = unique(long_dat$marker_number),
    Scenario = unique(long_dat$Scenario),
    Method = unique(long_dat$Method),
    fill = list(Value = 0)
  )

# ==============================================================================
# 2. PLOT THEME DEFINITION
# ==============================================================================

theme_publication <- theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", family = "Times New Roman", color = "black"),
    axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1, family = "Times New Roman", color = "black"),
    axis.title = element_text(size = 20, face = "bold", family = "Times New Roman", color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, family = "Times New Roman", color = "black"),
    legend.text = element_text(size = 12, face = "bold", family = "Times New Roman", color = "black"),
    legend.title = element_text(size = 12, face = "bold", family = "Times New Roman", color = "black"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  )

# ==============================================================================
# 3. GENERATE PLOTS FOR EACH MARKER NUMBER
# ==============================================================================

for (mk in sort(unique(all_complete$marker_number))) {
  sub_dat <- filter(all_complete, marker_number == mk)
  
  p <- ggplot(sub_dat, aes(x = Method, y = Value, fill = Scenario)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = percent(Value, accuracy = 1)),
              position = position_dodge(0.7),
              vjust = -0.5,
              family = "Times New Roman",
              fontface = "bold",
              size = 5) +
    scale_fill_manual(values = ratio_colors, name = "Mix ratio") +
    scale_y_continuous(limits = c(0, 1.1),
                       breaks = seq(0, 1, 0.2),
                       labels = percent_format(accuracy = 1)) +
    labs(
      title = paste("Mixed Cell Detection Accuracy (Markers =", mk, ")"),
      x = NULL,
      y = "Detection Rate"
    ) +
    theme_publication
  
  # Save in multiple formats
  ggsave(paste0("mixed_cell_detection_marker_", mk, ".pdf"), p, width = 10, height = 6, device = cairo_pdf)
  ggsave(paste0("mixed_cell_detection_marker_", mk, ".jpg"), p, width = 10, height = 6, dpi = 600)
  ggsave(paste0("mixed_cell_detection_marker_", mk, ".png"), p, width = 10, height = 6, dpi = 300, device = cairo_png)
}

# ==============================================================================
# 4. SOURCE DATA EXPORT
# ==============================================================================

# Format source data table
source_detect_marker <- all_complete %>%
  mutate(
    detection_rate_percent = round(Value * 100, 1)
  ) %>%
  select(marker_number, Scenario, Method, Value, detection_rate_percent) %>%
  arrange(marker_number, Scenario, Method)

# Export to CSV
write.csv(source_detect_marker, "SourceData_Fig_MixedCellDetection.csv", row.names = FALSE)

# Optional Excel export
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "Detection_By_Marker" = source_detect_marker
  ),
  file = "SourceData_Fig_MixedCellDetection.xlsx",
  asTable = TRUE
)

print("Mixed cell detection figures generated successfully!")