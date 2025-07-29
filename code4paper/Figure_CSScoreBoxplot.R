# ==============================================================================
# Figure: CS Score Distribution Boxplot
# Publication-ready R code for generating CS score boxplot visualization
# ==============================================================================

# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(extrafont)

# Load fonts
# font_import()  # Run once if needed
loadfonts(device = "win")

# ==============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load CS score data
csscore <- read_xlsx("csscore.xlsx")

# Rename columns
colnames(csscore) <- c("SimilarityScore", "ConsensusScore", "Evaluation")

# Calculate the minimum of the two scores
csscore$MinScore <- pmin(csscore$SimilarityScore, csscore$ConsensusScore/100)

# Convert Evaluation to factor for better visualization
csscore$EvaluationFactor <- factor(csscore$Evaluation)

# Basic summary statistics
summary(csscore)

# ==============================================================================
# 2. BOXPLOT GENERATION
# ==============================================================================

p3 <- ggplot(csscore, aes(x = EvaluationFactor, y = MinScore*100)) +
  geom_boxplot(fill = "#FF00FF", alpha = 0.7) +
  labs(
    x = "",
    y = "CS Score"
  ) +
  scale_x_discrete(labels = c("0" = "Incorrect", 
                              "0.5" = "Partially Correct", 
                              "1" = "Correct")) +
  scale_y_continuous(limits = c(20, 100)) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times New Roman"),
    axis.text = element_text(size = 20, color = "black", family = "Times New Roman", face = "bold"),
    axis.title = element_text(size = 20, color = "black", family = "Times New Roman", face = "bold"),
    panel.grid = element_blank(),                        # remove grid lines
    panel.background = element_rect(fill = "white", color = NA),  # white plot panel
    plot.background = element_rect(fill = "white", color = NA),   # white around the plot
    axis.line = element_line(color = "black")             # black x and y axis lines
  )

# Save to file
ggsave("cs_score_distribution.jpg", p3, width = 7, height = 6, dpi = 600)
ggsave("cs_score_distribution.pdf", p3, width = 7, height = 6)
ggsave("cs_score_distribution.png", p3, width = 7, height = 6, dpi = 300)

# Show plot
print(p3)

# ==============================================================================
# 3. SOURCE DATA EXPORT
# ==============================================================================

# Prepare the data for export
source_csscore <- csscore %>%
  select(Evaluation, MinScore) %>%
  mutate(
    Evaluation_Label = case_when(
      Evaluation == 0 ~ "Incorrect",
      Evaluation == 0.5 ~ "Partially Correct",
      Evaluation == 1 ~ "Correct",
      TRUE ~ as.character(Evaluation)
    ),
    CS_Score = MinScore * 100  # Convert to percentage scale
  ) %>%
  select(Evaluation, Evaluation_Label, CS_Score) %>%
  arrange(Evaluation)

# Export to CSV
write.csv(source_csscore, "SourceData_Fig_CSScoreBoxplot.csv", row.names = FALSE)

# Optional Excel export
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "CSScore_Boxplot" = source_csscore
  ),
  file = "SourceData_Fig_CSScoreBoxplot.xlsx",
  asTable = TRUE
)

print("CS score boxplot generated successfully!")