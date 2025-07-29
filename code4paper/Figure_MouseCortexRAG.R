# ==============================================================================
# Figure: Mouse Cortex Detailed Annotation with RAG Comparison
# Publication-ready R code for generating grouped bar plot
# ==============================================================================

# Load required libraries
library(ggplot2)

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

# Create data based on updated scores
data <- data.frame(
  CellType = rep(c("Projection Neurons", "Inhibitory Neurons", "Non Neurons"), 6),
  Method = c(
    rep("CASSIA with RAG", 3),
    rep("CASSIA without RAG", 3),
    rep("ScType", 3),
    rep("GPTcelltype4o", 3),
    rep("GPTcelltype4", 3),
    rep("SingleR", 3)
  ),
  Value = c(
    # CASSIA with RAG
    7, 7.71, 9.33,
    # CASSIA without RAG
    4.2, 7.14, 9.33, 
    # ScType
    2.6, 4.29, 7.83, 
    # GPTcelltype4o
    2.2, 5.14, 3, 
    # GPTcelltype4
    0.8, 3.71, 5.17,
    # SingleR
    0.0, 0.0, 0.5
  )
)

# Set factor levels
data$CellType <- factor(data$CellType, levels = c("Projection Neurons", "Inhibitory Neurons", "Non Neurons"))
data$Method <- factor(data$Method, levels = c("CASSIA with RAG", "CASSIA without RAG", "GPTcelltype4", "GPTcelltype4o", "ScType", "SingleR"))

# Convert score to percentage
data$Value <- data$Value * 10

# ==============================================================================
# 2. PLOT GENERATION
# ==============================================================================

p <- ggplot(data, aes(x = CellType, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "CASSIA with RAG" = "lightblue",
    "CASSIA without RAG" = "#2E6C97",
    "ScType" = "#9B95C9",
    "GPTcelltype4o" = "#F17C67",
    "GPTcelltype4" = "#CB1B45",
    "SingleR" = "#6B9CA3"
  )) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  scale_x_discrete(expand = expansion(add = 1)) +
  labs(
    title = "Detailed annotation of mouse cortex",
    x = "Celltype Group",
    y = "Average Score",
    fill = ""
  ) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.5,
      face = "bold", size = 14, family = "Times New Roman", color = "black"
    ),
    axis.text.y = element_text(
      face = "bold", size = 14, family = "Times New Roman", color = "black"
    ),
    axis.title.x = element_text(
      face = "bold", size = 14, family = "Times New Roman",
      margin = margin(t = 20, b = 5)
    ),
    axis.title.y = element_text(
      face = "bold", size = 14, family = "Times New Roman",
      margin = margin(r = 20, l = 5)
    ),
    legend.text = element_text(face = "bold", size = 14, family = "Times New Roman"),
    legend.title = element_text(face = "bold", size = 18, family = "Times New Roman"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 18,
      family = "Times New Roman", margin = margin(b = 20)
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# Show plot
print(p)

# Save as multiple formats
ggsave("mouse_cortex_rag_comparison.pdf", p, width = 10, height = 6, dpi = 300)
ggsave("mouse_cortex_rag_comparison.jpg", p, width = 10, height = 6, dpi = 600)
ggsave("mouse_cortex_rag_comparison.png", p, width = 10, height = 6, dpi = 300)

# ==============================================================================
# 3. SOURCE DATA EXPORT
# ==============================================================================

# Prepare the source data table
source_mouse_cortex <- data %>%
  mutate(
    Average_Score_Percent = round(Value, 2)
  ) %>%
  select(CellType, Method, Average_Score_Percent) %>%
  arrange(CellType, Method)

# Export to CSV
write.csv(source_mouse_cortex, "SourceData_Fig_MouseCortexAnnotation.csv", row.names = FALSE)

# Optional Excel export
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "MouseCortex_Annotation" = source_mouse_cortex
  ),
  file = "SourceData_Fig_MouseCortexAnnotation.xlsx",
  asTable = TRUE
)

# ==============================================================================
# 4. ADDITIONAL RAG PERFORMANCE COMPARISON DATA
# ==============================================================================

# Create data for broader RAG comparison
rag_comparison_data <- data.frame(
  CellType = rep(c("Cortex", "Cerebellum", "Lung"), 6),
  Method = c(rep("CASSIA with RAG", 3), 
            rep("CASSIA without RAG", 3),
            rep("ScType", 3),
            rep("GPTcelltype4o", 3),
            rep("GPTcelltype4", 3),
            rep("SingleR", 3)),
  Value = c(
    # CASSIA has RAG
    7.17, 7.89, 6.95,
    # CASSIA no RAG
    6.04, 6.61, 6.34, 
    # Sctype
    4.48, 4.39, 4.51, 
    # gpt celltype4o
    2.78, 5.38, 4.51, 
    # gpt celltype4
    1.96, 4.22, 4.76,
    # singleR
    0.0, 0.0, 3.41
  ) * 10  # Convert to percentage
)

# Set factor levels
rag_comparison_data$CellType <- factor(rag_comparison_data$CellType, 
                       levels = c("Cortex", "Cerebellum", "Lung"))
rag_comparison_data$Method <- factor(rag_comparison_data$Method, 
                     levels = c("CASSIA with RAG", "CASSIA without RAG", "GPTcelltype4", 
                              "GPTcelltype4o", "ScType", "SingleR"))

# Export RAG comparison source data
source_cassia_rag <- rag_comparison_data %>%
  mutate(
    Value_Percent = round(Value, 1)
  ) %>%
  select(CellType, Method, Value_Percent) %>%
  arrange(CellType, Method)

# Export to CSV
write.csv(source_cassia_rag, "SourceData_Fig_CASSIA_RAG_Performance.csv", row.names = FALSE)

# Export to Excel
openxlsx::write.xlsx(
  list(
    "CASSIA_RAG_Performance" = source_cassia_rag
  ),
  file = "SourceData_Fig_CASSIA_RAG_Performance.xlsx",
  asTable = TRUE
)

print("Mouse cortex RAG comparison figures generated successfully!")