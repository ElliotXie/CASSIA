# ==============================================================================
# Figure: Accuracy of Majority Voting by Number of Runs
# Publication-ready R code for generating majority voting accuracy curves
# ==============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)

# ==============================================================================
# 1. FUNCTIONS FOR CALCULATION
# ==============================================================================

# Function to calculate binomial probability
binom_prob <- function(n, k, p) {
  choose(n, k) * p^k * (1-p)^(n-k)
}

# Function to calculate probability of at least k successes in n trials
prob_at_least_k <- function(n, k, p) {
  sum(sapply(k:n, function(i) binom_prob(n, i, p)))
}

# Calculate accuracy for different number of runs
calculate_accuracies <- function(min_runs, max_runs, base_accuracy) {
  results <- data.frame(
    runs = min_runs:max_runs,
    accuracy = NA
  )
  
  for (i in 1:nrow(results)) {
    runs <- results$runs[i]
    # Standard majority threshold (more than half)
    majority_threshold <- ceiling(runs / 2)
    
    # Calculate accuracy with standard majority threshold
    results$accuracy[i] <- prob_at_least_k(runs, majority_threshold, base_accuracy)
  }
  
  return(results)
}

# ==============================================================================
# 2. DATA CALCULATION
# ==============================================================================

# Calculate results for p=0.6, p=0.7, and p=0.8
min_runs <- 5
max_runs <- 100
results_60 <- calculate_accuracies(min_runs, max_runs, 0.6)
results_70 <- calculate_accuracies(min_runs, max_runs, 0.7)
results_80 <- calculate_accuracies(min_runs, max_runs, 0.8)

# Combine results into a single data frame
results_60$base_accuracy <- "60%"
results_70$base_accuracy <- "70%"
results_80$base_accuracy <- "80%"
all_results <- rbind(results_60, results_70, results_80)

# Create a subset of points for a cleaner plot
# Keep all points from 5-20 runs, then every 5th point
subset_runs <- c(1:20, seq(25, 100, by = 5))
plot_data <- all_results %>% filter(runs %in% subset_runs)

# Create base accuracy reference data
base_acc_df <- data.frame(
  base_accuracy = c("60%", "70%", "80%"),
  value = c(0.6, 0.7, 0.8)
)

# ==============================================================================
# 3. PLOT GENERATION
# ==============================================================================

p <- ggplot() +
  # Add the accuracy lines
  geom_line(data = plot_data, 
            aes(x = runs, y = accuracy, color = base_accuracy, group = base_accuracy), 
            size = 1) +
  geom_point(data = plot_data, 
             aes(x = runs, y = accuracy, color = base_accuracy), 
             size = 2) +
  # Add horizontal lines for base accuracies
  geom_hline(data = base_acc_df, 
             aes(yintercept = value, color = base_accuracy), 
             linetype = "dashed", 
             alpha = 0.7) +
  # Set color scheme
  scale_color_manual(values = c("60%" = "#8884d8", "70%" = "#ff7300", "80%" = "#82ca9d"),
                     name = "Base Accuracy") +
  # Set axis labels and limits
  scale_y_continuous(labels = scales::percent, limits = c(0.5, 1)) +
  labs(title = "Accuracy of Majority Voting by Number of Runs",
       subtitle = "Comparing models with 60%, 70%, and 80% base accuracy",
       x = "Number of Runs",
       y = "Accuracy") +
  # Add theme elements
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

# Display the plot
print(p)

# Save the plot
ggsave("majority_voting_accuracy.jpg", p, width = 10, height = 7, dpi = 300)
ggsave("majority_voting_accuracy.pdf", p, width = 10, height = 7)
ggsave("majority_voting_accuracy.png", p, width = 10, height = 7, dpi = 300)

# ==============================================================================
# 4. SOURCE DATA EXPORT
# ==============================================================================

# Export full accuracy data
source_majority_all <- all_results %>%
  arrange(base_accuracy, runs) %>%
  mutate(accuracy_percent = round(accuracy * 100, 2))

# Export plot subset (cleaned version for selected points)
source_majority_plot <- plot_data %>%
  arrange(base_accuracy, runs) %>%
  mutate(accuracy_percent = round(accuracy * 100, 2))

# Write CSVs
write.csv(source_majority_all, "SourceData_Fig_MajorityVoting_All.csv", row.names = FALSE)
write.csv(source_majority_plot, "SourceData_Fig_MajorityVoting_PlotSubset.csv", row.names = FALSE)

# Write Excel export
library(openxlsx)

openxlsx::write.xlsx(
  list(
    "All_Runs" = source_majority_all,
    "Plot_Subset" = source_majority_plot
  ),
  file = "SourceData_Fig_MajorityVoting.xlsx",
  asTable = TRUE
)

print("Majority voting accuracy figure generated successfully!")