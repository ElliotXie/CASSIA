#!/usr/bin/env Rscript

# Example script demonstrating how to use the merge_annotations function
library(CASSIA)

# Create a sample CSV with test data
create_test_csv <- function(filepath) {
  data <- data.frame(
    "Cluster ID" = 1:8,
    "Predicted General Cell Type" = c(
      "macrophage",
      "CD4 T cell",
      "B cell",
      "dendritic cell",
      "CD8 T cell",
      "NK cell",
      "epithelial cell",
      "fibroblast"
    ),
    "Predicted Detailed Cell Type" = c(
      "inflammatory macrophage, resident macrophage", 
      "naive CD4 T cell, memory CD4 T cell", 
      "memory B cell, plasma cell", 
      "plasmacytoid dendritic cell, conventional dendritic cell", 
      "cytotoxic CD8 T cell, exhausted CD8 T cell", 
      "CD56bright NK cell, CD56dim NK cell", 
      "type II pneumocyte, type I pneumocyte", 
      "activated fibroblast, quiescent fibroblast"
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(data, filepath, row.names = FALSE)
  cat(sprintf("Created test CSV file at: %s\n", filepath))
  return(data)
}

# First, check if Python dependencies are installed
cat("Checking Python dependencies...\n")
deps <- test_python_dependencies(install = TRUE)
print(deps)

# Example 1: Using a single detail level
run_single_level_example <- function() {
  # Create test data
  test_csv <- "test_clusters.csv"
  create_test_csv(test_csv)
  
  # Additional context to help with annotations
  additional_context <- "
  Cell type reference information:
  - Macrophages and dendritic cells are types of myeloid cells
  - CD4 and CD8 T cells belong to T lymphocyte lineage
  - B cells are part of the B lymphocyte lineage
  - NK cells are natural killer cells and considered part of innate lymphoid cells
  - Epithelial cells form the tissue lining organs
  - Fibroblasts are connective tissue cells
  "
  
  # Set API key for OpenRouter (example - replace with your actual key)
  # set_llm_api_key("your_api_key_here", provider = "openrouter")
  
  # Process with a single detail level
  cat("\n=== Example 1: Single Detail Level (detailed) ===\n")
  result_df <- merge_annotations(
    csv_path = test_csv,
    output_path = "detailed_groupings.csv",
    provider = "openrouter",
    model = "deepseek/deepseek-chat-v3-0324",
    additional_context = additional_context,
    batch_size = 8,
    detail_level = "detailed"
  )
  
  # Display results
  cat("\nResults with detailed grouping level:\n")
  print(head(result_df))
  
  # Clean up
  unlink(test_csv)
  unlink("detailed_groupings.csv")
}

# Example 2: Process all levels in parallel
run_parallel_example <- function() {
  # Create test data
  test_csv <- "test_clusters.csv"
  create_test_csv(test_csv)
  
  # Additional context to help with annotations
  additional_context <- "
  Cell type reference information:
  - Macrophages and dendritic cells are types of myeloid cells
  - CD4 and CD8 T cells belong to T lymphocyte lineage
  - B cells are part of the B lymphocyte lineage
  - NK cells are natural killer cells and considered part of innate lymphoid cells
  - Epithelial cells form the tissue lining organs
  - Fibroblasts are connective tissue cells
  "
  
  # Process all detail levels in parallel
  cat("\n=== Example 2: All Detail Levels in Parallel ===\n")
  all_results_df <- merge_annotations(
    csv_path = test_csv,
    output_path = "all_groupings.csv",
    provider = "openrouter",
    model = "deepseek/deepseek-chat-v3-0324",
    additional_context = additional_context,
    batch_size = 8,
    process_all = TRUE
  )
  
  # Display results
  cat("\nResults from parallel processing (all three detail levels):\n")
  print(head(all_results_df))
  
  # Clean up
  unlink(test_csv)
  unlink("all_groupings.csv")
}

# Uncomment to run the examples
# run_single_level_example()
# run_parallel_example()

cat("\nPlease set your API key before running the examples:\n")
cat("set_llm_api_key('your_api_key_here', provider = 'openrouter')\n")
cat("Then uncomment and run one of the example functions\n") 