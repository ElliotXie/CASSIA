# Simple CASSIA R Package Test Script
# This script tests core functionality and stores results in a folder

# Create results directory
results_dir <- "test_results_simple"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Set working directory to results folder for outputs
original_wd <- getwd()
setwd(results_dir)

cat("=== Simple CASSIA R Package Test ===\n")
cat("Results stored in:", file.path(original_wd, results_dir), "\n")
cat("Test started at:", as.character(Sys.time()), "\n\n")

# Load library and set environment
library(CASSIA)
reticulate::use_condaenv("cassia_env")

# Load marker data
cat("Loading marker data...\n")
unprocessed_markers <- CASSIA::loadBuiltinMarkers(marker_type = "unprocessed")
subcluster_markers_data <- CASSIA::loadBuiltinMarkers(marker_type = "subcluster_results")

cat("✓ Loaded", nrow(unprocessed_markers), "unprocessed markers\n")
cat("✓ Loaded", nrow(subcluster_markers_data), "subcluster markers\n")

# Test parameters
tissue_param <- "large intestine"
species_param <- "human"
model_param <- "google/gemini-2.5-flash-preview-05-20"
provider_param <- "openrouter"

# --- Test 1: Basic Batch Analysis ---
cat("\n--- Test 1: Basic Batch Analysis ---\n")
batch_output <- "simple_batch_test"

tryCatch({
  start_time <- Sys.time()
  CASSIA::runCASSIA_batch(
    marker = unprocessed_markers,
    output_name = batch_output,
    model = model_param,
    tissue = tissue_param,
    species = species_param,
    max_workers = 3,  # Reduced for stability
    n_genes = 50,
    provider = provider_param
  )
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat("✓ Batch analysis completed in", round(duration, 1), "seconds\n")
}, error = function(e) {
  cat("✗ Error in batch analysis:", e$message, "\n")
})

# --- Test 2: Quality Scoring ---
cat("\n--- Test 2: Quality Scoring ---\n")
input_file <- paste0(batch_output, "_full.csv")
scored_file <- paste0(batch_output, "_scored.csv")

tryCatch({
  if (file.exists(input_file)) {
    CASSIA::runCASSIA_score_batch(
      input_file = input_file,
      output_file = scored_file,
      max_workers = 3,
      model = "deepseek/deepseek-chat-v3-0324",
      provider = provider_param
    )
    cat("✓ Quality scoring completed\n")
    
    # Generate report
    CASSIA::runCASSIA_generate_score_report(
      csv_path = scored_file,
      output_name = paste0(batch_output, "_report.html")
    )
    cat("✓ Score report generated\n")
  } else {
    cat("⚠ Skipping scoring - input file not found\n")
  }
}, error = function(e) {
  cat("✗ Error in quality scoring:", e$message, "\n")
})

# --- Test 3: Cell Type Comparison ---
cat("\n--- Test 3: Cell Type Comparison ---\n")
marker_string <- "IGLL5, IGLV6-57, JCHAIN, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1"

tryCatch({
  CASSIA::compareCelltypes(
    tissue = tissue_param,
    celltypes = c("Plasma Cells", "IgA-secreting Plasma Cells"),
    marker = marker_string,
    species = species_param,
    output_file = "cell_comparison_test.html"
  )
  cat("✓ Cell type comparison completed\n")
}, error = function(e) {
  cat("✗ Error in cell type comparison:", e$message, "\n")
})

# --- Test 4: Annotation Boost ---
cat("\n--- Test 4: Annotation Boost ---\n")
tryCatch({
  if (file.exists(input_file)) {
    CASSIA::runCASSIA_annotationboost(
      full_result_path = input_file,
      marker = unprocessed_markers,
      output_name = "annotation_boost_test",
      cluster_name = "monocyte",
      major_cluster_info = "human large intestine",
      num_iterations = 3,  # Reduced for testing
      model = model_param,
      provider = provider_param
    )
    cat("✓ Annotation boost completed\n")
  } else {
    cat("⚠ Skipping annotation boost - input file not found\n")
  }
}, error = function(e) {
  cat("✗ Error in annotation boost:", e$message, "\n")
})

# --- Test 5: Subclustering ---
cat("\n--- Test 5: Subclustering Analysis ---\n")
tryCatch({
  CASSIA::runCASSIA_subclusters(
    marker = subcluster_markers_data,
    major_cluster_info = "cd8 t cell",
    output_name = "subcluster_test",
    model = model_param,
    provider = provider_param,
    n_genes = 50L
  )
  cat("✓ Subclustering analysis completed\n")
}, error = function(e) {
  cat("✗ Error in subclustering:", e$message, "\n")
})

# --- Test Summary ---
cat("\n--- Test Summary ---\n")
all_files <- list.files(".", recursive = TRUE)
cat("Generated files:\n")
for(file in all_files) {
  cat("  -", file, "\n")
}

# Create summary
summary_data <- data.frame(
  Test = c("Batch Analysis", "Quality Scoring", "Cell Type Comparison", 
           "Annotation Boost", "Subclustering"),
  Status = c("Tested", "Tested", "Tested", "Tested", "Tested"),
  Files = c("CSV files", "HTML reports", "HTML comparison", 
            "HTML/TXT", "CSV files")
)

write.csv(summary_data, "test_summary.csv", row.names = FALSE)
cat("✓ Test summary saved\n")

# Return to original directory
setwd(original_wd)

cat("\n=== Simple CASSIA test completed ===\n")
cat("Check the", results_dir, "folder for all outputs\n") 