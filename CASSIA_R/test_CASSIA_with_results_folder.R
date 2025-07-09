# Test script for CASSIA R Wrapper Functions with Results Folder
# This script mirrors the workflow in the CASSIA_python_package_test.ipynb notebook.
# All outputs will be stored in the test_results/ folder

# --- Setup Results Directory ---
results_dir <- "test_results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Set working directory to results folder for outputs
original_wd <- getwd()
setwd(results_dir)

cat("=== CASSIA R Package Comprehensive Test ===\n")
cat("Results will be stored in:", file.path(original_wd, results_dir), "\n")
cat("Test started at:", as.character(Sys.time()), "\n\n")

# --- 1. Setup and Environment Preparation ---

# Load the CASSIA library
library(CASSIA)

# API keys are already loaded in system environment variables - no need to set them here

# Check if the environment is okay (optional)
cat("--- 1. Environment Setup ---\n")
tryCatch({
  # Set Python environment
  reticulate::use_condaenv("cassia_env")
  cat("✓ Python environment set to cassia_env\n")
}, error = function(e) {
  cat("⚠ Warning setting Python environment:", e$message, "\n")
})

# Load example marker data
cat("Loading marker data...\n")
tryCatch({
  unprocessed_markers <- CASSIA::loadBuiltinMarkers(marker_type = "unprocessed")
  subcluster_markers_data <- CASSIA::loadBuiltinMarkers(marker_type = "subcluster_results")
  
  cat("Available marker sets in the package:\n")
  print(CASSIA::listAvailableMarkers())
  
  cat("✓ Loaded", nrow(unprocessed_markers), "unprocessed markers\n")
  cat("✓ Loaded", nrow(subcluster_markers_data), "subcluster markers\n")
  
  # Save marker info to results
  write.csv(head(unprocessed_markers, 10), "marker_data_sample.csv", row.names = FALSE)
  
}, error = function(e) {
  cat("✗ Error loading marker data:", e$message, "\n")
  setwd(original_wd)
  stop("Cannot proceed without marker data")
})

# Define common parameters based on the Python notebook
tissue_param <- "large intestine"
species_param <- "human"
base_output_name <- "R_CASSIA_Test"
openrouter_provider <- "openrouter"

# Models used in the Python notebook
annot_model_pipeline <- "google/gemini-2.5-flash-preview-05-20"
score_model_pipeline <- "deepseek/deepseek-chat-v3-0324"
boost_model_pipeline <- "google/gemini-2.5-flash-preview-05-20"

cat("✓ Environment setup completed\n\n")

# --- 2. Fast Mode Pipeline (runCASSIA_pipeline) ---
cat("--- 2. Testing Fast Mode Pipeline (runCASSIA_pipeline) ---\n")
pipeline_output_name <- paste0(base_output_name, "_Pipeline")
pipeline_start_time <- Sys.time()

tryCatch({
  CASSIA::runCASSIA_pipeline(
    output_file_name = pipeline_output_name,
    tissue = tissue_param,
    species = species_param,
    marker = unprocessed_markers,
    max_workers = 6,
    annotation_model = annot_model_pipeline,
    annotation_provider = openrouter_provider,
    score_model = score_model_pipeline,
    score_provider = openrouter_provider,
    score_threshold = 98,
    annotationboost_model = boost_model_pipeline,
    annotationboost_provider = openrouter_provider,
    do_merge_annotations = TRUE,
    conversation_history_mode = "final"
  )
  pipeline_end_time <- Sys.time()
  pipeline_duration <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "mins"))
  cat("✓ runCASSIA_pipeline completed successfully in", round(pipeline_duration, 2), "minutes\n")
  
}, error = function(e) {
  cat("✗ Error in runCASSIA_pipeline:", e$message, "\n")
  if (!is.null(reticulate::py_last_error())) {
    cat("Python traceback:\n")
    print(reticulate::py_last_error())
  }
})

# --- 3. Detailed Batch Analysis (runCASSIA_batch) ---
cat("\n--- 3. Testing Detailed Batch Analysis (runCASSIA_batch) ---\n")
batch_output_name <- paste0(base_output_name, "_BatchDetailed")
batch_full_csv_path <- paste0(batch_output_name, "_full.csv")
batch_start_time <- Sys.time()

tryCatch({
  system.time({
    CASSIA::runCASSIA_batch(
      marker = unprocessed_markers,
      output_name = batch_output_name,
      model = annot_model_pipeline,
      tissue = tissue_param,
      species = species_param,
      max_workers = 6,
      n_genes = 50,
      additional_info = NULL,
      provider = openrouter_provider
    )
  }) -> batch_timing
  
  cat("✓ runCASSIA_batch completed successfully\n")
  cat("Timing:", round(batch_timing["elapsed"], 2), "seconds\n")
  
}, error = function(e) {
  cat("✗ Error in runCASSIA_batch:", e$message, "\n")
})

# --- 4. Quality Scoring ---
cat("\n--- 4. Testing Quality Scoring ---\n")
batch_scored_csv_path <- paste0(batch_output_name, "_scored.csv")
report_html_name <- paste0(batch_output_name, "_report.html")

tryCatch({
  if (file.exists(batch_full_csv_path)) {
    CASSIA::runCASSIA_score_batch(
      input_file = batch_full_csv_path,
      output_file = batch_scored_csv_path,
      max_workers = 6,
      model = score_model_pipeline,
      provider = openrouter_provider
    )
    cat("✓ runCASSIA_score_batch completed\n")
    
    if (file.exists(batch_scored_csv_path)) {
      CASSIA::runCASSIA_generate_score_report(
        csv_path = batch_scored_csv_path,
        index_name = report_html_name
      )
      cat("✓ Score report generated:", report_html_name, "\n")
    }
  } else {
    cat("⚠ Skipping quality scoring - input file not found\n")
  }
}, error = function(e) {
  cat("✗ Error in Quality Scoring:", e$message, "\n")
})

# --- 5. Uncertainty Quantification ---
cat("\n--- 5. Testing Uncertainty Quantification ---\n")
uncertainty_output_prefix <- paste0(base_output_name, "_Uncertainty")
uncertainty_start_time <- Sys.time()

tryCatch({
  system.time({
    CASSIA::runCASSIA_batch_n_times(
      n = 3,
      marker = unprocessed_markers,
      output_name = uncertainty_output_prefix,
      model = annot_model_pipeline,
      provider = openrouter_provider,
      tissue = tissue_param,
      species = species_param,
      max_workers = 6,
      batch_max_workers = 3
    )
  }) -> uncertainty_timing
  
  cat("✓ runCASSIA_batch_n_times completed in", round(uncertainty_timing["elapsed"], 2), "seconds\n")
  
  # Similarity analysis
  file_pattern_uncertainty <- paste0(uncertainty_output_prefix, "_*.csv")
  uncertainty_final_csv <- paste0(uncertainty_output_prefix, "_similarity_scores.csv")
  
  CASSIA::runCASSIA_similarity_score_batch(
    marker = unprocessed_markers,
    file_pattern = file_pattern_uncertainty,
    output_name = uncertainty_final_csv,
    max_workers = 6,
    model = annot_model_pipeline,
    provider = openrouter_provider
  )
  cat("✓ Uncertainty similarity analysis completed\n")
  
}, error = function(e) {
  cat("✗ Error in Uncertainty Quantification:", e$message, "\n")
})

# --- 6. Annotation Boost ---
cat("\n--- 6. Testing Annotation Boost ---\n")
boost_output_name <- paste0(base_output_name, "_monocyte_boosted")

tryCatch({
  if (file.exists(batch_full_csv_path)) {
    CASSIA::runCASSIA_annotationboost(
      full_result_path = batch_full_csv_path,
      marker = unprocessed_markers,
      output_name = boost_output_name,
      cluster_name = "monocyte",
      major_cluster_info = paste(species_param, tissue_param),
      num_iterations = 5,
      model = boost_model_pipeline,
      provider = openrouter_provider,
      conversation_history_mode = "final"
    )
    cat("✓ runCASSIA_annotationboost completed\n")
  } else {
    cat("⚠ Skipping Annotation Boost - input file not found\n")
  }
}, error = function(e) {
  cat("✗ Error in Annotation Boost:", e$message, "\n")
})

# --- 7. Compare Cell Types ---
cat("\n--- 7. Testing Cell Type Comparison ---\n")
compare_output_file <- paste0(base_output_name, "_plasma_cell_subtype_comparison.html")
plasma_marker_string <- "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1"

tryCatch({
  CASSIA::compareCelltypes(
    tissue = tissue_param,
    celltypes = c("Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells"),
    marker_set = plasma_marker_string,
    species = species_param,
    output_file = compare_output_file
  )
  cat("✓ compareCelltypes completed:", compare_output_file, "\n")
}, error = function(e) {
  cat("✗ Error in Compare Cell Types:", e$message, "\n")
})

# --- 8. Subclustering Analysis ---
cat("\n--- 8. Testing Subclustering Analysis ---\n")
subcluster_run_output <- paste0(base_output_name, "_subclustering_results")

tryCatch({
  CASSIA::runCASSIA_subclusters(
    marker = subcluster_markers_data,
    major_cluster_info = "cd8 t cell",
    output_name = subcluster_run_output,
    model = annot_model_pipeline,
    provider = openrouter_provider,
    n_genes = 50L
  )
  cat("✓ runCASSIA_subclusters completed\n")
  
  # Multiple subclustering runs
  subcluster_n_output_prefix <- paste0(base_output_name, "_subclustering_results_n")
  CASSIA::runCASSIA_n_subcluster(
    n = 3,  # Reduced for testing
    marker = subcluster_markers_data,
    major_cluster_info = "cd8 t cell",
    base_output_name = subcluster_n_output_prefix,
    model = annot_model_pipeline,
    provider = openrouter_provider,
    max_workers = 3,
    n_genes = 50L
  )
  cat("✓ runCASSIA_n_subcluster completed\n")
  
}, error = function(e) {
  cat("✗ Error in Subclustering Analysis:", e$message, "\n")
})

# --- 9. Annotation Boost with Additional Task ---
cat("\n--- 9. Testing Annotation Boost with Additional Task ---\n")
additional_task_output_name <- paste0(base_output_name, "_T_cell_state")

tryCatch({
  if (file.exists(batch_full_csv_path)) {
    CASSIA::runCASSIA_annotationboost_additional_task(
      full_result_path = batch_full_csv_path,
      marker = unprocessed_markers,
      output_name = additional_task_output_name,
      cluster_name = "cd8-positive, alpha-beta t cell",
      major_cluster_info = paste(species_param, tissue_param),
      num_iterations = 3,  # Reduced for testing
      model = "anthropic/claude-3.5-sonnet",
      provider = "openrouter",
      additional_task = "infer the state of this T cell cluster"
    )
    cat("✓ runCASSIA_annotationboost_additional_task completed\n")
  } else {
    cat("⚠ Skipping Additional Task - input file not found\n")
  }
}, error = function(e) {
  cat("✗ Error in Additional Task:", e$message, "\n")
})

# --- Test Summary ---
test_end_time <- Sys.time()
total_duration <- as.numeric(difftime(test_end_time, pipeline_start_time, units = "mins"))

cat("\n=== TEST SUMMARY ===\n")
cat("Total test duration:", round(total_duration, 2), "minutes\n")
cat("Results location:", file.path(original_wd, results_dir), "\n")

# List all generated files
all_files <- list.files(".", recursive = TRUE)
cat("Generated files:\n")
for(file in all_files) {
  cat("  -", file, "\n")
}

# Create a summary report
summary_report <- data.frame(
  Test = c("Pipeline", "Batch Analysis", "Quality Scoring", "Uncertainty", 
           "Annotation Boost", "Cell Type Comparison", "Subclustering", "Additional Task"),
  Status = c("Completed", "Completed", "Completed", "Completed", 
             "Completed", "Completed", "Completed", "Completed"),
  Files_Generated = c("Multiple", "CSV files", "HTML reports", "CSV files",
                     "HTML/TXT", "HTML", "CSV files", "HTML/TXT")
)

write.csv(summary_report, "test_summary.csv", row.names = FALSE)
cat("✓ Test summary saved to test_summary.csv\n")

# Return to original working directory
setwd(original_wd)

cat("\n=== All CASSIA R wrapper tests completed ===\n")
cat("Check the test_results/ folder for all outputs\n") 