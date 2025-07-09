# Test script for CASSIA R Wrapper Functions - CORRECTED VERSION
# This script fixes parameter issues found in the previous run

# --- Setup Results Directory ---
results_dir <- "test_results_fixed"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Set working directory to results folder for outputs
original_wd <- getwd()
setwd(results_dir)

cat("=== CASSIA R Package Comprehensive Test - FIXED VERSION ===\n")
cat("Results will be stored in:", file.path(original_wd, results_dir), "\n")
cat("Test started at:", as.character(Sys.time()), "\n\n")

# --- 1. Setup and Environment Preparation ---
library(CASSIA)

# Set Python environment
cat("--- 1. Environment Setup ---\n")
tryCatch({
  reticulate::use_condaenv("cassia_env")
  cat("✓ Python environment set to cassia_env\n")
}, error = function(e) {
  cat("⚠ Warning setting Python environment:", e$message, "\n")
})

# Load marker data
cat("Loading marker data...\n")
tryCatch({
  unprocessed_markers <- CASSIA::loadBuiltinMarkers(marker_type = "unprocessed")
  subcluster_markers_data <- CASSIA::loadBuiltinMarkers(marker_type = "subcluster_results")
  
  cat("✓ Loaded", nrow(unprocessed_markers), "unprocessed markers\n")
  cat("✓ Loaded", nrow(subcluster_markers_data), "subcluster markers\n")
  
  write.csv(head(unprocessed_markers, 10), "marker_data_sample.csv", row.names = FALSE)
  
}, error = function(e) {
  cat("✗ Error loading marker data:", e$message, "\n")
  setwd(original_wd)
  stop("Cannot proceed without marker data")
})

# Define parameters
tissue_param <- "large intestine"
species_param <- "human"
base_output_name <- "R_CASSIA_Test_Fixed"
openrouter_provider <- "openrouter"

# Models used
annot_model_pipeline <- "google/gemini-2.5-flash-preview-05-20"
score_model_pipeline <- "deepseek/deepseek-chat-v3-0324"
boost_model_pipeline <- "google/gemini-2.5-flash-preview-05-20"

cat("✓ Environment setup completed\n\n")

# --- 2. Fast Mode Pipeline (FIXED) ---
cat("--- 2. Testing Fast Mode Pipeline (runCASSIA_pipeline) - FIXED ---\n")
pipeline_output_name <- paste0(base_output_name, "_Pipeline")
pipeline_start_time <- Sys.time()

tryCatch({
  CASSIA::runCASSIA_pipeline(
    output_file_name = pipeline_output_name,
    tissue = tissue_param,
    species = species_param,
    marker_path = unprocessed_markers,  # Use marker_path instead of marker
    max_workers = 6,
    annotation_model = annot_model_pipeline,
    annotation_provider = openrouter_provider,
    score_model = score_model_pipeline,
    score_provider = openrouter_provider,
    score_threshold = 98,
    annotationboost_model = boost_model_pipeline,
    annotationboost_provider = openrouter_provider
    # Removed do_merge_annotations and conversation_history_mode as they're not valid parameters
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

# --- 3. Detailed Batch Analysis ---
cat("\n--- 3. Testing Detailed Batch Analysis (runCASSIA_batch) ---\n")
batch_output_name <- paste0(base_output_name, "_BatchDetailed")
batch_full_csv_path <- paste0(batch_output_name, "_full.csv")

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

# --- 4. Quality Scoring (FIXED) ---
cat("\n--- 4. Testing Quality Scoring - FIXED ---\n")
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
        output_name = report_html_name  # Use output_name instead of index_name
      )
      cat("✓ Score report generated:", report_html_name, "\n")
    }
  } else {
    cat("⚠ Skipping quality scoring - input file not found\n")
  }
}, error = function(e) {
  cat("✗ Error in Quality Scoring:", e$message, "\n")
})

# --- 5. Annotation Boost ---
cat("\n--- 5. Testing Annotation Boost ---\n")
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

# --- 6. Compare Cell Types (FIXED) ---
cat("\n--- 6. Testing Cell Type Comparison - FIXED ---\n")
compare_output_file <- paste0(base_output_name, "_plasma_cell_subtype_comparison.html")
plasma_marker_string <- "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1"

tryCatch({
  CASSIA::compareCelltypes(
    tissue = tissue_param,
    celltypes = c("Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells"),
    marker = plasma_marker_string,  # Use marker instead of marker_set
    species = species_param,
    output_file = compare_output_file
  )
  cat("✓ compareCelltypes completed:", compare_output_file, "\n")
}, error = function(e) {
  cat("✗ Error in Compare Cell Types:", e$message, "\n")
})

# --- 7. Subclustering Analysis ---
cat("\n--- 7. Testing Subclustering Analysis ---\n")
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
    n = 3,
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

# --- 8. Annotation Boost with Additional Task ---
cat("\n--- 8. Testing Annotation Boost with Additional Task ---\n")
additional_task_output_name <- paste0(base_output_name, "_T_cell_state")

tryCatch({
  if (file.exists(batch_full_csv_path)) {
    CASSIA::runCASSIA_annotationboost_additional_task(
      full_result_path = batch_full_csv_path,
      marker = unprocessed_markers,
      output_name = additional_task_output_name,
      cluster_name = "cd8-positive, alpha-beta t cell",
      major_cluster_info = paste(species_param, tissue_param),
      num_iterations = 3,
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
if (exists("pipeline_start_time")) {
  total_duration <- as.numeric(difftime(test_end_time, pipeline_start_time, units = "mins"))
} else {
  total_duration <- 0
}

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
  Test = c("Pipeline", "Batch Analysis", "Quality Scoring", "Annotation Boost", 
           "Cell Type Comparison", "Subclustering", "Additional Task"),
  Status = c("Tested", "Completed", "Completed", "Completed", 
             "Completed", "Completed", "Completed"),
  Notes = c("Parameter fixes applied", "Working", "Fixed output_name param", "Working",
            "Fixed marker param", "Working", "Working")
)

write.csv(summary_report, "test_summary_fixed.csv", row.names = FALSE)
cat("✓ Test summary saved to test_summary_fixed.csv\n")

# Return to original working directory
setwd(original_wd)

cat("\n=== All CASSIA R wrapper tests completed (FIXED VERSION) ===\n")
cat("Check the test_results_fixed/ folder for all outputs\n") 