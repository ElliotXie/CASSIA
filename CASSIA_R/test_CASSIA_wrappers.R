# Test script for CASSIA R Wrapper Functions
# This script mirrors the workflow in the CASSIA_python_package_test.ipynb notebook.

# --- 1. Setup and Environment Preparation ---

# Load the CASSIA library
library(CASSIA)

# API keys are already loaded in system environment variables - no need to set them here

# It's good practice to ensure the Python environment is set up.
# If this is the first time or you have issues, uncomment and run:
# CASSIA::setup_cassia_env()
# CASSIA::set_python_env("cassia_env") # or your chosen environment name

# Check if the environment is okay (optional)
# print(CASSIA::check_python_env())

# Load example marker data
# The Python notebook uses CASSIA.loadmarker().
# The R package has loadBuiltinMarkers() or loadExampleMarkers().
# We'll use loadBuiltinMarkers as it's more analogous to the Python notebook's direct CASSIA module calls.
message("Loading marker data...")
tryCatch({
  unprocessed_markers <- CASSIA::loadBuiltinMarkers(marker_type = "unprocessed")
  subcluster_markers_data <- CASSIA::loadBuiltinMarkers(marker_type = "subcluster_results")
  
  print("Available marker sets in the package:")
  print(CASSIA::listAvailableMarkers())
  
  message("Unprocessed markers head:")
  print(head(unprocessed_markers))
  message("Subcluster markers head:")
  print(head(subcluster_markers_data))
}, error = function(e) {
  message("Error loading marker data: ", e$message)
})

# Define common parameters based on the Python notebook
tissue_param <- "large intestine"
species_param <- "human"
base_output_name <- "R_CASSIA_Test"
openrouter_provider <- "openrouter" # Default provider in many notebook calls

# Models used in the Python notebook (ensure these are available via your OpenRouter key)
annot_model_pipeline <- "google/gemini-2.5-flash-preview-05-20"
score_model_pipeline <- "deepseek/deepseek-chat-v3-0324"
boost_model_pipeline <- "google/gemini-2.5-flash-preview-05-20"

# --- 2. Fast Mode Pipeline (runCASSIA_pipeline) ---
message("\n--- Starting Test: 2. Fast Mode Pipeline (runCASSIA_pipeline) ---")
pipeline_output_name <- paste0(base_output_name, "_Pipeline")
tryCatch({
  CASSIA::runCASSIA_pipeline(
    output_file_name = pipeline_output_name,
    tissue = tissue_param,
    species = species_param,
    marker = unprocessed_markers, # Using the loaded R data frame
    max_workers = 6,
    annotation_model = annot_model_pipeline,
    annotation_provider = openrouter_provider,
    score_model = score_model_pipeline,
    score_provider = openrouter_provider,
    score_threshold = 98,
    annotationboost_model = boost_model_pipeline,
    annotationboost_provider = openrouter_provider,
    # ranking_method = "p_val_adj" # This parameter is in Python but not in R wrapper for pipeline
    do_merge_annotations = TRUE, # Explicitly set, default is TRUE
    conversation_history_mode = "final" # Explicitly set, default is "final"
  )
  message("runCASSIA_pipeline test completed successfully!")
}, error = function(e) {
  message("Error in runCASSIA_pipeline: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# --- 3. Detailed Batch Analysis (runCASSIA_batch) ---
message("\n--- Starting Test: 3. Detailed Batch Analysis (runCASSIA_batch) ---")
batch_output_name <- paste0(base_output_name, "_BatchDetailed")
batch_full_csv_path <- paste0(batch_output_name, "_full.csv") # Path for the next step
tryCatch({
  CASSIA::runCASSIA_batch(
    marker = unprocessed_markers,
    output_name = batch_output_name,
    model = annot_model_pipeline, # Using a model from the pipeline for consistency
    tissue = tissue_param,
    species = species_param,
    max_workers = 6,
    n_genes = 50,
    additional_info = NULL,
    provider = openrouter_provider
  )
  message("runCASSIA_batch test completed. Check files starting with: ", batch_output_name)
}, error = function(e) {
  message("Error in runCASSIA_batch: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# --- 4. Quality Scoring (runCASSIA_score_batch & runCASSIA_generate_score_report) ---
message("\n--- Starting Test: 4. Quality Scoring ---")
batch_scored_csv_path <- paste0(batch_output_name, "_scored.csv")
report_html_name <- paste0(batch_output_name, "_report.html")
tryCatch({
  if (file.exists(batch_full_csv_path)) {
    CASSIA::runCASSIA_score_batch(
      input_file = batch_full_csv_path,
      output_file = batch_scored_csv_path,
      max_workers = 6,
      model = score_model_pipeline, # Using score model from pipeline
      provider = openrouter_provider
    )
    message("runCASSIA_score_batch test completed.")
    
    if (file.exists(batch_scored_csv_path)) {
      CASSIA::runCASSIA_generate_score_report(
        csv_path = batch_scored_csv_path,
        output_name = report_html_name # R function uses output_name, Python uses index_name
      )
      message("runCASSIA_generate_score_report test completed. Check file: ", report_html_name)
    } else {
      message("Skipping report generation, scored CSV not found: ", batch_scored_csv_path)
    }
  } else {
    message("Skipping quality scoring, input CSV not found: ", batch_full_csv_path)
  }
}, error = function(e) {
  message("Error in Quality Scoring: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# --- 5. Uncertainty Quantification (runCASSIA_batch_n_times & runCASSIA_similarity_score_batch) ---
message("\n--- Starting Test: 5. Uncertainty Quantification ---")
uncertainty_output_prefix <- paste0(base_output_name, "_Uncertainty")
uncertainty_final_csv <- paste0(uncertainty_output_prefix, "_similarity_scores.csv") # Example, actual output name might differ based on Python

tryCatch({
  CASSIA::runCASSIA_batch_n_times(
    n = 3,
    marker = unprocessed_markers,
    output_name = uncertainty_output_prefix, # R function uses output_name (prefix)
    model = annot_model_pipeline,
    provider = openrouter_provider,
    tissue = tissue_param,
    species = species_param,
    max_workers = 6,
    batch_max_workers = 3
  )
  message("runCASSIA_batch_n_times test completed.")
  
  # The file_pattern needs to match the output of runCASSIA_batch_n_times
  # Python notebook: output_name + "_Uncertainty_*_full.csv"
  # R function creates files like: uncertainty_output_prefix_0_full.csv, uncertainty_output_prefix_1_full.csv etc.
  file_pattern_uncertainty <- paste0(uncertainty_output_prefix, "_*_full.csv")
  
  CASSIA::runCASSIA_similarity_score_batch(
    marker = unprocessed_markers,
    file_pattern = file_pattern_uncertainty,
    output_name = uncertainty_final_csv, # R output_name is the direct file name
    max_workers = 6,
    model = annot_model_pipeline, # Python notebook uses a specific model here
    provider = openrouter_provider
  )
  message("runCASSIA_similarity_score_batch for uncertainty test completed. Check file: ", uncertainty_final_csv)
  
}, error = function(e) {
  message("Error in Uncertainty Quantification: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# --- 6. Annotation Boost (runCASSIA_annotationboost) ---
message("\n--- Starting Test: 6. Annotation Boost ---")
# Needs an input CSV from a batch run, e.g., batch_full_csv_path
# The Python notebook boosts "monocyte"
boost_output_name <- paste0(base_output_name, "_monocyte_boosted")
tryCatch({
  if (file.exists(batch_full_csv_path)) {
    CASSIA::runCASSIA_annotationboost(
      full_result_path = batch_full_csv_path,
      marker = unprocessed_markers,
      output_name = boost_output_name,
      cluster_name = "monocyte", # Example from Python notebook
      major_cluster_info = paste(species_param, tissue_param), # e.g. "Human Large Intestine"
      num_iterations = 5,
      model = boost_model_pipeline,
      provider = openrouter_provider,
      conversation_history_mode = "final"
    )
    message("runCASSIA_annotationboost test completed. Check files starting with: ", boost_output_name)
  } else {
    message("Skipping Annotation Boost, input CSV not found: ", batch_full_csv_path)
  }
}, error = function(e) {
  message("Error in Annotation Boost: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# --- 7. Compare Cell Types (compareCelltypes) ---
message("\n--- Starting Test: 7. Compare Cell Types ---")
compare_output_file <- paste0(base_output_name, "_plasma_cell_subtype_comparison.html")
# Marker string from Python notebook for plasma cell comparison
plasma_marker_string <- "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

tryCatch({
  CASSIA::compareCelltypes(
    tissue = tissue_param,
    celltypes = c("Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"),
    marker = plasma_marker_string, # R function takes 'marker', Python takes 'marker_set'
    species = species_param,
    output_file = compare_output_file,
    model_list = NULL # Using default models in the Python function
  )
  message("compareCelltypes test completed. Check file: ", compare_output_file)
}, error = function(e) {
  message("Error in Compare Cell Types: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# --- 8. Subclustering Analysis ---
message("\n--- Starting Test: 8. Subclustering Analysis ---")
# runCASSIA_subclusters
subcluster_run_output <- paste0(base_output_name, "_subclustering_results")
subcluster_run_csv <- paste0(subcluster_run_output, ".csv") # Expected output
tryCatch({
  CASSIA::runCASSIA_subclusters(
    marker = subcluster_markers_data, # Using the loaded R data frame for subclusters
    major_cluster_info = "cd8 t cell", # Example from Python notebook
    output_name = subcluster_run_output,
    model = annot_model_pipeline,
    provider = openrouter_provider,
    n_genes = 50L
  )
  message("runCASSIA_subclusters test completed. Check file: ", subcluster_run_csv)
}, error = function(e) {
  message("Error in runCASSIA_subclusters: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# runCASSIA_n_subcluster
subcluster_n_output_prefix <- paste0(base_output_name, "_subclustering_results_n")
tryCatch({
  CASSIA::runCASSIA_n_subcluster(
    n = 5, # Python notebook uses n=5
    marker = subcluster_markers_data,
    major_cluster_info = "cd8 t cell",
    base_output_name = subcluster_n_output_prefix, # R uses base_output_name
    model = annot_model_pipeline,
    provider = openrouter_provider,
    max_workers = 5,
    n_genes = 50L
  )
  message("runCASSIA_n_subcluster test completed. Check files starting with: ", subcluster_n_output_prefix)
}, error = function(e) {
  message("Error in runCASSIA_n_subcluster: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# runCASSIA_similarity_score_batch for subclusters
subcluster_similarity_output <- paste0(base_output_name, "_subclustering_uncertainty.csv")
# Python notebook: "subclustering_results_n_*.csv"
file_pattern_subcluster_n <- paste0(subcluster_n_output_prefix, "_*_full.csv") # Match expected output of R's n_subcluster

tryCatch({
  CASSIA::runCASSIA_similarity_score_batch(
    marker = subcluster_markers_data,
    file_pattern = file_pattern_subcluster_n,
    output_name = subcluster_similarity_output,
    max_workers = 6,
    model = annot_model_pipeline,
    provider = openrouter_provider
  )
  message("runCASSIA_similarity_score_batch for subclusters test completed. Check file: ", subcluster_similarity_output)
}, error = function(e) {
  message("Error in runCASSIA_similarity_score_batch for subclusters: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

# --- 9. Annotation Boost with Additional Task ---
message("\n--- Starting Test: 9. Annotation Boost with Additional Task ---")
# Needs an input CSV, e.g., batch_full_csv_path
additional_task_output_name <- paste0(base_output_name, "_T_cell_state")
tryCatch({
  if (file.exists(batch_full_csv_path)) {
    CASSIA::runCASSIA_annotationboost_additional_task(
      full_result_path = batch_full_csv_path,
      marker = unprocessed_markers,
      output_name = additional_task_output_name,
      cluster_name = "cd8-positive, alpha-beta t cell", # Example from Python notebook
      major_cluster_info = paste(species_param, tissue_param),
      num_iterations = 5, # Python notebook default is 5
      model = "anthropic/claude-3.5-sonnet", # Specific model from notebook
      provider = "openrouter", # Python notebook does not specify provider here, assumes it's handled
      additional_task = "infer the state of this T cell cluster"
    )
    message("runCASSIA_annotationboost_additional_task test completed. Check files: ", additional_task_output_name)
  } else {
    message("Skipping Annotation Boost with Additional Task, input CSV not found: ", batch_full_csv_path)
  }
}, error = function(e) {
  message("Error in Annotation Boost with Additional Task: ", e$message)
  if (!is.null(reticulate::py_last_error())) {
    message("Python traceback: ")
    print(reticulate::py_last_error())
  }
})

message("\n--- All CASSIA R wrapper tests concluded. ---")
message("Please check the console output for messages and the working directory for generated files.")
# The placeholder 'δείτε το timestamp' in the pipeline message should be replaced with actual timestamp logic if needed,
# but typically the Python functions handle the exact folder naming. This R script focuses on calling the wrappers. 