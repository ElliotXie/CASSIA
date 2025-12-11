# CASSIA Test 17: Provider Defaults (R) - INSTALL MODE
# =====================================================
# Tests the overall_provider parameter via R package.
# Uses devtools::install_local() for full package installation testing.
#
# Usage:
#     Rscript test_provider_defaults_install.R
#
# Functions tested:
# - runCASSIA_pipeline() with overall_provider parameter

# Get script directory (works with Rscript and source())
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg), winslash = "/")))
  }
  for (i in sys.nframe():1) {
    if (!is.null(sys.frame(i)$ofile)) {
      return(dirname(normalizePath(sys.frame(i)$ofile, winslash = "/")))
    }
  }
  return(normalizePath(getwd(), winslash = "/"))
}
script_dir <- get_script_dir()

# Source shared utilities
source(file.path(script_dir, "..", "shared", "r", "test_utils.R"))
source(file.path(script_dir, "..", "shared", "r", "fixtures.R"))
source(file.path(script_dir, "..", "shared", "r", "result_manager.R"))
source(file.path(script_dir, "..", "shared", "r", "logging_manager.R"))

run_provider_defaults_install_test <- function() {
  print_test_header("17 - Provider Defaults (R) [INSTALL MODE]")

  # Load configuration
  config <- load_config()
  print_config_summary(config)

  # Setup API keys
  setup_api_keys()

  # Setup CASSIA R package - INSTALL MODE (full package installation)
  tryCatch({
    setup_cassia_install(force = TRUE)
    message("CASSIA R package installed and loaded successfully")
  }, error = function(e) {
    message("Error installing CASSIA: ", e$message)
    return(FALSE)
  })

  # Get settings
  llm_config <- config$llm
  data_config <- config$data

  # Create results directory
  results_dirs <- create_results_dir("17_provider_defaults", get_test_mode())

  start_logging(results_dirs$logs)
  log_msg("Results will be saved to:", results_dirs$base)

  # Use only 2 clusters for faster testing
  test_clusters <- c("monocyte", "plasma cell")

  # Load raw marker data
  data_folder <- file.path(script_dir, "..", "16_cassia_pipeline", "data")
  raw_markers_path <- file.path(data_folder, "scanpy_rank_genes_df.csv")

  if (!file.exists(raw_markers_path)) {
    stop(paste("Raw marker file not found:", raw_markers_path))
  }

  # Load and filter markers
  raw_markers <- read.csv(raw_markers_path)
  marker_df <- raw_markers[raw_markers$group %in% test_clusters, ]

  log_msg("\nTesting pipeline with overall_provider parameter")
  log_msg("Testing for", length(test_clusters), "clusters:")
  for (cluster in test_clusters) {
    log_msg("  -", cluster)
  }

  # Test overall_provider parameter - using openrouter as specified
  test_provider <- "openrouter"
  output_name <- file.path(results_dirs$outputs, "provider_defaults_test")

  log_msg("\n--- Testing overall_provider='", test_provider, "' ---")
  log_msg("Using provider defaults automatically")

  # Change to results directory
  original_dir <- getwd()
  setwd(results_dirs$outputs)

  # Run tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"

  tryCatch({
    log_msg("\nRunning runCASSIA_pipeline with overall_provider...")

    CASSIA::runCASSIA_pipeline(
      output_file_name = output_name,
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      marker = marker_df,
      max_workers = as.integer(llm_config$max_workers %||% 3),
      overall_provider = test_provider,
      # All model params are NULL - should use provider defaults
      annotation_model = NULL,
      annotation_provider = NULL,
      score_model = NULL,
      score_provider = NULL,
      annotationboost_model = NULL,
      annotationboost_provider = NULL,
      do_merge_annotations = TRUE,
      merge_model = NULL,
      merge_provider = NULL,
      score_threshold = 99,
      validator_involvement = config$validator$default %||% "v1"
    )

    # Find the pipeline output directory
    output_items <- list.files(results_dirs$outputs, pattern = "^CASSIA_", full.names = TRUE)
    output_dirs <- output_items[dir.exists(output_items)]

    if (length(output_dirs) > 0) {
      pipeline_output_dir <- output_dirs[1]
      log_msg("\nPipeline output directory:", basename(pipeline_output_dir))

      csv_dir <- file.path(pipeline_output_dir, "03_csv_files")

      if (dir.exists(csv_dir)) {
        log_msg("  [OK] 03_csv_files exists")
        final_results_files <- list.files(csv_dir, pattern = "_FINAL_RESULTS\\.csv$", full.names = TRUE)
        if (length(final_results_files) > 0) {
          log_msg("  [OK] FINAL_RESULTS.csv exists:", basename(final_results_files[1]))
          results_df <- read.csv(final_results_files[1])
          log_msg("       - Contains", nrow(results_df), "rows")
          status <- "passed"
        } else {
          log_msg("  [WARN] FINAL_RESULTS.csv not found")
          status <- "failed"
          errors <- list("FINAL_RESULTS.csv not found")
        }
      } else {
        log_msg("  [FAIL] 03_csv_files missing")
        status <- "failed"
        errors <- list("03_csv_files directory missing")
      }
    } else {
      status <- "failed"
      errors <- list("Pipeline output directory not created")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    log_msg("\nError:", e$message)
  })

  # Change back to original directory
  setwd(original_dir)

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "provider_defaults_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(test_clusters),
    errors = errors
  )
  metadata$test_provider <- test_provider
  save_test_metadata(results_dirs$outputs, metadata)

  save_test_results(results_dirs$outputs, list(
    overall_provider = test_provider,
    clusters_tested = test_clusters,
    mode = "install"
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  stop_logging()
  return(success)
}

# Run test
success <- run_provider_defaults_install_test()
quit(status = if (success) 0 else 1)
