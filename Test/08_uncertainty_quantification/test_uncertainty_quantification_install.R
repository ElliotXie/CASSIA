# CASSIA Test 08: Uncertainty Quantification (R) - INSTALL MODE
# ==============================================================
# Tests the uncertainty quantification functions via R package.
# Uses devtools::install_local() for full package installation testing.
#
# Usage:
#     Rscript test_uncertainty_quantification_install.R
#
# Functions tested:
# - runCASSIA_n_times_similarity_score(): Run n single analyses with similarity score

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

run_uncertainty_quantification_install_test <- function() {
  print_test_header("08 - Uncertainty Quantification (R) [INSTALL MODE]")

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
  results_dirs <- create_results_dir("08_uncertainty_quantification", get_test_mode())

  start_logging(results_dirs$logs)
  log_msg("Results will be saved to:", results_dirs$base)

  # Test cluster
  test_cluster <- "plasma cell"
  markers <- get_cluster_markers(test_cluster)

  # Run tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  uq_results <- list()

  tryCatch({
    log_msg("\n--- Test: runCASSIA_n_times_similarity_score ---")
    log_msg("  Cluster:", test_cluster)
    log_msg("  Model:", llm_config$model %||% "google/gemini-2.5-flash")
    log_msg("  Provider:", llm_config$provider %||% "openrouter")
    log_msg("  N iterations: 3 (reduced for testing)")

    result <- CASSIA::runCASSIA_n_times_similarity_score(
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      additional_info = NULL,
      temperature = llm_config$temperature %||% 0.3,
      marker_list = markers,
      model = llm_config$model %||% "google/gemini-2.5-flash",
      max_workers = llm_config$max_workers %||% 3,
      n = 3,  # Reduced for faster testing
      provider = llm_config$provider %||% "openrouter",
      validator_involvement = config$validator$default %||% "v1"
    )

    # Check result
    if (is.list(result)) {
      log_msg("\nResults:")
      log_msg("  General cell type (LLM):", result$general_celltype_llm %||% "N/A")
      log_msg("  Sub cell type (LLM):", result$sub_celltype_llm %||% "N/A")
      log_msg("  Similarity score:", result$similarity_score %||% "N/A")

      uq_results <- list(
        general_celltype_llm = result$general_celltype_llm,
        sub_celltype_llm = result$sub_celltype_llm,
        similarity_score = result$similarity_score,
        consensus_types = result$consensus_types,
        mixed_celltypes = result$Possible_mixed_celltypes_llm
      )

      # Validate results
      if (!is.null(result$general_celltype_llm) && !is.null(result$similarity_score)) {
        status <- "passed"
      } else {
        status <- "failed"
        errors <- list("Missing expected result fields")
      }
    } else {
      status <- "failed"
      errors <- list("Unexpected result format")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    log_msg("\nError:", e$message)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "uncertainty_quantification_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(test_cluster),
    errors = errors
  )
  save_test_metadata(results_dirs$outputs, metadata)

  save_test_results(results_dirs$outputs, list(
    test_cluster = test_cluster,
    n_iterations = 3,
    results = uq_results,
    mode = "install"
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_uncertainty_quantification_install_test()
quit(status = if (success) 0 else 1)
