# CASSIA Test 01: Single Cluster Annotation (R) - INSTALL MODE
# =============================================================
# Tests the runCASSIA function on a single cell cluster via R package.
# Uses devtools::install_local() for full package installation testing.
#
# Usage:
#     Rscript test_single_annotation_install.R

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

run_single_annotation_test <- function() {
  print_test_header("01 - Single Cluster Annotation (R) [INSTALL MODE]")

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

  # Test cluster
  test_cluster <- "plasma cell"
  n_genes <- data_config$n_genes %||% 30

  cat("\nTesting single annotation for:", test_cluster, "\n")
  cat("Using top", n_genes, "markers\n")

  # Get markers for the test cluster
  marker_df <- get_marker_dataframe_for_cluster(test_cluster, n_genes)
  cat("Loaded", nrow(marker_df), "markers\n")

  # Create results directory
  results_dir <- create_results_dir("01_single_annotation")
  cat("Results will be saved to:", results_dir, "\n")

  # Run the test
  start_time <- Sys.time()
  errors <- list()
  result <- NULL
  status <- "error"

  tryCatch({
    cat("\nRunning runCASSIA via R package (install mode)...\n")

    # Call CASSIA R function
    result <- CASSIA::runCASSIA(
      model = llm_config$model %||% "google/gemini-2.5-flash",
      temperature = llm_config$temperature %||% 0.3,
      marker_list = marker_df,
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      provider = llm_config$provider %||% "openrouter",
      validator_involvement = config$validator$default %||% "v1"
    )

    # Validate result
    validation <- validate_result(result)

    if (validation$valid) {
      cat("\nAnnotation Result:\n")
      cat("  Main cell type:", result$main_cell_type, "\n")
      cat("  Sub cell types:", paste(result$sub_cell_types, collapse = ", "), "\n")
      status <- "passed"
    } else {
      errors <- validation$errors
      status <- "failed"
      cat("\nValidation errors:", paste(errors, collapse = ", "), "\n")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    cat("\nError:", e$message, "\n")
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save results
  metadata <- create_test_metadata(
    test_name = "single_annotation_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(test_cluster),
    errors = errors
  )
  save_test_metadata(results_dir, metadata)

  if (!is.null(result)) {
    save_test_results(results_dir, list(
      cluster = test_cluster,
      result = result,
      mode = "install"
    ))
  }

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  return(success)
}

# Run test
success <- run_single_annotation_test()
quit(status = if (success) 0 else 1)
