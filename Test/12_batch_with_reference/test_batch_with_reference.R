# CASSIA Test 12: Batch Annotation with Reference (R)
# =====================================================
# Tests the runCASSIA_batch_with_reference function via R/reticulate.
# Verifies the two-step ReAct reference workflow.
#
# Usage:
#     Rscript test_batch_with_reference.R

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

run_batch_with_reference_test <- function() {
  print_test_header("12 - Batch Annotation with Reference (R)")

  # Load configuration
  config <- load_config()
  print_config_summary(config)

  # Setup API keys
  setup_api_keys()

  # Setup CASSIA R package (for Python interop via reticulate)
  tryCatch({
    setup_cassia()
    message("CASSIA R package loaded successfully")
  }, error = function(e) {
    message("Error loading CASSIA: ", e$message)
    return(FALSE)
  })

  # Get settings
  llm_config <- config$llm
  data_config <- config$data

  # Create results directory
  results <- create_results_dir("12_batch_with_reference", get_test_mode())

  start_logging(results$logs)
  log_msg("Results will be saved to:", results$base)

  # Get marker data
  marker_df <- get_full_marker_dataframe()
  all_clusters <- get_all_clusters()
  log_msg("\nMarker data loaded:", nrow(marker_df), "clusters")

  # Run tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  test_results <- list()

  # Import Python function via reticulate
  tryCatch({
    # Import the function directly from CASSIA
    cassia_module <- reticulate::import("CASSIA")
    run_batch_with_ref <- cassia_module$runCASSIA_batch_with_reference

    log_msg("\n--- Running runCASSIA_batch_with_reference with use_reference=TRUE ---")
    log_msg("  Model:", llm_config$model %||% "google/gemini-2.5-flash")
    log_msg("  Provider:", llm_config$provider %||% "openrouter")
    log_msg("  Clusters:", length(all_clusters))

    output_name <- file.path(results$outputs, "batch_ref_results")

    # Run batch with reference
    run_batch_with_ref(
      marker = marker_df,
      output_name = output_name,
      n_genes = as.integer(data_config$n_genes %||% 30),
      model = llm_config$model %||% "google/gemini-2.5-flash",
      temperature = llm_config$temperature %||% 0.3,
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      max_workers = as.integer(llm_config$max_workers %||% 3),
      provider = llm_config$provider %||% "openrouter",
      validator_involvement = config$validator$default %||% "v1",
      use_reference = TRUE,
      verbose = TRUE
    )

    # Check output files
    summary_csv <- paste0(output_name, "_summary.csv")

    if (file.exists(summary_csv)) {
      results_df <- read.csv(summary_csv, stringsAsFactors = FALSE)
      clusters_annotated <- nrow(results_df)

      # Verify reference columns exist
      has_ref_column <- "Reference.Used" %in% names(results_df)

      log_msg("\nBatch with Reference Results:")
      log_msg("  Clusters annotated:", clusters_annotated, "/", length(all_clusters))
      log_msg("  Reference column present:", has_ref_column)

      if (has_ref_column) {
        ref_usage <- table(results_df$Reference.Used)
        log_msg("  Reference usage:")
        print(ref_usage)
      }

      test_results <- list(
        clusters_annotated = clusters_annotated,
        total_clusters = length(all_clusters),
        has_reference_column = has_ref_column,
        output_file = summary_csv
      )

      if (clusters_annotated == length(all_clusters) && has_ref_column) {
        status <- "passed"
      } else {
        status <- "failed"
        if (!has_ref_column) {
          errors <- list("Reference Used column not found")
        }
      }
    } else {
      status <- "failed"
      errors <- list("Output file not created")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    log_msg("\nError:", e$message)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "batch_with_reference",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(all_clusters),
    errors = errors
  )
  metadata$use_reference <- TRUE
  save_test_metadata(results$outputs, metadata)

  save_test_results(results$outputs, list(
    test_results = test_results
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_batch_with_reference_test()
quit(status = if (success) 0 else 1)
