# CASSIA Test 02: Batch Annotation (R)
# =====================================
# Tests the runCASSIA_batch function on all 6 cell clusters via R package.
#
# Usage:
#     Rscript test_batch_annotation.R

# Get script directory (works with Rscript and source())
get_script_dir <- function() {
  # Method 1: commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg), winslash = "/")))
  }
  # Method 2: sys.frame (works with source())
  for (i in sys.nframe():1) {
    if (!is.null(sys.frame(i)$ofile)) {
      return(dirname(normalizePath(sys.frame(i)$ofile, winslash = "/")))
    }
  }
  # Fallback: current working directory
  return(normalizePath(getwd(), winslash = "/"))
}
script_dir <- get_script_dir()

# Source shared utilities
source(file.path(script_dir, "..", "shared", "r", "test_utils.R"))
source(file.path(script_dir, "..", "shared", "r", "fixtures.R"))
source(file.path(script_dir, "..", "shared", "r", "result_manager.R"))
source(file.path(script_dir, "..", "shared", "r", "logging_manager.R"))

run_batch_annotation_test <- function() {
  print_test_header("02 - Batch Annotation (R)")

  # Load configuration
  config <- load_config()
  print_config_summary(config)

  # Setup API keys
  setup_api_keys()

  # Setup CASSIA R package
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

  # Get all clusters
  all_clusters <- get_all_clusters()

  # Create results directory (BEFORE any logging)
  results <- create_results_dir("02_batch_annotation", get_test_mode())
  start_logging(results$logs)

  log_msg("Testing batch annotation for", length(all_clusters), "clusters:")
  for (cluster in all_clusters) {
    log_msg("  -", cluster)
  }

  # Get full marker dataframe
  marker_df <- get_full_marker_dataframe()
  log_msg("Loaded marker data:", nrow(marker_df), "rows")

  output_name <- file.path(results$outputs, "batch_results")
  log_msg("Results will be saved to:", results$base)

  # Run the test
  start_time <- Sys.time()
  errors <- list()
  status <- "error"

  tryCatch({
    log_msg("Running runCASSIA_batch via R package...")

    # Call CASSIA R function
    CASSIA::runCASSIA_batch(
      marker = marker_df,
      output_name = output_name,
      n_genes = data_config$n_genes %||% 30,
      model = llm_config$model %||% "google/gemini-2.5-flash",
      temperature = llm_config$temperature %||% 0.3,
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      max_workers = llm_config$max_workers %||% 3,
      provider = llm_config$provider %||% "openrouter",
      validator_involvement = config$validator$default %||% "v1"
    )

    # Check output files
    full_csv <- paste0(output_name, "_full.csv")
    summary_csv <- paste0(output_name, "_summary.csv")

    if (file.exists(full_csv)) {
      results_df <- read.csv(full_csv)
      clusters_annotated <- nrow(results_df)

      log_msg("Batch Results:")
      log_msg("  Clusters annotated:", clusters_annotated, "/", length(all_clusters))
      log_msg("  Output files created:")
      log_msg("    -", basename(full_csv))
      if (file.exists(summary_csv)) {
        log_msg("    -", basename(summary_csv))
      }

      if (clusters_annotated == length(all_clusters)) {
        status <- "passed"
      } else {
        status <- "failed"
        errors <- list(paste("Only", clusters_annotated, "/", length(all_clusters), "clusters annotated"))
      }
    } else {
      status <- "failed"
      errors <- list("Output file not created")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    log_error(e)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata
  metadata <- create_test_metadata(
    test_name = "batch_annotation",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(all_clusters),
    errors = errors
  )
  save_test_metadata(results$outputs, metadata)

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  stop_logging()
  return(success)
}

# Run test
success <- run_batch_annotation_test()
quit(status = if (success) 0 else 1)
