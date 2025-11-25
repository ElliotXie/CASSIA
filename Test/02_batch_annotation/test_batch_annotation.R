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
  cat("\nTesting batch annotation for", length(all_clusters), "clusters:\n")
  for (cluster in all_clusters) {
    cat("  -", cluster, "\n")
  }

  # Get full marker dataframe
  marker_df <- get_full_marker_dataframe()
  cat("\nLoaded marker data:", nrow(marker_df), "rows\n")

  # Create results directory
  results_dir <- create_results_dir("02_batch_annotation")
  output_name <- file.path(results_dir, "batch_results")
  cat("Results will be saved to:", results_dir, "\n")

  # Run the test
  start_time <- Sys.time()
  errors <- list()
  status <- "error"

  tryCatch({
    cat("\nRunning runCASSIA_batch via R package...\n")

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

      cat("\nBatch Results:\n")
      cat("  Clusters annotated:", clusters_annotated, "/", length(all_clusters), "\n")
      cat("  Output files created:\n")
      cat("    -", basename(full_csv), "\n")
      if (file.exists(summary_csv)) {
        cat("    -", basename(summary_csv), "\n")
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
    cat("\nError:", e$message, "\n")
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
  save_test_metadata(results_dir, metadata)

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  return(success)
}

# Run test
success <- run_batch_annotation_test()
quit(status = if (success) 0 else 1)
