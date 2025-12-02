# CASSIA Test 03: Validator Comparison (R)
# =========================================
# Compares v0 (strict) and v1 (moderate) validators via R package.
#
# Usage:
#     Rscript test_validator_comparison.R

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

run_validator_comparison_test <- function() {
  print_test_header("03 - Validator Comparison (R)")

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

  # Get marker data
  marker_df <- get_full_marker_dataframe()
  all_clusters <- get_all_clusters()

  # Create results directory (BEFORE any logging)
  results <- create_results_dir("03_validator_comparison", get_test_mode())
  start_logging(results$logs)

  log_msg("Comparing validators on", length(all_clusters), "clusters")
  log_msg("Results will be saved to:", results$base)

  # Run with both validators
  validators <- c("v0", "v1")
  validator_results <- list()
  errors <- list()

  start_time <- Sys.time()

  for (validator in validators) {
    log_separator()
    log_msg("Testing", validator, "validator",
        if (validator == "v0") "(strict)" else "(moderate)")
    log_separator()

    output_name <- file.path(results$outputs, paste0("results_", validator))

    tryCatch({
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
        validator_involvement = validator
      )

      # Load results
      full_csv <- paste0(output_name, "_full.csv")
      if (file.exists(full_csv)) {
        df <- read.csv(full_csv)
        validator_results[[validator]] <- list(
          clusters_annotated = nrow(df),
          results_file = full_csv
        )
        log_msg("  Annotated", nrow(df), "clusters")
      } else {
        errors <- c(errors, paste0(validator, ": Output file not created"))
      }

    }, error = function(e) {
      errors <<- c(errors, paste0(validator, ": ", e$message))
      log_error(e)
    })
  }

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Compare results
  log_separator()
  log_msg("Comparison Results")
  log_separator()

  comparison <- list()
  status <- "failed"

  if ("v0" %in% names(validator_results) && "v1" %in% names(validator_results)) {
    v0_count <- validator_results$v0$clusters_annotated
    v1_count <- validator_results$v1$clusters_annotated
    comparison <- list(
      v0_clusters = v0_count,
      v1_clusters = v1_count,
      both_complete = (v0_count == length(all_clusters) && v1_count == length(all_clusters))
    )
    log_msg("  v0 (strict):", v0_count, "clusters")
    log_msg("  v1 (moderate):", v1_count, "clusters")
    status <- if (comparison$both_complete) "passed" else "failed"
  } else {
    log_msg("  Could not compare - one or both validators failed")
  }

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "validator_comparison",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(all_clusters),
    errors = as.list(errors)
  )
  save_test_metadata(results$outputs, metadata)

  save_test_results(results$outputs, list(
    validator_results = validator_results,
    comparison = comparison
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  stop_logging()
  return(success)
}

# Run test
success <- run_validator_comparison_test()
quit(status = if (success) 0 else 1)
