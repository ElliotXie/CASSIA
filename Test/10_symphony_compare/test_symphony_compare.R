# CASSIA Test 10: Symphony Compare (R)
# =====================================
# Tests the symphonyCompare function via R package for multi-model
# cell type comparison with AI consensus building.
#
# Usage:
#     Rscript test_symphony_compare.R
#
# Functions tested:
# - symphonyCompare(): Orchestrate multiple AI models to compare cell types

# Get script directory
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") {
  script_dir <- "."
}

# Source shared utilities
source(file.path(script_dir, "..", "shared", "r", "test_utils.R"))
source(file.path(script_dir, "..", "shared", "r", "fixtures.R"))
source(file.path(script_dir, "..", "shared", "r", "result_manager.R"))

run_symphony_compare_test <- function() {
  print_test_header("10 - Symphony Compare (R)")

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

  # Create results directory
  results_dir <- create_results_dir("10_symphony_compare")
  cat("Results will be saved to:", results_dir, "\n")

  # Test parameters
  test_cluster <- "plasma cell"
  markers <- get_cluster_markers(test_cluster)
  marker_set <- paste(markers[1:15], collapse = ", ")

  # Cell types to compare
  celltypes <- c("Plasma cell", "B cell", "T cell")

  # Run tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  symphony_results <- list()

  tryCatch({
    cat("\n--- Test: symphonyCompare ---\n")
    cat("  Tissue:", data_config$tissue %||% "large intestine", "\n")
    cat("  Species:", data_config$species %||% "human", "\n")
    cat("  Cell types to compare:", paste(celltypes, collapse = ", "), "\n")
    cat("  Marker set:", substr(marker_set, 1, 50), "...\n")
    cat("  Model preset: budget (cost-effective for testing)\n")
    cat("  Discussion: Disabled (for faster testing)\n")

    result <- CASSIA::symphonyCompare(
      tissue = data_config$tissue %||% "large intestine",
      celltypes = celltypes,
      marker_set = marker_set,
      species = data_config$species %||% "human",
      model_preset = "budget",
      output_dir = results_dir,
      output_basename = "symphony_test",
      enable_discussion = FALSE,
      max_discussion_rounds = 0,
      consensus_threshold = 0.6,
      generate_report = TRUE,
      verbose = TRUE
    )

    # Check result
    if (is.list(result)) {
      cat("\n--- Symphony Compare Results ---\n")
      cat("  Consensus:", result$consensus %||% "No consensus", "\n")
      cat("  Confidence:", sprintf("%.1f%%", (result$confidence %||% 0) * 100), "\n")
      cat("  CSV file:", basename(result$csv_file %||% "N/A"), "\n")
      cat("  HTML report:", basename(result$html_file %||% "N/A"), "\n")

      # Summary statistics
      summary_data <- result$summary
      if (!is.null(summary_data)) {
        cat("\n  Summary:\n")
        cat("    Models used:", summary_data$models_used %||% "N/A", "\n")
        cat("    Total rounds:", summary_data$total_rounds %||% "N/A", "\n")
        cat("    Consensus reached:", summary_data$consensus_reached %||% FALSE, "\n")
      }

      symphony_results <- list(
        consensus = result$consensus,
        confidence = result$confidence,
        csv_file = result$csv_file,
        html_file = result$html_file,
        summary = summary_data
      )

      # Validate results
      if (!is.null(result$csv_file) && file.exists(result$csv_file)) {
        status <- "passed"
      } else {
        status <- "failed"
        errors <- list("CSV file not created")
      }
    } else {
      status <- "failed"
      errors <- list("Unexpected result format")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    cat("\nError:", e$message, "\n")
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "symphony_compare",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(celltypes),
    errors = errors
  )
  save_test_metadata(results_dir, metadata)

  save_test_results(results_dir, list(
    celltypes_compared = celltypes,
    marker_set = marker_set,
    model_preset = "budget",
    results = symphony_results
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  return(success)
}

# Run test
success <- run_symphony_compare_test()
quit(status = if (success) 0 else 1)
