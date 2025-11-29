# CASSIA Test 04: Quality Scoring (R) - INSTALL MODE
# ===================================================
# Tests the runCASSIA_score_batch function to score annotation results.
# Uses devtools::install_local() for full package installation testing.
#
# Usage:
#     Rscript test_quality_scoring_install.R
#
# Note: This test requires batch annotation results from Test 02.

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

run_quality_scoring_test <- function() {
  print_test_header("04 - Quality Scoring (R) [INSTALL MODE]")

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
  results_dir <- create_results_dir("04_quality_scoring")
  cat("Results will be saved to:", results_dir, "\n")

  # Check for existing batch results from Test 02
  batch_results_dir <- get_latest_results("02_batch_annotation")
  input_file <- NULL

  if (!is.null(batch_results_dir)) {
    potential_file <- file.path(batch_results_dir, "batch_results_full.csv")
    if (file.exists(potential_file)) {
      input_file <- potential_file
      cat("\nUsing existing batch results:", input_file, "\n")
    }
  }

  # If no existing results, run a quick batch annotation
  if (is.null(input_file)) {
    cat("\nNo existing batch results found. Running batch annotation first...\n")
    marker_df <- get_full_marker_dataframe()
    batch_output <- file.path(results_dir, "batch_for_scoring")

    CASSIA::runCASSIA_batch(
      marker = marker_df,
      output_name = batch_output,
      n_genes = data_config$n_genes %||% 30,
      model = llm_config$model %||% "google/gemini-2.5-flash",
      temperature = llm_config$temperature %||% 0.3,
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      max_workers = llm_config$max_workers %||% 3,
      provider = llm_config$provider %||% "openrouter",
      validator_involvement = config$validator$default %||% "v1"
    )
    input_file <- paste0(batch_output, "_full.csv")
  }

  # Run quality scoring
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  scoring_results <- list()

  output_file <- file.path(results_dir, "scored_results.csv")

  tryCatch({
    cat("\nRunning quality scoring (install mode)...\n")
    cat("Input:", input_file, "\n")
    cat("Output:", output_file, "\n")

    CASSIA::runCASSIA_score_batch(
      input_file = input_file,
      output_file = output_file,
      max_workers = llm_config$max_workers %||% 3,
      model = llm_config$model %||% "google/gemini-2.5-flash",
      provider = llm_config$provider %||% "openrouter"
    )

    # Check results
    if (file.exists(output_file)) {
      scored_df <- read.csv(output_file)

      # Check for score column
      score_col <- NULL
      for (col in c("score", "Score", "SCORE")) {
        if (col %in% names(scored_df)) {
          score_col <- col
          break
        }
      }

      if (!is.null(score_col)) {
        scores <- scored_df[[score_col]]
        scores <- scores[!is.na(scores)]

        scoring_results <- list(
          clusters_scored = length(scores),
          avg_score = mean(scores),
          min_score = min(scores),
          max_score = max(scores),
          scores = as.list(scores)
        )

        cat("\nScoring Results:\n")
        cat("  Clusters scored:", length(scores), "\n")
        cat("  Average score:", round(mean(scores), 1), "\n")
        cat("  Score range:", min(scores), "-", max(scores), "\n")

        # Validate scores are in expected range
        if (min(scores) >= 0 && max(scores) <= 100) {
          status <- "passed"
        } else {
          status <- "failed"
          errors <- list("Scores outside expected range (0-100)")
        }
      } else {
        status <- "failed"
        errors <- list("No score column found in output")
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

  # Save metadata and results
  all_clusters <- get_all_clusters()
  metadata <- create_test_metadata(
    test_name = "quality_scoring_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(all_clusters),
    errors = errors
  )
  save_test_metadata(results_dir, metadata)

  save_test_results(results_dir, list(
    input_file = input_file,
    output_file = output_file,
    scoring_results = scoring_results,
    mode = "install"
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  return(success)
}

# Run test
success <- run_quality_scoring_test()
quit(status = if (success) 0 else 1)
