# CASSIA Test 04: Quality Scoring (R)
# ====================================
# Tests the runCASSIA_score_batch function to score annotation results.
#
# Usage:
#     Rscript test_quality_scoring.R
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
source(file.path(script_dir, "..", "shared", "r", "logging_manager.R"))

run_quality_scoring_test <- function() {
  print_test_header("04 - Quality Scoring (R)")

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

  # Create results directory (BEFORE any logging)
  results <- create_results_dir("04_quality_scoring", get_test_mode())
  start_logging(results$logs)

  log_msg("Results will be saved to:", results$base)

  # Check for existing batch results from Test 02
  batch_results_dir <- get_latest_results("02_batch_annotation")
  input_file <- NULL

  if (!is.null(batch_results_dir)) {
    potential_file <- file.path(batch_results_dir, "batch_results_full.csv")
    if (file.exists(potential_file)) {
      input_file <- potential_file
      log_msg("Using existing batch results:", input_file)
    }
  }

  # If no existing results, run a quick batch annotation
  if (is.null(input_file)) {
    log_msg("No existing batch results found. Running batch annotation first...")
    marker_df <- get_full_marker_dataframe()
    batch_output <- file.path(results$outputs, "batch_for_scoring")

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

  output_file <- file.path(results$outputs, "scored_results.csv")

  tryCatch({
    log_msg("Running quality scoring...")
    log_msg("Input:", input_file)
    log_msg("Output:", output_file)

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

        log_msg("Scoring Results:")
        log_msg("  Clusters scored:", length(scores))
        log_msg("  Average score:", round(mean(scores), 1))
        log_msg("  Score range:", min(scores), "-", max(scores))

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
    log_error(e)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  all_clusters <- get_all_clusters()
  metadata <- create_test_metadata(
    test_name = "quality_scoring",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(all_clusters),
    errors = errors
  )
  save_test_metadata(results$outputs, metadata)

  save_test_results(results$outputs, list(
    input_file = input_file,
    output_file = output_file,
    scoring_results = scoring_results
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  stop_logging()
  return(success)
}

# Run test
success <- run_quality_scoring_test()
quit(status = if (success) 0 else 1)
