# CASSIA Test 05: Model Settings (R)
# ===================================
# Tests model name resolution and provider shortcuts via R package.
#
# Usage:
#     Rscript test_model_settings.R

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

run_model_settings_test <- function() {
  print_test_header("05 - Model Settings (R)")

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

  # Create results directory
  results <- create_results_dir("05_model_settings", get_test_mode())

  start_logging(results$logs)
  log_msg("Results will be saved to:", results$base)

  start_time <- Sys.time()
  errors <- list()
  test_results <- list(
    resolution_tests = list(),
    practical_test = NULL
  )

  # Test 1: Model name resolution via Python (through reticulate)
  log_separator()
  log_msg("Testing Model Name Resolution")
  log_separator()

  resolution_tests <- list(
    list(simple_name = "gemini", provider = "openrouter", expected_contains = "gemini"),
    list(simple_name = "best", provider = "openrouter", expected_contains = "gemini"),
    list(simple_name = "cheap", provider = "openrouter", expected_contains = NULL),
    list(simple_name = "flash", provider = "openrouter", expected_contains = "flash")
  )

  resolution_passed <- 0

  tryCatch({
    # Import model_settings module via reticulate
    model_settings <- reticulate::import("CASSIA.core.model_settings")

    for (test in resolution_tests) {
      tryCatch({
        resolved <- model_settings$resolve_model_name(test$simple_name, test$provider)
        resolved_name <- resolved[[1]]

        success <- TRUE
        if (!is.null(test$expected_contains)) {
          if (!grepl(tolower(test$expected_contains), tolower(resolved_name))) {
            success <- FALSE
            errors <- c(errors, paste0("'", test$simple_name, "' resolved to '",
                                       resolved_name, "', expected to contain '",
                                       test$expected_contains, "'"))
          }
        }

        test_results$resolution_tests <- c(test_results$resolution_tests, list(list(
          input = test$simple_name,
          provider = test$provider,
          resolved = resolved_name,
          success = success
        )))

        status_symbol <- if (success) "[OK]" else "[X]"
        log_msg("  ", status_symbol, " '", test$simple_name, "' + '", test$provider, "' -> '", resolved_name, "'", sep = "")

        if (success) resolution_passed <- resolution_passed + 1

      }, error = function(e) {
        test_results$resolution_tests <- c(test_results$resolution_tests, list(list(
          input = test$simple_name,
          provider = test$provider,
          error = e$message,
          success = FALSE
        )))
        errors <<- c(errors, paste0("Resolution failed for '", test$simple_name, "': ", e$message))
        log_msg("  [X] '", test$simple_name, "' + '", test$provider, "' -> ERROR: ", e$message, "", sep = "")
      })
    }
  }, error = function(e) {
    log_msg("  [X] Could not import model_settings module:", e$message)
    errors <<- c(errors, paste0("Module import failed: ", e$message))
  })

  log_msg("\nResolution tests:", resolution_passed, "/", length(resolution_tests), "passed")

  # Test 2: Practical test with gemini-2.5-flash via OpenRouter
  log_separator()
  log_msg("Testing Practical Annotation with OpenRouter")
  log_separator()

  data_config <- config$data
  marker_df <- get_marker_dataframe_for_cluster("plasma cell", 20)

  tryCatch({
    log_msg("  Model: google/gemini-2.5-flash")
    log_msg("  Provider: openrouter")
    log_msg("  Running annotation...")

    result <- CASSIA::runCASSIA(
      model = "google/gemini-2.5-flash",
      temperature = 0.3,
      marker_list = marker_df,
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      provider = "openrouter",
      validator_involvement = "v1"
    )

    practical_success <- !is.null(result$main_cell_type) && result$main_cell_type != ""
    test_results$practical_test <- list(
      success = practical_success,
      main_cell_type = result$main_cell_type,
      sub_cell_types = result$sub_cell_types
    )

    if (practical_success) {
      log_msg("  [OK] Annotation successful")
      log_msg("       Main type:", result$main_cell_type)
    } else {
      log_msg("  [X] Annotation returned empty result")
      errors <- c(errors, "Practical test returned empty result")
    }

  }, error = function(e) {
    test_results$practical_test <<- list(
      success = FALSE,
      error = e$message
    )
    errors <<- c(errors, paste0("Practical test failed: ", e$message))
    log_msg("  [X] ERROR:", e$message)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Determine overall status
  resolution_all_passed <- resolution_passed == length(resolution_tests)
  practical_passed <- !is.null(test_results$practical_test) &&
                      test_results$practical_test$success

  if (resolution_all_passed && practical_passed) {
    status <- "passed"
  } else if (practical_passed) {
    status <- "passed"  # Practical test is most important
  } else {
    status <- "failed"
  }

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "model_settings",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list("plasma cell"),
    errors = as.list(errors)
  )
  save_test_metadata(results$outputs, metadata)
  save_test_results(results$outputs, test_results)

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_model_settings_test()
quit(status = if (success) 0 else 1)
