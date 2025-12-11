# CASSIA Test 14: Customized Provider Testing (R)
# ================================================
# Tests that call_llm correctly connects to customized/custom API providers.
#
# This test demonstrates how to use a custom OpenAI-compatible API endpoint
# with CASSIA, using DeepSeek as an example.
#
# Usage:
#     Rscript test_customized_provider.R

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
source(file.path(script_dir, "..", "shared", "r", "result_manager.R"))
source(file.path(script_dir, "..", "shared", "r", "logging_manager.R"))

run_customized_provider_test <- function() {
  print_test_header("14 - Customized Provider Testing (R)")

  # Load configuration
  config <- load_config()
  print_config_summary(config)

  # Setup API keys
  setup_api_keys()

  # Setup CASSIA
  tryCatch({
    setup_cassia()
    message("CASSIA R package loaded successfully")
  }, error = function(e) {
    message("Error loading CASSIA: ", e$message)
    return(FALSE)
  })

  # Import CASSIA Python module
  tryCatch({
    CASSIA <- import_cassia()
    message("CASSIA Python module imported successfully")
  }, error = function(e) {
    message("Error importing CASSIA: ", e$message)
    return(FALSE)
  })

  # Create results directory
  results_dirs <- create_results_dir("14_customized_provider")
  results_dir <- results_dirs$base
  message(sprintf("Results will be saved to: %s", results_dir))

  # Test cases for customized providers: list(provider_url, model, env_var_name, display_name)
  # DeepSeek is used as an example of a custom OpenAI-compatible API
  test_cases <- list(
    list(
      provider = "https://api.deepseek.com",
      model = "deepseek-chat",
      env_var = "CUSTOMIZED_API_KEY",
      display_name = "DeepSeek Custom API"
    )
  )

  # Simple test prompt
  test_prompt <- "Reply with exactly one word: 'success'"

  start_time <- Sys.time()
  errors <- list()
  test_results <- list(provider_tests = list())

  passed_count <- 0
  skipped_count <- 0
  failed_count <- 0
  total_tests <- length(test_cases)

  for (test_case in test_cases) {
    provider <- test_case$provider
    model <- test_case$model
    env_var <- test_case$env_var
    display_name <- test_case$display_name

    cat("\n", strrep("=", 50), "\n")
    cat(sprintf("Provider: %s\n", display_name))
    cat(sprintf("  Endpoint: %s\n", provider))
    cat(sprintf("  Model: %s\n", model))
    cat(strrep("=", 50), "\n")

    test_result <- list(
      provider = provider,
      display_name = display_name,
      model = model,
      env_var = env_var,
      status = "pending",
      response_length = 0,
      error = NULL
    )

    # Check if API key is set
    api_key <- Sys.getenv(env_var)
    if (api_key == "") {
      cat(sprintf("  Status: [SKIP] %s not set in environment\n", env_var))
      test_result$status <- "skipped"
      skipped_count <- skipped_count + 1
      test_results$provider_tests <- append(test_results$provider_tests, list(test_result))
      next
    }

    cat("  Status: API key found\n")
    cat("  Sending test prompt...\n")

    tryCatch({
      response <- CASSIA$call_llm(
        prompt = test_prompt,
        provider = provider,
        model = model,
        temperature = 0.3,
        max_tokens = 50L
      )

      if (!is.null(response) && nchar(response) > 0) {
        test_result$status <- "passed"
        test_result$response_length <- nchar(response)
        test_result$response_preview <- if (nchar(response) > 100) substr(response, 1, 100) else response
        cat(sprintf("  [OK] Response received (%d chars)\n", nchar(response)))
        display_response <- if (nchar(response) > 80) paste0(substr(response, 1, 80), "...") else response
        cat(sprintf("  Response: %s\n", display_response))
        passed_count <- passed_count + 1
      } else {
        test_result$status <- "failed"
        test_result$error <- "Empty response received"
        cat("  [X] FAILED: Empty response received\n")
        errors <- append(errors, sprintf("%s: Empty response", display_name))
        failed_count <- failed_count + 1
      }

    }, error = function(e) {
      test_result$status <<- "failed"
      test_result$error <<- e$message
      error_msg <- e$message
      # Truncate long error messages for display
      if (nchar(error_msg) > 200) {
        error_msg <- paste0(substr(error_msg, 1, 200), "...")
      }
      cat(sprintf("  [X] ERROR: %s\n", error_msg))
      errors <<- append(errors, sprintf("%s: %s", display_name, substr(e$message, 1, 100)))
      failed_count <<- failed_count + 1
    })

    test_results$provider_tests <- append(test_results$provider_tests, list(test_result))
  }

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Summary
  cat("\n", strrep("=", 50), "\n")
  cat("SUMMARY\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("  Passed:  %d/%d\n", passed_count, total_tests))
  cat(sprintf("  Failed:  %d/%d\n", failed_count, total_tests))
  cat(sprintf("  Skipped: %d/%d\n", skipped_count, total_tests))
  cat(sprintf("  Duration: %.2fs\n", duration))

  # Determine overall status
  # Test passes if all non-skipped tests passed
  tests_run <- total_tests - skipped_count
  if (tests_run == 0) {
    status <- "skipped"
    cat("\n  Note: All tests skipped (no API keys configured)\n")
  } else if (passed_count == tests_run) {
    status <- "passed"
  } else {
    status <- "failed"
  }

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "customized_provider_test_r",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(),  # Not applicable for this test
    errors = errors
  )
  metadata$test_summary <- list(
    passed = passed_count,
    failed = failed_count,
    skipped = skipped_count,
    total = total_tests
  )
  save_test_metadata(results_dir, metadata)
  save_test_results(results_dir, test_results)

  # Print final result
  success <- status == "passed"
  print_test_result(success, sprintf("Duration: %.2fs", duration))

  return(success)
}

# Run the test
success <- tryCatch({
  run_customized_provider_test()
}, error = function(e) {
  message("Test failed with error: ", e$message)
  FALSE
})

# Exit with appropriate code
quit(status = if (success) 0 else 1)
