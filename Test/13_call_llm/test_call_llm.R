# CASSIA Test 13: call_llm Provider Testing (R)
# ===============================================
# Tests that call_llm correctly connects to different LLM providers.
#
# This test verifies API connectivity for:
# - OpenAI (gpt-4o)
# - Anthropic (claude-sonnet-4-5)
# - OpenRouter (google/gemini-2.5-flash)
# - Custom/DeepSeek (deepseek-chat via https://api.deepseek.com)
#
# Usage:
#     Rscript test_call_llm.R

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

run_call_llm_test <- function() {
  print_test_header("13 - call_llm Provider Testing (R)")

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

  # Create results directory
  results_dir <- create_results_dir("13_call_llm")

  start_logging(results_dir)
  log_msg("Results will be saved to:", results_dir)

  # Import Python call_llm function via reticulate
  cassia_module <- reticulate::import("CASSIA")
  call_llm <- cassia_module$call_llm

  # Test cases: list(provider, model, env_var_name, display_name)
  test_cases <- list(
    list(provider = "openai", model = "gpt-4o", env_var = "OPENAI_API_KEY", display_name = "OpenAI"),
    list(provider = "anthropic", model = "claude-sonnet-4-5", env_var = "ANTHROPIC_API_KEY", display_name = "Anthropic"),
    list(provider = "openrouter", model = "google/gemini-2.5-flash", env_var = "OPENROUTER_API_KEY", display_name = "OpenRouter"),
    list(provider = "https://api.deepseek.com", model = "deepseek-chat", env_var = "CUSTOMIZED_API_KEY", display_name = "Custom (DeepSeek)")
  )

  # Simple test prompt
  test_prompt <- "Reply with exactly one word: 'success'"

  start_time <- Sys.time()
  errors <- list()
  test_results <- list()

  passed_count <- 0
  skipped_count <- 0
  failed_count <- 0
  total_tests <- length(test_cases)

  for (tc in test_cases) {
    log_separator()
    log_msg("Provider:", tc$display_name)
    log_msg("  Model:", tc$model)
    log_separator()

    test_result <- list(
      provider = tc$provider,
      display_name = tc$display_name,
      model = tc$model,
      env_var = tc$env_var,
      status = "pending",
      response_length = 0,
      error = NULL
    )

    # Check if API key is set
    api_key <- Sys.getenv(tc$env_var, unset = "")
    if (api_key == "") {
      log_msg("  Status: [SKIP]", tc$env_var, "not set in environment")
      test_result$status <- "skipped"
      skipped_count <- skipped_count + 1
      test_results <- c(test_results, list(test_result))
      next
    }

    log_msg("  Status: API key found")
    log_msg("  Sending test prompt...")

    tryCatch({
      response <- call_llm(
        prompt = test_prompt,
        provider = tc$provider,
        model = tc$model,
        temperature = 0.3,
        max_tokens = as.integer(50)
      )

      if (!is.null(response) && nchar(response) > 0) {
        test_result$status <- "passed"
        test_result$response_length <- nchar(response)
        test_result$response_preview <- substr(response, 1, 100)
        log_msg("  [OK] Response received (", nchar(response), " chars)", sep = "")
        log_msg("  Response:", substr(response, 1, 80), if (nchar(response) > 80) "..." else "")
        passed_count <- passed_count + 1
      } else {
        test_result$status <- "failed"
        test_result$error <- "Empty response received"
        log_msg("  [X] FAILED: Empty response received")
        errors <- c(errors, paste0(tc$display_name, ": Empty response"))
        failed_count <- failed_count + 1
      }
    }, error = function(e) {
      test_result$status <<- "failed"
      test_result$error <<- e$message
      error_msg <- substr(e$message, 1, 200)
      log_msg("  [X] ERROR:", error_msg)
      errors <<- c(errors, paste0(tc$display_name, ": ", substr(e$message, 1, 100)))
      failed_count <<- failed_count + 1
    })

    test_results <- c(test_results, list(test_result))
  }

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Summary
  log_separator()
  log_msg("SUMMARY")
  log_separator()
  log_msg("  Passed: ", passed_count, "/", total_tests, "", sep = "")
  log_msg("  Failed: ", failed_count, "/", total_tests, "", sep = "")
  log_msg("  Skipped:", skipped_count, "/", total_tests, "", sep = "")
  log_msg("  Duration:", round(duration, 2), "s")

  # Determine overall status
  tests_run <- total_tests - skipped_count
  if (tests_run == 0) {
    status <- "skipped"
    log_msg("\n  Note: All tests skipped (no API keys configured)")
  } else if (passed_count == tests_run) {
    status <- "passed"
  } else {
    status <- "failed"
  }

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "call_llm_provider_test",
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

  save_test_results(results_dir, list(
    provider_tests = test_results
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_call_llm_test()
quit(status = if (success) 0 else 1)
