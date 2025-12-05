# CASSIA Test 17: API Key Validation Testing (R)
# ================================================
# Tests the validate_api_keys() function for different LLM providers.
#
# This test verifies:
# - API key validation for OpenAI, Anthropic, and OpenRouter
# - Caching behavior (second validation should be instant)
# - Force revalidation
# - Error handling for invalid/missing keys
#
# Usage:
#     Rscript test_api_validation.R

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

run_api_validation_test <- function() {
  print_test_header("17 - API Key Validation Testing (R)")

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

  # Import CASSIA
  tryCatch({
    CASSIA <- import_cassia()
    message("CASSIA Python module imported successfully")
  }, error = function(e) {
    message("Error importing CASSIA: ", e$message)
    return(FALSE)
  })

  # Create results directory
  results_dir <- create_results_dir("17_api_validation")
  message(sprintf("Results will be saved to: %s", results_dir))

  start_time <- Sys.time()
  errors <- list()
  all_passed <- TRUE

  # Test 1: Validate all providers
  cat("\n", strrep("=", 60), "\n")
  cat("Test 1: Validate All Providers\n")
  cat(strrep("=", 60), "\n")

  tryCatch({
    result <- CASSIA$validate_api_keys(verbose=TRUE)
    message("Validation result: ", toString(result))

    if (is.null(result) || (!is.logical(result) && length(result) == 0)) {
      message("  [SKIP] No API keys configured")
    } else if (is.logical(result)) {
      if (result) {
        message("  [OK] API key is valid")
      } else {
        message("  [FAIL] API key is invalid")
        all_passed <- FALSE
      }
    } else {
      # Dictionary result
      message("  [OK] Validated multiple providers")
    }
  }, error = function(e) {
    message("  [ERROR] ", e$message)
    errors <- append(errors, sprintf("validate_all: %s", e$message))
    all_passed <- FALSE
  })

  # Test 2: Validate specific provider
  cat("\n", strrep("=", 60), "\n")
  cat("Test 2: Validate Specific Provider (OpenAI)\n")
  cat(strrep("=", 60), "\n")

  tryCatch({
    if (Sys.getenv("OPENAI_API_KEY") != "") {
      result <- CASSIA$validate_api_keys(provider="openai", verbose=TRUE)
      message("OpenAI validation: ", result)

      if (result) {
        message("  [OK] OpenAI key is valid")
      } else {
        message("  [FAIL] OpenAI key is invalid")
        all_passed <- FALSE
      }
    } else {
      message("  [SKIP] OPENAI_API_KEY not set")
    }
  }, error = function(e) {
    message("  [ERROR] ", e$message)
    errors <- append(errors, sprintf("validate_openai: %s", e$message))
    all_passed <- FALSE
  })

  # Test 3: Test caching behavior
  cat("\n", strrep("=", 60), "\n")
  cat("Test 3: Caching Behavior\n")
  cat(strrep("=", 60), "\n")

  tryCatch({
    if (Sys.getenv("OPENAI_API_KEY") != "") {
      # Clear cache
      CASSIA$clear_validation_cache(provider="openai")

      # First validation
      message("First validation (should make API call)...")
      start1 <- Sys.time()
      result1 <- CASSIA$validate_api_keys(provider="openai", verbose=FALSE)
      duration1 <- as.numeric(difftime(Sys.time(), start1, units="secs"))
      message(sprintf("  Duration: %.3f seconds", duration1))

      # Second validation (should use cache)
      message("Second validation (should use cache)...")
      start2 <- Sys.time()
      result2 <- CASSIA$validate_api_keys(provider="openai", verbose=FALSE)
      duration2 <- as.numeric(difftime(Sys.time(), start2, units="secs"))
      message(sprintf("  Duration: %.6f seconds", duration2))

      if (duration2 < 0.01) {
        message("  [OK] Cache is working (second call < 10ms)")
      } else {
        message("  [WARN] Cache may not be working optimally")
      }
    } else {
      message("  [SKIP] OPENAI_API_KEY not set")
    }
  }, error = function(e) {
    message("  [ERROR] ", e$message)
    errors <- append(errors, sprintf("caching: %s", e$message))
  })

  # Test 4: Invalid API key
  cat("\n", strrep("=", 60), "\n")
  cat("Test 4: Invalid API Key\n")
  cat(strrep("=", 60), "\n")

  tryCatch({
    fake_key <- "sk-fake-key-for-testing-12345"
    result <- CASSIA$validate_api_keys(provider="openai", api_key=fake_key, verbose=TRUE)

    if (!result) {
      message("  [OK] Invalid key correctly detected")
    } else {
      message("  [FAIL] Invalid key reported as valid!")
      all_passed <- FALSE
    }
  }, error = function(e) {
    # Exception is also acceptable for invalid key
    message("  [OK] Exception raised for invalid key")
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units="secs"))

  # Summary
  cat("\n", strrep("=", 60), "\n")
  cat("SUMMARY\n")
  cat(strrep("=", 60), "\n")
  cat(sprintf("  Status: %s\n", if (all_passed) "PASSED" else "FAILED"))
  cat(sprintf("  Duration: %.2f seconds\n", duration))
  cat(sprintf("  Errors: %d\n", length(errors)))

  # Save metadata
  metadata <- create_test_metadata(
    test_name = "api_validation_test_r",
    config = config,
    duration_seconds = duration,
    status = if (all_passed) "passed" else "failed",
    clusters_tested = list(),
    errors = errors
  )
  save_test_metadata(results_dir, metadata)

  # Print final result
  print_test_result(all_passed, sprintf("Duration: %.2f seconds", duration))

  return(all_passed)
}

# Run the test
success <- tryCatch({
  run_api_validation_test()
}, error = function(e) {
  message("Test failed with error: ", e$message)
  FALSE
})

# Exit with appropriate code
quit(status = if (success) 0 else 1)
