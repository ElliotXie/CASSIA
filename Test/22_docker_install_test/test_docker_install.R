#!/usr/bin/env Rscript
# =============================================================================
# CASSIA Test 22: Docker Installation Test (R Wrapper)
# =============================================================================
# This R script wraps the Python Docker test orchestrator.
# It checks Docker availability and calls the Python script.
#
# Usage:
#   Rscript test_docker_install.R
# =============================================================================

# Get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
test_root <- dirname(script_dir)

# Source shared utilities
source(file.path(test_root, "shared", "r", "test_utils.R"))
source(file.path(test_root, "shared", "r", "result_manager.R"))
source(file.path(test_root, "shared", "r", "logging_manager.R"))

# =============================================================================
# Main test function
# =============================================================================

run_docker_test <- function() {
  print_test_header("22 - Docker Installation Test (R Wrapper)")

  # Check Docker availability
  docker_check <- tryCatch({
    result <- system2("docker", "info", stdout = TRUE, stderr = TRUE)
    TRUE
  }, error = function(e) {
    FALSE
  }, warning = function(w) {
    FALSE
  })

  if (!docker_check) {
    cat("\nERROR: Docker is not available or not running.\n")
    cat("Please start Docker Desktop and try again.\n")
    return(FALSE)
  }

  cat("\nDocker is available.\n")

  # Load configuration
  config <- load_config()

  # Setup API keys
  setup_api_keys()

  # Create results directory
  results <- create_results_dir("22_docker_install_test", get_test_mode())
  cat(sprintf("Results will be saved to: %s\n", results$base))

  # Start logging
  start_logging(results$logs)

  # Run Python orchestrator
  start_time <- Sys.time()
  errors <- list()
  test_results <- list(success = FALSE)

  tryCatch({
    script_path <- file.path(script_dir, "scripts", "run_scenarios.py")

    # Build command
    cmd_args <- c(
      script_path,
      "--output", results$base
    )

    # Add API key if available
    api_key <- Sys.getenv("OPENROUTER_API_KEY", "")
    if (nchar(api_key) > 0) {
      cmd_args <- c(cmd_args, "--api-key", api_key)
      cat("\nAPI key found - will include annotation test\n")
    } else {
      cat("\nNo API key - skipping annotation test\n")
    }

    cat("\nRunning Docker scenarios via Python orchestrator...\n")
    cat("This may take 10-30 minutes on first run (building images)\n")
    cat(strrep("-", 50), "\n")

    # Run Python script
    exit_code <- system2("python", cmd_args)

    test_results$exit_code <- exit_code
    test_results$success <- (exit_code == 0)

  }, error = function(e) {
    errors <<- c(errors, paste("Test failed:", e$message))
    test_results$success <<- FALSE
  })

  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Stop logging
  stop_logging()

  # Save metadata
  metadata <- create_test_metadata(
    test_name = "22_docker_install_test",
    config = config,
    duration_seconds = duration,
    status = if(test_results$success) "passed" else "failed",
    clusters_tested = list(),
    errors = errors
  )
  save_test_metadata(results$base, metadata)

  # Save results
  save_test_results(results$base, test_results, "docker_test_results.json")

  # Print final result
  cat("\n", strrep("=", 50), "\n", sep = "")
  if (test_results$success) {
    print_test_result(TRUE, sprintf("Docker installation tests completed in %.1fs", duration))
    return(TRUE)
  } else {
    print_test_result(FALSE, sprintf("Docker installation tests failed after %.1fs", duration))
    for (err in errors) {
      cat(sprintf("  Error: %s\n", err))
    }
    return(FALSE)
  }
}

# =============================================================================
# Run test
# =============================================================================

success <- run_docker_test()
quit(status = if(success) 0 else 1, save = "no")
