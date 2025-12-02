# CASSIA Test 14: Input Validation (R) - INSTALL MODE
# ====================================================
# Tests input validation for runCASSIA and runCASSIA_batch functions.
# Uses devtools::install_local() for full package installation testing.
#
# This test suite verifies that:
# 1. Invalid inputs are rejected with clear error messages
# 2. Valid inputs are accepted and processed correctly
# 3. Edge cases are handled appropriately (warnings vs errors)
#
# Usage:
#     Rscript test_input_validation_install.R

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

# Helper function to print individual test results
print_subtest_result <- function(test_name, passed, message = NULL) {
  status <- if (passed) "[OK]" else "[X]"
  log_msg(" ", status, test_name)
  if (!is.null(message) && !passed) {
    log_msg("      ", message)
  }
}

# Import Python validation modules
validation <- NULL
exceptions <- NULL

setup_validation_modules <- function() {
  # Import validation module
  validation <<- reticulate::import("CASSIA.core.validation")
  exceptions <<- reticulate::import("CASSIA.core.exceptions")
}

test_marker_list_validation <- function() {
  log_separator()
  log_msg("Testing Marker List Validation")
  log_separator()

  tests_passed <- 0
  tests_failed <- 0

  # Test 1: Valid list of strings
  tryCatch({
    result <- validation$validate_marker_list(c("CD4", "CD8", "FOXP3"))
    if (length(result) == 3 && all(result == c("CD4", "CD8", "FOXP3"))) {
      print_subtest_result("Valid list of strings", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Valid list of strings", FALSE, paste("Got:", paste(result, collapse=", ")))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Valid list of strings", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 2: Valid comma-separated string
  tryCatch({
    result <- validation$validate_marker_list("CD4, CD8, FOXP3")
    if (length(result) == 3 && all(result == c("CD4", "CD8", "FOXP3"))) {
      print_subtest_result("Valid comma-separated string", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Valid comma-separated string", FALSE, paste("Got:", paste(result, collapse=", ")))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Valid comma-separated string", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 3: None should raise error
  tryCatch({
    validation$validate_marker_list(NULL)
    print_subtest_result("None raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    if (grepl("cannot be None|None|null", e$message, ignore.case = TRUE)) {
      print_subtest_result("None raises error", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("None raises error", FALSE, paste("Wrong error:", e$message))
      tests_failed <- tests_failed + 1
    }
  })

  # Test 4: Empty list should raise error
  tryCatch({
    validation$validate_marker_list(character(0))
    print_subtest_result("Empty list raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    if (grepl("empty", e$message, ignore.case = TRUE)) {
      print_subtest_result("Empty list raises error", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Empty list raises error", FALSE, paste("Wrong error:", e$message))
      tests_failed <- tests_failed + 1
    }
  })

  # Test 5: Too many markers should raise error
  tryCatch({
    many_markers <- paste0("GENE", 1:600)
    validation$validate_marker_list(many_markers)
    print_subtest_result("Too many markers raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    if (grepl("Too many markers|500", e$message, ignore.case = TRUE)) {
      print_subtest_result("Too many markers raises error (>500)", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Too many markers raises error", FALSE, paste("Wrong error:", e$message))
      tests_failed <- tests_failed + 1
    }
  })

  # Test 6: Ensembl IDs should raise error
  tryCatch({
    validation$validate_marker_list(c("ENSG00000141510", "ENSG00000134644", "ENSG00000087086"))
    print_subtest_result("Ensembl IDs raise error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    if (grepl("Ensembl", e$message, ignore.case = TRUE)) {
      print_subtest_result("Ensembl IDs raise error", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Ensembl IDs raise error", FALSE, paste("Wrong error:", e$message))
      tests_failed <- tests_failed + 1
    }
  })

  return(c(tests_passed, tests_failed))
}

test_temperature_validation <- function() {
  log_separator()
  log_msg("Testing Temperature Validation")
  log_separator()

  tests_passed <- 0
  tests_failed <- 0

  # Test 1: Valid temperature (0)
  tryCatch({
    result <- validation$validate_temperature(0)
    if (result == 0.0) {
      print_subtest_result("Valid temperature (0)", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Valid temperature (0)", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Valid temperature (0)", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 2: Valid temperature (1.5)
  tryCatch({
    result <- validation$validate_temperature(1.5)
    if (result == 1.5) {
      print_subtest_result("Valid temperature (1.5)", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Valid temperature (1.5)", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Valid temperature (1.5)", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 3: Negative temperature should raise error
  tryCatch({
    validation$validate_temperature(-0.5)
    print_subtest_result("Negative raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    if (grepl(">= 0|negative|must be", e$message, ignore.case = TRUE)) {
      print_subtest_result("Negative raises error", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Negative raises error", FALSE, paste("Wrong error:", e$message))
      tests_failed <- tests_failed + 1
    }
  })

  # Test 4: None should raise error
  tryCatch({
    validation$validate_temperature(NULL)
    print_subtest_result("None raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    if (grepl("cannot be None|None|null", e$message, ignore.case = TRUE)) {
      print_subtest_result("None raises error", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("None raises error", FALSE, paste("Wrong error:", e$message))
      tests_failed <- tests_failed + 1
    }
  })

  # Test 5: High temperature (5.0) should be accepted
  tryCatch({
    result <- validation$validate_temperature(5.0)
    if (result == 5.0) {
      print_subtest_result("High temperature accepted (no upper bound)", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("High temperature accepted", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("High temperature accepted", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  return(c(tests_passed, tests_failed))
}

test_tissue_species_validation <- function() {
  log_separator()
  log_msg("Testing Tissue/Species Validation")
  log_separator()

  tests_passed <- 0
  tests_failed <- 0

  # Test 1: Valid tissue
  tryCatch({
    result <- validation$validate_tissue("lung")
    if (result == "lung") {
      print_subtest_result("Valid tissue", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Valid tissue", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Valid tissue", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 2: Special tissue value 'none'
  tryCatch({
    result <- validation$validate_tissue("none")
    if (result == "none") {
      print_subtest_result("Special tissue 'none'", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Special tissue 'none'", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Special tissue 'none'", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 3: Valid species
  tryCatch({
    result <- validation$validate_species("human")
    if (result == "human") {
      print_subtest_result("Valid species", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Valid species", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Valid species", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  return(c(tests_passed, tests_failed))
}

test_provider_validation <- function() {
  log_separator()
  log_msg("Testing Provider Validation")
  log_separator()

  tests_passed <- 0
  tests_failed <- 0

  # Test 1: Valid providers
  for (provider in c("openai", "anthropic", "openrouter")) {
    tryCatch({
      result <- validation$validate_provider(provider)
      if (result == provider) {
        print_subtest_result(paste0("Valid provider '", provider, "'"), TRUE)
        tests_passed <- tests_passed + 1
      } else {
        print_subtest_result(paste0("Valid provider '", provider, "'"), FALSE, paste("Got:", result))
        tests_failed <- tests_failed + 1
      }
    }, error = function(e) {
      print_subtest_result(paste0("Valid provider '", provider, "'"), FALSE, e$message)
      tests_failed <<- tests_failed + 1
    })
  }

  # Test 2: Case insensitive
  tryCatch({
    result <- validation$validate_provider("OpenAI")
    if (result == "openai") {
      print_subtest_result("Case insensitive provider", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Case insensitive provider", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Case insensitive provider", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 3: HTTP URL accepted
  tryCatch({
    url <- "http://localhost:8000/v1"
    result <- validation$validate_provider(url)
    if (result == url) {
      print_subtest_result("HTTP URL accepted", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("HTTP URL accepted", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("HTTP URL accepted", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 4: Unknown provider raises error
  tryCatch({
    validation$validate_provider("gemini")
    print_subtest_result("Unknown provider raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    if (grepl("Unknown provider", e$message, ignore.case = TRUE)) {
      print_subtest_result("Unknown provider raises error", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Unknown provider raises error", FALSE, paste("Wrong error:", e$message))
      tests_failed <- tests_failed + 1
    }
  })

  # Test 5: None raises error
  tryCatch({
    validation$validate_provider(NULL)
    print_subtest_result("None provider raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("None provider raises error", TRUE)
    tests_passed <- tests_passed + 1
  })

  return(c(tests_passed, tests_failed))
}

test_model_validation <- function() {
  log_separator()
  log_msg("Testing Model Validation")
  log_separator()

  tests_passed <- 0
  tests_failed <- 0

  # Test 1: Valid model names
  for (model in c("gpt-4", "claude-3", "gemini-pro")) {
    tryCatch({
      result <- validation$validate_model(model)
      if (result == model) {
        print_subtest_result(paste0("Valid model '", model, "'"), TRUE)
        tests_passed <- tests_passed + 1
      } else {
        print_subtest_result(paste0("Valid model '", model, "'"), FALSE, paste("Got:", result))
        tests_failed <- tests_failed + 1
      }
    }, error = function(e) {
      print_subtest_result(paste0("Valid model '", model, "'"), FALSE, e$message)
      tests_failed <<- tests_failed + 1
    })
  }

  # Test 2: None raises error
  tryCatch({
    validation$validate_model(NULL)
    print_subtest_result("None model raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("None model raises error", TRUE)
    tests_passed <- tests_passed + 1
  })

  # Test 3: Empty string raises error
  tryCatch({
    validation$validate_model("")
    print_subtest_result("Empty model raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("Empty model raises error", TRUE)
    tests_passed <- tests_passed + 1
  })

  return(c(tests_passed, tests_failed))
}

test_batch_parameters <- function() {
  log_separator()
  log_msg("Testing Batch Parameter Validation")
  log_separator()

  tests_passed <- 0
  tests_failed <- 0

  # Test 1: Valid positive integers
  tryCatch({
    result1 <- validation$validate_positive_int(as.integer(10), "n_genes")
    result2 <- validation$validate_positive_int(as.integer(5), "max_workers")
    if (result1 == 10 && result2 == 5) {
      print_subtest_result("Valid positive integers", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Valid positive integers", FALSE, paste("Got:", result1, result2))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Valid positive integers", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 2: Zero allowed with flag
  tryCatch({
    result <- validation$validate_positive_int(as.integer(0), "max_retries", allow_zero = TRUE)
    if (result == 0) {
      print_subtest_result("Zero allowed with allow_zero=TRUE", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("Zero allowed with allow_zero=TRUE", FALSE, paste("Got:", result))
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("Zero allowed with allow_zero=TRUE", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 3: Zero rejected without flag
  tryCatch({
    validation$validate_positive_int(as.integer(0), "n_genes", allow_zero = FALSE)
    print_subtest_result("Zero rejected without allow_zero", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("Zero rejected without allow_zero", TRUE)
    tests_passed <- tests_passed + 1
  })

  # Test 4: Valid ranking methods
  for (method in c("avg_log2FC", "p_val_adj", "pct_diff", "Score")) {
    tryCatch({
      result <- validation$validate_ranking_method(method)
      if (result == method) {
        print_subtest_result(paste0("Valid ranking method '", method, "'"), TRUE)
        tests_passed <- tests_passed + 1
      } else {
        print_subtest_result(paste0("Valid ranking method '", method, "'"), FALSE, paste("Got:", result))
        tests_failed <- tests_failed + 1
      }
    }, error = function(e) {
      print_subtest_result(paste0("Valid ranking method '", method, "'"), FALSE, e$message)
      tests_failed <<- tests_failed + 1
    })
  }

  # Test 5: Invalid ranking method raises error
  tryCatch({
    validation$validate_ranking_method("invalid_method")
    print_subtest_result("Invalid ranking method raises error", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("Invalid ranking method raises error", TRUE)
    tests_passed <- tests_passed + 1
  })

  return(c(tests_passed, tests_failed))
}

test_integrated_validation <- function() {
  log_separator()
  log_msg("Testing Integrated runCASSIA Validation")
  log_separator()

  tests_passed <- 0
  tests_failed <- 0

  # Test 1: All valid inputs
  tryCatch({
    markers <- c("CD4", "CD8", "FOXP3", "CD3", "CD19", "CD20", "CD56", "CD14", "CD68", "HLA-DR")
    result <- validation$validate_runCASSIA_inputs(
      model = "gpt-4",
      temperature = 0,
      marker_list = markers,
      tissue = "lung",
      species = "human",
      provider = "openai",
      additional_info = "Test analysis"
    )
    if ("marker_list" %in% names(result) && "temperature" %in% names(result)) {
      print_subtest_result("All valid inputs accepted", TRUE)
      tests_passed <- tests_passed + 1
    } else {
      print_subtest_result("All valid inputs accepted", FALSE, "Missing expected keys in result")
      tests_failed <- tests_failed + 1
    }
  }, error = function(e) {
    print_subtest_result("All valid inputs accepted", FALSE, e$message)
    tests_failed <<- tests_failed + 1
  })

  # Test 2: Invalid marker_list fails early
  tryCatch({
    validation$validate_runCASSIA_inputs(
      model = "gpt-4",
      temperature = 0,
      marker_list = NULL,  # Invalid
      tissue = "lung",
      species = "human",
      provider = "openai"
    )
    print_subtest_result("Invalid marker_list caught", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("Invalid marker_list caught", TRUE)
    tests_passed <- tests_passed + 1
  })

  # Test 3: Invalid temperature fails
  tryCatch({
    validation$validate_runCASSIA_inputs(
      model = "gpt-4",
      temperature = -1,  # Invalid
      marker_list = c("CD4", "CD8", "FOXP3", "CD3", "CD19", "CD20", "CD56", "CD14", "CD68", "HLA-DR"),
      tissue = "lung",
      species = "human",
      provider = "openai"
    )
    print_subtest_result("Invalid temperature caught", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("Invalid temperature caught", TRUE)
    tests_passed <- tests_passed + 1
  })

  # Test 4: Invalid provider fails
  tryCatch({
    validation$validate_runCASSIA_inputs(
      model = "gpt-4",
      temperature = 0,
      marker_list = c("CD4", "CD8", "FOXP3", "CD3", "CD19", "CD20", "CD56", "CD14", "CD68", "HLA-DR"),
      tissue = "lung",
      species = "human",
      provider = "invalid"  # Invalid
    )
    print_subtest_result("Invalid provider caught", FALSE, "Expected error but got none")
    tests_failed <- tests_failed + 1
  }, error = function(e) {
    print_subtest_result("Invalid provider caught", TRUE)
    tests_passed <- tests_passed + 1
  })

  return(c(tests_passed, tests_failed))
}

run_all_tests <- function() {
  print_test_header("14 - Input Validation (R) [INSTALL MODE]")

  # Load configuration
  config <- load_config()
  print_config_summary(config)

  # Setup CASSIA R package - INSTALL MODE (full package installation)
  tryCatch({
    setup_cassia_install(force = TRUE)
    message("CASSIA R package installed and loaded successfully")
  }, error = function(e) {
    message("Error installing CASSIA: ", e$message)
    return(FALSE)
  })

  # Create results directory
  results_dir <- create_results_dir("14_input_validation")

  start_logging(results_dir)
  log_msg("Results will be saved to:", results_dir)

  # Setup validation modules
  setup_validation_modules()

  start_time <- Sys.time()
  total_passed <- 0
  total_failed <- 0
  test_group_results <- list()

  # Run all test groups
  results <- test_marker_list_validation()
  total_passed <- total_passed + results[1]
  total_failed <- total_failed + results[2]
  test_group_results$marker_list <- list(passed = results[1], failed = results[2])

  results <- test_temperature_validation()
  total_passed <- total_passed + results[1]
  total_failed <- total_failed + results[2]
  test_group_results$temperature <- list(passed = results[1], failed = results[2])

  results <- test_tissue_species_validation()
  total_passed <- total_passed + results[1]
  total_failed <- total_failed + results[2]
  test_group_results$tissue_species <- list(passed = results[1], failed = results[2])

  results <- test_provider_validation()
  total_passed <- total_passed + results[1]
  total_failed <- total_failed + results[2]
  test_group_results$provider <- list(passed = results[1], failed = results[2])

  results <- test_model_validation()
  total_passed <- total_passed + results[1]
  total_failed <- total_failed + results[2]
  test_group_results$model <- list(passed = results[1], failed = results[2])

  results <- test_batch_parameters()
  total_passed <- total_passed + results[1]
  total_failed <- total_failed + results[2]
  test_group_results$batch_parameters <- list(passed = results[1], failed = results[2])

  results <- test_integrated_validation()
  total_passed <- total_passed + results[1]
  total_failed <- total_failed + results[2]
  test_group_results$integrated <- list(passed = results[1], failed = results[2])

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Print summary
  log_separator()
  log_msg("TEST SUMMARY")
  log_separator()
  log_msg("Tests passed:", total_passed)
  log_msg("Tests failed:", total_failed)
  log_msg("Total tests: ", total_passed + total_failed)

  # Determine status
  if (total_failed == 0) {
    log_msg("\n*** ALL TESTS PASSED ***")
    status <- "passed"
  } else {
    log_msg("\n***", total_failed, "TEST(S) FAILED ***")
    status <- "failed"
  }

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "input_validation_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(),
    errors = if (total_failed > 0) list(paste(total_failed, "validation tests failed")) else list()
  )
  metadata$test_summary <- list(
    passed = total_passed,
    failed = total_failed,
    total = total_passed + total_failed
  )
  save_test_metadata(results_dir, metadata)

  save_test_results(results_dir, list(
    test_groups = test_group_results,
    mode = "install"
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_all_tests()
quit(status = if (success) 0 else 1)
