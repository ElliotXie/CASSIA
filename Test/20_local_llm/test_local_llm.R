# CASSIA Test 20: Local LLM Support (R)
# =====================================
# Tests that CASSIA correctly handles local LLM endpoints
# without requiring an API key.
#
# Includes a full batch annotation test using real marker data.
#
# Usage:
#     Rscript test_local_llm.R

# Get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg), winslash = "/")))
  }
  return(normalizePath(getwd(), winslash = "/"))
}
script_dir <- get_script_dir()

# Source shared utilities
source(file.path(script_dir, "..", "shared", "r", "test_utils.R"))
source(file.path(script_dir, "..", "shared", "r", "fixtures.R"))
source(file.path(script_dir, "..", "shared", "r", "result_manager.R"))
source(file.path(script_dir, "..", "shared", "r", "logging_manager.R"))


check_ollama_running <- function(host = "localhost", port = 11434, timeout = 2) {
  tryCatch({
    con <- socketConnection(host = host, port = port,
                           open = "r", blocking = TRUE, timeout = timeout)
    close(con)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}


get_ollama_model <- function() {
  # Get the first available Ollama model
  tryCatch({
    response <- httr::GET("http://localhost:11434/api/tags", httr::timeout(5))
    if (httr::status_code(response) == 200) {
      data <- httr::content(response, as = "parsed")
      models <- sapply(data$models, function(m) m$name)
      if (length(models) > 0) {
        return(models[1])
      }
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}


test_localhost_detection <- function() {
  # Test localhost detection in R
  test_cases <- list(
    list(url = "http://localhost:11434/v1", expected = TRUE),
    list(url = "http://localhost:1234/v1", expected = TRUE),
    list(url = "http://127.0.0.1:8000/v1", expected = TRUE),
    list(url = "https://api.deepseek.com", expected = FALSE)
  )

  errors <- list()
  for (tc in test_cases) {
    is_localhost <- grepl("localhost|127\\.0\\.0\\.1", tc$url, ignore.case = TRUE)
    if (is_localhost != tc$expected) {
      errors <- c(errors, paste0("URL '", tc$url, "': expected ", tc$expected, ", got ", is_localhost))
    }
  }

  list(success = length(errors) == 0, errors = errors)
}


test_api_key_bypass <- function() {
  # Test that setLLMApiKey works without key for localhost
  tryCatch({
    # Clear any existing key
    Sys.unsetenv("CUSTOMIZED_API_KEY")

    # This should NOT error for localhost
    setLLMApiKey(provider = "http://localhost:11434/v1")

    # Check that placeholder was set
    key <- Sys.getenv("CUSTOMIZED_API_KEY")
    if (key == "ollama") {
      return(list(success = TRUE, errors = list()))
    } else {
      return(list(success = FALSE, errors = list(paste("Expected 'ollama' placeholder, got:", key))))
    }
  }, error = function(e) {
    list(success = FALSE, errors = list(e$message))
  })
}


test_remote_requires_key <- function() {
  # Test that remote URLs still require API key
  tryCatch({
    # Clear any existing key
    Sys.unsetenv("CUSTOMIZED_API_KEY")

    # This should error
    setLLMApiKey(provider = "https://api.deepseek.com")

    # If we get here, it didn't error - that's wrong
    list(success = FALSE, errors = list("Remote URL should require API key"))
  }, error = function(e) {
    if (grepl("API key", e$message, ignore.case = TRUE)) {
      list(success = TRUE, errors = list())
    } else {
      list(success = FALSE, errors = list(paste("Wrong error:", e$message)))
    }
  })
}


test_ollama_connection <- function() {
  # Test actual Ollama connection if running
  # Note: This test verifies Ollama is reachable via HTTP API
  # The actual LLM call is tested in Test 5 (batch annotation)
  if (!check_ollama_running()) {
    return(list(status = "skipped", errors = list("Ollama not running on localhost:11434")))
  }

  model_to_use <- get_ollama_model()
  if (is.null(model_to_use)) {
    return(list(status = "skipped", errors = list("No Ollama models available")))
  }

  message("  Using Ollama model: ", model_to_use)
  message("  Ollama API is reachable and has models available")
  message("  (Full LLM call tested in Test 5 - Batch Annotation)")

  # Just verify we can set up the API key for localhost
  tryCatch({
    setLLMApiKey(provider = "http://localhost:11434/v1")
    return(list(status = TRUE, errors = list()))
  }, error = function(e) {
    return(list(status = FALSE, errors = list(paste("Failed to set up localhost:", e$message))))
  })
}


test_batch_annotation_ollama <- function(output_dir) {
  # Test full batch annotation with local Ollama - like test 02
  if (!check_ollama_running()) {
    return(list(status = "skipped", errors = list("Ollama not running on localhost:11434")))
  }

  model_to_use <- get_ollama_model()
  if (is.null(model_to_use)) {
    return(list(status = "skipped", errors = list("No Ollama models available")))
  }

  message("  Using Ollama model: ", model_to_use)
  message("  Provider: http://localhost:11434/v1")

  # Load marker data - use 2 clusters like test 02
  test_clusters <- c("monocyte", "plasma cell")
  full_df <- get_full_marker_dataframe()
  marker_df <- full_df[full_df$Broad.cell.type %in% test_clusters, ]

  message("  Testing annotation for: ", paste(test_clusters, collapse = ", "))
  message("  Marker data shape: ", nrow(marker_df), " rows x ", ncol(marker_df), " cols")

  output_name <- file.path(output_dir, "ollama_batch_results")

  tryCatch({
    message("  Running runCASSIA_batch with local Ollama...")

    # Set up API key for localhost
    setLLMApiKey(provider = "http://localhost:11434/v1")

    CASSIA::runCASSIA_batch(
      marker = marker_df,
      output_name = output_name,
      n_genes = 30,
      model = model_to_use,
      temperature = 0.3,
      tissue = "large intestine",
      species = "human",
      max_workers = 1,  # Single worker for local LLM
      provider = "http://localhost:11434/v1",
      validator_involvement = "v1"
    )

    # Check output files
    summary_csv <- paste0(output_name, "_summary.csv")
    if (file.exists(summary_csv)) {
      results_df <- read.csv(summary_csv)
      clusters_annotated <- nrow(results_df)
      message("  Clusters annotated: ", clusters_annotated, "/", length(test_clusters))

      if (clusters_annotated == length(test_clusters)) {
        # Check if we got actual annotations
        if ("Predicted.General.Cell.Type" %in% names(results_df)) {
          cell_type <- results_df$Predicted.General.Cell.Type[1]
          message("  First annotation result: ", cell_type)
        }
        return(list(status = TRUE, errors = list()))
      } else {
        return(list(status = FALSE, errors = list(paste("Only", clusters_annotated, "/", length(test_clusters), "clusters annotated"))))
      }
    } else {
      return(list(status = FALSE, errors = list("Summary CSV not created")))
    }
  }, error = function(e) {
    error_msg <- e$message
    # Handle known issues gracefully
    if (grepl("timeout", error_msg, ignore.case = TRUE)) {
      return(list(status = "skipped", errors = list(paste("Ollama timeout:", substr(error_msg, 1, 100)))))
    }
    return(list(status = FALSE, errors = list(paste("Batch annotation failed:", substr(error_msg, 1, 200)))))
  })
}


run_local_llm_test <- function() {
  print_test_header("20 - Local LLM Support Testing (R)")

  config <- load_config()
  print_config_summary(config)

  # Setup CASSIA
  tryCatch({
    setup_cassia()
    message("CASSIA R package loaded successfully")
  }, error = function(e) {
    message("Error loading CASSIA: ", e$message)
    return(FALSE)
  })

  # Create results directory
  results <- create_results_dir("20_local_llm", get_test_mode())
  start_logging(results$logs)

  start_time <- Sys.time()
  test_results <- list(tests = list())
  passed <- 0
  failed <- 0
  skipped <- 0
  errors <- list()

  # Test 1: Localhost detection
  message("\n", paste(rep("=", 50), collapse = ""))
  message("Test 1: Localhost Detection Logic")
  message(paste(rep("=", 50), collapse = ""))
  result <- test_localhost_detection()
  test_results$tests <- c(test_results$tests, list(list(
    name = "localhost_detection",
    status = if (result$success) "passed" else "failed",
    errors = result$errors
  )))
  if (result$success) {
    message("  [OK] Localhost detection working correctly")
    passed <- passed + 1
  } else {
    message("  [X] FAILED: ", paste(result$errors, collapse = ", "))
    failed <- failed + 1
    errors <- c(errors, result$errors)
  }

  # Test 2: API key bypass
  message("\n", paste(rep("=", 50), collapse = ""))
  message("Test 2: API Key Bypass for Localhost")
  message(paste(rep("=", 50), collapse = ""))
  result <- test_api_key_bypass()
  test_results$tests <- c(test_results$tests, list(list(
    name = "api_key_bypass",
    status = if (result$success) "passed" else "failed",
    errors = result$errors
  )))
  if (result$success) {
    message("  [OK] API key not required for localhost")
    passed <- passed + 1
  } else {
    message("  [X] FAILED: ", paste(result$errors, collapse = ", "))
    failed <- failed + 1
    errors <- c(errors, result$errors)
  }

  # Test 3: Remote requires key
  message("\n", paste(rep("=", 50), collapse = ""))
  message("Test 3: Remote URL Requires API Key")
  message(paste(rep("=", 50), collapse = ""))
  result <- test_remote_requires_key()
  test_results$tests <- c(test_results$tests, list(list(
    name = "remote_requires_key",
    status = if (result$success) "passed" else "failed",
    errors = result$errors
  )))
  if (result$success) {
    message("  [OK] Remote URLs correctly require API key")
    passed <- passed + 1
  } else {
    message("  [X] FAILED: ", paste(result$errors, collapse = ", "))
    failed <- failed + 1
    errors <- c(errors, result$errors)
  }

  # Test 4: Ollama connection (optional)
  message("\n", paste(rep("=", 50), collapse = ""))
  message("Test 4: Ollama Connection (Optional)")
  message(paste(rep("=", 50), collapse = ""))
  result <- test_ollama_connection()
  if (identical(result$status, "skipped")) {
    message("  [SKIP] ", result$errors[[1]])
    skipped <- skipped + 1
    test_results$tests <- c(test_results$tests, list(list(
      name = "ollama_connection",
      status = "skipped",
      errors = result$errors
    )))
  } else if (isTRUE(result$status)) {
    message("  [OK] Successfully connected to Ollama")
    passed <- passed + 1
    test_results$tests <- c(test_results$tests, list(list(
      name = "ollama_connection",
      status = "passed",
      errors = list()
    )))
  } else {
    message("  [X] FAILED: ", paste(result$errors, collapse = ", "))
    failed <- failed + 1
    errors <- c(errors, result$errors)
    test_results$tests <- c(test_results$tests, list(list(
      name = "ollama_connection",
      status = "failed",
      errors = result$errors
    )))
  }

  # Test 5: Full batch annotation with Ollama (like test 02)
  message("\n", paste(rep("=", 50), collapse = ""))
  message("Test 5: Batch Annotation with Local Ollama")
  message(paste(rep("=", 50), collapse = ""))
  result <- test_batch_annotation_ollama(results$outputs)
  if (identical(result$status, "skipped")) {
    message("  [SKIP] ", result$errors[[1]])
    skipped <- skipped + 1
    test_results$tests <- c(test_results$tests, list(list(
      name = "batch_annotation_ollama",
      status = "skipped",
      errors = result$errors
    )))
  } else if (isTRUE(result$status)) {
    message("  [OK] Batch annotation completed successfully")
    passed <- passed + 1
    test_results$tests <- c(test_results$tests, list(list(
      name = "batch_annotation_ollama",
      status = "passed",
      errors = list()
    )))
  } else {
    message("  [X] FAILED: ", paste(result$errors, collapse = ", "))
    failed <- failed + 1
    errors <- c(errors, result$errors)
    test_results$tests <- c(test_results$tests, list(list(
      name = "batch_annotation_ollama",
      status = "failed",
      errors = result$errors
    )))
  }

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Summary
  message("\n", paste(rep("=", 50), collapse = ""))
  message("SUMMARY")
  message(paste(rep("=", 50), collapse = ""))
  message("  Passed:  ", passed)
  message("  Failed:  ", failed)
  message("  Skipped: ", skipped)
  message("  Duration: ", round(duration, 2), "s")

  # Determine status
  if (failed > 0) {
    status <- "failed"
  } else if (passed == 0) {
    status <- "skipped"
  } else {
    status <- "passed"
  }

  # Save metadata
  metadata <- create_test_metadata(
    test_name = "local_llm_test",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(),
    errors = errors
  )
  metadata$test_summary <- list(passed = passed, failed = failed, skipped = skipped)
  save_test_metadata(results$outputs, metadata)

  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  stop_logging()
  return(success)
}

# Run test
success <- run_local_llm_test()
quit(status = if (success) 0 else 1)
