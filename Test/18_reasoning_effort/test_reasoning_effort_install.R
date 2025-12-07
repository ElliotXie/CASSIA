# CASSIA Test 18: Reasoning Effort (R INSTALLED MODE)
# ====================================================
# Tests the reasoning parameter with runCASSIA for:
# 1. Direct OpenAI - GPT-4o (no reasoning) - Chat Completions API
# 2. Direct OpenAI - GPT-5.1 (no reasoning) - Chat Completions API
# 3. Direct OpenAI - GPT-5 (with reasoning) - Responses API [SKIPPED - org verification]
# 4. Batch Annotation - GPT-5.1 (no reasoning) - CSV/HTML outputs
# 5. OpenRouter GPT-5.1 - reasoning effort LOW
#
# Usage:
#     Rscript test_reasoning_effort_install.R

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

# Test helper functions
test_single_annotation <- function(model, provider, reasoning = NULL, test_name = "test") {
  log_msg("\n", paste(rep("=", 60), collapse = ""))
  log_msg("TEST:", test_name)
  log_msg("  Model:", model)
  log_msg("  Provider:", provider)
  log_msg("  Reasoning:", if (is.null(reasoning)) "None" else reasoning)
  log_msg(paste(rep("=", 60), collapse = ""))

  # Get monocyte markers for quick test
  marker_list <- get_cluster_markers("monocyte", n_genes = 15)
  log_msg("Markers:", length(marker_list), "genes")

  start_time <- Sys.time()
  errors <- list()
  status <- "error"

  tryCatch({
    log_msg("Running runCASSIA...")

    result <- CASSIA::runCASSIA(
      model = model,
      temperature = 0.3,
      marker_list = marker_list,
      tissue = "blood",
      species = "human",
      provider = provider,
      validator_involvement = "v1",
      reasoning = reasoning
    )

    if (!is.null(result) && !is.null(result$main_cell_type)) {
      log_msg("Results:")
      log_msg("  Main cell type:", result$main_cell_type)
      log_msg("  Iterations:", result$iterations %||% "N/A")
      status <- "passed"
    } else {
      status <- "failed"
      errors <- list("No valid result returned")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    log_error(e)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  log_msg("Duration:", round(duration, 2), "s")
  log_msg("Status:", toupper(status))

  return(list(passed = status == "passed", duration = duration, errors = errors))
}

test_batch_annotation <- function(output_dir, model, provider, reasoning = NULL, test_name = "batch_test") {
  log_msg("\n", paste(rep("=", 60), collapse = ""))
  log_msg("TEST:", test_name)
  log_msg("  Model:", model)
  log_msg("  Provider:", provider)
  log_msg("  Reasoning:", if (is.null(reasoning)) "None" else reasoning)
  log_msg(paste(rep("=", 60), collapse = ""))

  # Load marker data for 2 clusters
  full_df <- load_markers()
  test_clusters <- c("monocyte", "plasma cell")
  marker_df <- full_df[full_df$Broad.cell.type %in% test_clusters, ]

  log_msg("Clusters:", paste(test_clusters, collapse = ", "))

  suffix <- if (is.null(reasoning)) "no_reasoning" else reasoning
  output_name <- file.path(output_dir, paste0("batch_", suffix))

  start_time <- Sys.time()
  errors <- list()
  status <- "error"

  tryCatch({
    log_msg("Running runCASSIA_batch...")

    CASSIA::runCASSIA_batch(
      marker = marker_df,
      output_name = output_name,
      n_genes = 15,
      model = model,
      temperature = 0.3,
      tissue = "blood",
      species = "human",
      max_workers = 2,
      provider = provider,
      validator_involvement = "v1",
      reasoning = reasoning
    )

    # Check output files
    full_csv <- paste0(output_name, "_full.csv")

    if (file.exists(full_csv)) {
      results_df <- read.csv(full_csv)
      clusters_annotated <- nrow(results_df)

      log_msg("Batch Results:")
      log_msg("  Clusters annotated:", clusters_annotated, "/", length(test_clusters))
      log_msg("  Output file:", basename(full_csv))

      if (clusters_annotated == length(test_clusters)) {
        status <- "passed"
      } else {
        status <- "failed"
        errors <- list(paste("Only", clusters_annotated, "/", length(test_clusters), "clusters annotated"))
      }
    } else {
      status <- "failed"
      errors <- list("Output CSV file not created")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    log_error(e)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  log_msg("Duration:", round(duration, 2), "s")
  log_msg("Status:", toupper(status))

  return(list(passed = status == "passed", duration = duration, errors = errors))
}

run_reasoning_effort_test <- function() {
  print_test_header("18 - Reasoning Effort (R INSTALLED MODE)")

  # Load configuration
  config <- load_config()

  # Setup API keys
  setup_api_keys()

  # Setup CASSIA R package (INSTALLED MODE)
  tryCatch({
    setup_cassia_install(force = TRUE)
    message("CASSIA R package loaded successfully (installed mode)")

    # Get package info
    pkg_version <- as.character(packageVersion("CASSIA"))
    message("CASSIA version: ", pkg_version)
  }, error = function(e) {
    message("Error loading CASSIA: ", e$message)
    return(FALSE)
  })

  # Create results directory (installed mode)
  results <- create_results_dir("18_reasoning_effort", "installed")
  start_logging(results$logs)

  total_start <- Sys.time()
  test_results <- list()
  all_errors <- list()

  # Test 1: GPT-4o without reasoning (Chat Completions API)
  log_msg("\n", paste(rep("#", 70), collapse = ""))
  log_msg("# RUNNING TEST 1: Direct OpenAI GPT-4o (Chat Completions API)")
  log_msg(paste(rep("#", 70), collapse = ""))
  result1 <- test_single_annotation("gpt-4o", "openai", NULL, "GPT-4o (no reasoning)")
  test_results$gpt4o <- result1
  all_errors <- c(all_errors, result1$errors)

  # Test 2: GPT-5.1 without reasoning (Chat Completions API)
  log_msg("\n", paste(rep("#", 70), collapse = ""))
  log_msg("# RUNNING TEST 2: Direct OpenAI GPT-5.1 (Chat Completions API)")
  log_msg(paste(rep("#", 70), collapse = ""))
  result2 <- test_single_annotation("gpt-5.1", "openai", NULL, "GPT-5.1 (no reasoning)")
  test_results$gpt5_no_reasoning <- result2
  all_errors <- c(all_errors, result2$errors)

  # Test 3: GPT-5 with reasoning (Responses API) - SKIPPED
  log_msg("\n", paste(rep("#", 70), collapse = ""))
  log_msg("# SKIPPING TEST 3: Direct OpenAI GPT-5 (Responses API)")
  log_msg("# Reason: Requires organization verification")
  log_msg(paste(rep("#", 70), collapse = ""))
  result3 <- list(passed = TRUE, duration = 0, errors = list())

  # Test 4: Batch annotation with GPT-5.1 (no reasoning)
  log_msg("\n", paste(rep("#", 70), collapse = ""))
  log_msg("# RUNNING TEST 4: Batch Annotation GPT-5.1 (Chat Completions API)")
  log_msg(paste(rep("#", 70), collapse = ""))
  result4 <- test_batch_annotation(results$outputs, "gpt-5.1", "openai", NULL, "Batch GPT-5.1 (no reasoning)")
  test_results$batch_gpt5 <- result4
  all_errors <- c(all_errors, result4$errors)

  # Test 5: OpenRouter GPT-5.1 with reasoning effort LOW
  log_msg("\n", paste(rep("#", 70), collapse = ""))
  log_msg("# RUNNING TEST 5: OpenRouter GPT-5.1 (Reasoning: LOW)")
  log_msg(paste(rep("#", 70), collapse = ""))
  result5 <- test_batch_annotation(results$outputs, "openai/gpt-5.1", "openrouter", "low", "OpenRouter (effort=LOW)")
  test_results$openrouter_low <- result5
  all_errors <- c(all_errors, result5$errors)

  total_duration <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))

  # Summary
  log_msg("\n", paste(rep("=", 70), collapse = ""))
  log_msg("SUMMARY (INSTALLED MODE)")
  log_msg(paste(rep("=", 70), collapse = ""))
  log_msg("Test 1 (GPT-4o, Chat Completions):", if (result1$passed) "PASSED" else "FAILED", "(", round(result1$duration, 2), "s)")
  log_msg("Test 2 (GPT-5.1, Chat Completions):", if (result2$passed) "PASSED" else "FAILED", "(", round(result2$duration, 2), "s)")
  log_msg("Test 3 (GPT-5, Responses API): SKIPPED (org verification)")
  log_msg("Test 4 (Batch GPT-5.1, no reasoning):", if (result4$passed) "PASSED" else "FAILED", "(", round(result4$duration, 2), "s)")
  log_msg("Test 5 (OpenRouter, effort=LOW):", if (result5$passed) "PASSED" else "FAILED", "(", round(result5$duration, 2), "s)")
  log_msg("\nTotal Duration:", round(total_duration, 2), "s")

  # Overall status
  all_passed <- result1$passed && result2$passed && result4$passed && result5$passed
  status <- if (all_passed) "passed" else "failed"

  # Save metadata with install info
  metadata <- create_test_metadata(
    test_name = "reasoning_effort_r_install",
    config = list(
      install_mode = TRUE,
      cassia_version = tryCatch(as.character(packageVersion("CASSIA")), error = function(e) "unknown"),
      tests = list(
        gpt4o = list(model = "gpt-4o", reasoning = NULL, api = "Chat Completions"),
        gpt5_no_reasoning = list(model = "gpt-5.1", reasoning = NULL, api = "Chat Completions"),
        gpt5_with_reasoning = list(model = "gpt-5", reasoning = "medium", api = "Responses", status = "skipped"),
        batch_gpt5 = list(model = "gpt-5.1", reasoning = NULL, api = "Chat Completions", type = "batch"),
        openrouter_low = list(model = "openai/gpt-5.1", reasoning = "low", provider = "openrouter")
      )
    ),
    duration_seconds = total_duration,
    status = status,
    errors = all_errors
  )
  save_test_metadata(results$outputs, metadata)

  # Print final result
  print_test_result(all_passed, paste("Total Duration:", round(total_duration, 2), "s"))

  stop_logging()
  return(all_passed)
}

# Run test
success <- run_reasoning_effort_test()
quit(status = if (success) 0 else 1)
