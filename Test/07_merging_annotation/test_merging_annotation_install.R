# CASSIA Test 07: Merging Annotation (R) - INSTALL MODE
# =====================================================
# Tests the merge_annotations and merge_annotations_all functions via R package.
# Uses devtools::install_local() for full package installation testing.
#
# Usage:
#     Rscript test_merging_annotation_install.R
#
# Note: This test requires batch annotation results. If none exist, it will
# run a batch annotation first.

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

run_merging_annotation_install_test <- function() {
  print_test_header("07 - Merging Annotation (R) [INSTALL MODE]")

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
  results_dirs <- create_results_dir("07_merging_annotation", get_test_mode())

  start_logging(results_dirs$logs)
  log_msg("Results will be saved to:", results_dirs$base)

  # Check for existing batch results from Test 02
  batch_results_dir <- get_latest_results("02_batch_annotation")
  batch_results_file <- NULL

  if (!is.null(batch_results_dir)) {
    potential_file <- file.path(batch_results_dir, "batch_results_full.csv")
    if (file.exists(potential_file)) {
      batch_results_file <- potential_file
      log_msg("\nUsing existing batch results:", batch_results_file)
    }
  }

  # If no existing results, run batch annotation first
  if (is.null(batch_results_file)) {
    log_msg("\nNo existing batch results found. Running batch annotation first...")
    marker_df <- get_full_marker_dataframe()
    batch_output <- file.path(results_dirs$outputs, "batch_for_merge")

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
    batch_results_file <- paste0(batch_output, "_full.csv")
  }

  # Run merging tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  merge_results <- list()

  tryCatch({
    # Test 1: Single detail level merge (broad)
    log_msg("\n--- Test 1: Single detail level merge (broad) ---")
    broad_output <- file.path(results_dirs$outputs, "merge_broad.csv")

    result_broad <- CASSIA::merge_annotations(
      csv_path = batch_results_file,
      output_path = broad_output,
      provider = llm_config$provider %||% "openrouter",
      model = llm_config$model %||% "google/gemini-2.5-flash",
      detail_level = "broad",
      batch_size = 10
    )

    if (!is.null(result_broad) && "Merged_Grouping_1" %in% names(result_broad)) {
      log_msg("  Broad merge: SUCCESS")
      log_msg("  Output rows:", nrow(result_broad))
      log_msg("  Sample groupings:")
      cluster_col <- if ("Cluster ID" %in% names(result_broad)) "Cluster ID" else "True Cell Type"
      for (i in 1:min(3, nrow(result_broad))) {
        log_msg("   ", result_broad[[cluster_col]][i], "->", result_broad$Merged_Grouping_1[i])
      }
      merge_results$broad <- list(status = "success", output_file = broad_output, num_rows = nrow(result_broad))
    } else {
      errors <- c(errors, "Broad merge failed - missing Merged_Grouping_1 column")
      merge_results$broad <- list(status = "failed")
    }

    # Test 2: Single detail level merge (detailed)
    log_msg("\n--- Test 2: Single detail level merge (detailed) ---")
    detailed_output <- file.path(results_dirs$outputs, "merge_detailed.csv")

    result_detailed <- CASSIA::merge_annotations(
      csv_path = batch_results_file,
      output_path = detailed_output,
      provider = llm_config$provider %||% "openrouter",
      model = llm_config$model %||% "google/gemini-2.5-flash",
      detail_level = "detailed",
      batch_size = 10
    )

    if (!is.null(result_detailed) && "Merged_Grouping_2" %in% names(result_detailed)) {
      log_msg("  Detailed merge: SUCCESS")
      log_msg("  Output rows:", nrow(result_detailed))
      log_msg("  Sample groupings:")
      cluster_col <- if ("Cluster ID" %in% names(result_detailed)) "Cluster ID" else "True Cell Type"
      for (i in 1:min(3, nrow(result_detailed))) {
        log_msg("   ", result_detailed[[cluster_col]][i], "->", result_detailed$Merged_Grouping_2[i])
      }
      merge_results$detailed <- list(status = "success", output_file = detailed_output, num_rows = nrow(result_detailed))
    } else {
      errors <- c(errors, "Detailed merge failed - missing Merged_Grouping_2 column")
      merge_results$detailed <- list(status = "failed")
    }

    # Test 3: All detail levels merge (parallel)
    log_msg("\n--- Test 3: All detail levels merge (parallel) ---")
    all_output <- file.path(results_dirs$outputs, "merge_all.csv")

    result_all <- CASSIA::merge_annotations_all(
      csv_path = batch_results_file,
      output_path = all_output,
      provider = llm_config$provider %||% "openrouter",
      model = llm_config$model %||% "google/gemini-2.5-flash",
      batch_size = 10
    )

    expected_columns <- c("Merged_Grouping_1", "Merged_Grouping_2", "Merged_Grouping_3")
    if (!is.null(result_all) && all(expected_columns %in% names(result_all))) {
      log_msg("  All levels merge: SUCCESS")
      log_msg("  Output rows:", nrow(result_all))
      log_msg("  Columns:", paste(expected_columns, collapse = ", "))
      log_msg("\n  Sample comparison:")
      cluster_col <- if ("Cluster ID" %in% names(result_all)) "Cluster ID" else "True Cell Type"
      for (i in 1:min(3, nrow(result_all))) {
        log_msg("   ", result_all[[cluster_col]][i], ":")
        log_msg("      Broad:        ", result_all$Merged_Grouping_1[i])
        log_msg("      Detailed:     ", result_all$Merged_Grouping_2[i])
        log_msg("      Very Detailed:", result_all$Merged_Grouping_3[i])
      }
      merge_results$all <- list(status = "success", output_file = all_output, num_rows = nrow(result_all), columns = expected_columns)
    } else {
      missing <- expected_columns[!expected_columns %in% names(result_all)]
      errors <- c(errors, paste("All levels merge failed - missing columns:", paste(missing, collapse = ", ")))
      merge_results$all <- list(status = "failed")
    }

    # Determine overall status
    statuses <- sapply(merge_results, function(r) r$status)
    if (all(statuses == "success")) {
      status <- "passed"
    } else if (any(statuses == "success")) {
      status <- "partial"
    } else {
      status <- "failed"
    }

  }, error = function(e) {
    errors <<- c(errors, e$message)
    status <<- "error"
    log_msg("\nError:", e$message)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "merging_annotation_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list("all"),
    errors = errors
  )
  save_test_metadata(results_dirs$outputs, metadata)

  save_test_results(results_dirs$outputs, list(
    batch_results_file = batch_results_file,
    merge_results = merge_results,
    detail_levels_tested = c("broad", "detailed", "very_detailed", "all"),
    mode = "install"
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_merging_annotation_install_test()
quit(status = if (success) 0 else 1)
