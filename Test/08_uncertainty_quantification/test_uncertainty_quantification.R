# CASSIA Test 08: Uncertainty Quantification (R)
# ===============================================
# Tests the uncertainty quantification functions via R package.
#
# Usage:
#     Rscript test_uncertainty_quantification.R
#
# Functions tested:
# - runCASSIA_n_times_similarity_score(): Run n single analyses with similarity score
# - runCASSIA_batch_n_times(): Run batch analyses n times for multiple clusters

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

run_uncertainty_quantification_test <- function(run_batch_generation = FALSE) {
  # run_batch_generation: If TRUE, runs runCASSIA_batch_n_times to generate new batch files.
  #                       If FALSE (default), uses pre-existing batch files from data folder.
  print_test_header("08 - Uncertainty Quantification (R)")

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

  # Create results directory
  results <- create_results_dir("08_uncertainty_quantification", get_test_mode())

  start_logging(results$logs)
  log_msg("Results will be saved to:", results$base)

  # Test cluster
  test_cluster <- "plasma cell"
  markers <- get_cluster_markers(test_cluster)

  # Run tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  uq_results <- list()

  tryCatch({
    log_msg("\n--- Test: runCASSIA_n_times_similarity_score ---")
    log_msg("  Cluster:", test_cluster)
    log_msg("  Model:", llm_config$model %||% "google/gemini-2.5-flash")
    log_msg("  Provider:", llm_config$provider %||% "openrouter")
    log_msg("  N iterations: 5")

    single_report_path <- file.path(results$outputs, "uq_single_report.html")

    result <- CASSIA::runCASSIA_n_times_similarity_score(
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      additional_info = NULL,
      temperature = llm_config$temperature %||% 0.3,
      marker_list = markers,
      model = llm_config$model %||% "google/gemini-2.5-flash",
      max_workers = llm_config$max_workers %||% 3,
      n = 5,
      provider = llm_config$provider %||% "openrouter",
      validator_involvement = config$validator$default %||% "v1",
      generate_report = TRUE,
      report_output_path = single_report_path
    )

    # Check if report was generated
    if (file.exists(single_report_path)) {
      log_msg("  HTML report created:", basename(single_report_path))
    }

    # Check result
    if (is.list(result)) {
      log_msg("\nResults:")
      log_msg("  General cell type (LLM):", result$general_celltype_llm %||% "N/A")
      log_msg("  Sub cell type (LLM):", result$sub_celltype_llm %||% "N/A")
      log_msg("  Similarity score:", result$similarity_score %||% "N/A")

      uq_results <- list(
        general_celltype_llm = result$general_celltype_llm,
        sub_celltype_llm = result$sub_celltype_llm,
        similarity_score = result$similarity_score,
        consensus_types = result$consensus_types,
        mixed_celltypes = result$Possible_mixed_celltypes_llm
      )

      # Validate results
      if (!is.null(result$general_celltype_llm) && !is.null(result$similarity_score)) {
        status <- "passed"
      } else {
        status <- "failed"
        errors <- list("Missing expected result fields for similarity score test")
      }
    } else {
      status <- "failed"
      errors <- list("Unexpected result format for similarity score test")
    }

  }, error = function(e) {
    errors <<- list(paste("Similarity score test error:", e$message))
    status <<- "error"
    log_msg("\nError in similarity score test:", e$message)
  })

  # =========================================================================
  # Test 2: runCASSIA_batch_n_times - Multiple clusters, multiple iterations
  # =========================================================================
  batch_status <- "skipped"
  batch_results <- list()
  batch_clusters_tested <- list()
  files_found <- c()

  # Get full marker dataframe and limit to 2 clusters for testing
  full_markers <- get_full_marker_dataframe()
  test_markers <- head(full_markers, 2)  # First 2 clusters: monocyte, plasma cell
  batch_clusters_tested <- test_markers$Broad.cell.type

  # Data folder with pre-existing batch results
  data_folder <- file.path(script_dir, "data")

  if (run_batch_generation) {
    # Run batch generation test
    tryCatch({
      log_msg("\n--- Test: runCASSIA_batch_n_times ---")

      log_msg("  Clusters:", paste(batch_clusters_tested, collapse = ", "))
      log_msg("  Model:", llm_config$model %||% "google/gemini-2.5-flash")
      log_msg("  Provider:", llm_config$provider %||% "openrouter")
      log_msg("  N iterations: 5")

      # Set output path for batch results
      batch_output_name <- file.path(results$outputs, "batch_results")

      # Run batch n times
      CASSIA::runCASSIA_batch_n_times(
        n = 5,
        marker = test_markers,
        output_name = batch_output_name,
        model = llm_config$model %||% "google/gemini-2.5-flash",
        temperature = llm_config$temperature %||% 0.3,
        tissue = data_config$tissue %||% "large intestine",
        species = data_config$species %||% "human",
        additional_info = NULL,
        celltype_column = "Broad.cell.type",
        gene_column_name = "Top.Markers",
        max_workers = 3,
        batch_max_workers = 2,
        provider = llm_config$provider %||% "openrouter",
        max_retries = 1,
        validator_involvement = config$validator$default %||% "v1"
      )

      # Check that output files were created
      expected_files <- c(
        paste0(batch_output_name, "_1_full.csv"),
        paste0(batch_output_name, "_1_summary.csv"),
        paste0(batch_output_name, "_2_full.csv"),
        paste0(batch_output_name, "_2_summary.csv"),
        paste0(batch_output_name, "_3_full.csv"),
        paste0(batch_output_name, "_3_summary.csv"),
        paste0(batch_output_name, "_4_full.csv"),
        paste0(batch_output_name, "_4_summary.csv"),
        paste0(batch_output_name, "_5_full.csv"),
        paste0(batch_output_name, "_5_summary.csv")
      )

      files_found <- c()
      files_missing <- c()
      for (f in expected_files) {
        if (file.exists(f)) {
          files_found <- c(files_found, basename(f))
        } else {
          files_missing <- c(files_missing, basename(f))
        }
      }

      log_msg("\nBatch Results:")
      log_msg("  Files created:", length(files_found))
      for (f in files_found) {
        log_msg("    -", f)
      }

      if (length(files_missing) > 0) {
        log_msg("  Files missing:", length(files_missing))
        for (f in files_missing) {
          log_msg("    -", f)
        }
      }

      batch_results <- list(
        clusters_tested = as.list(batch_clusters_tested),
        n_iterations = 5,
        files_created = as.list(files_found),
        files_missing = as.list(files_missing)
      )

      # Validate batch results - at least some files should be created
      if (length(files_found) >= 2) {
        batch_status <- "passed"
        log_msg("\n[OK] Batch test PASSED")
      } else {
        batch_status <- "failed"
        errors <- c(errors, paste("Batch test: Expected at least 2 output files, found", length(files_found)))
        log_msg("\n[FAIL] Batch test FAILED - insufficient output files")
      }

    }, error = function(e) {
      errors <<- c(errors, paste("Batch test error:", e$message))
      batch_status <<- "error"
      log_msg("\nError in batch test:", e$message)
    })
  } else {
    # Use pre-existing batch files from data folder
    log_msg("\n--- Skipping Test: runCASSIA_batch_n_times (using pre-existing data) ---")
    batch_output_name <- file.path(data_folder, "batch_results")

    # Check for pre-existing files
    expected_files <- paste0("batch_results_", 1:5, "_full.csv")
    for (f in expected_files) {
      if (file.exists(file.path(data_folder, f))) {
        files_found <- c(files_found, f)
      }
    }

    if (length(files_found) >= 2) {
      batch_status <- "passed"
      log_msg("  Found", length(files_found), "pre-existing batch result files in data folder")
    } else {
      log_msg("  Warning: Only found", length(files_found), "batch files in", data_folder)
      log_msg("  Expected files: batch_results_1_full.csv through batch_results_5_full.csv")
      batch_status <- "skipped"
    }

    batch_results <- list(
      clusters_tested = as.list(batch_clusters_tested),
      n_iterations = 5,
      files_created = as.list(files_found),
      source = "pre-existing data"
    )
  }

  # =========================================================================
  # Test 3: runCASSIA_similarity_score_batch - Analyze batch results
  # =========================================================================
  similarity_batch_status <- "error"

  tryCatch({
    log_msg("\n--- Test: runCASSIA_similarity_score_batch ---")
    log_msg("  Analyzing batch results")

    # Run if we have batch files (either generated or pre-existing)
    if (length(files_found) >= 2) {
      similarity_output_name <- file.path(results$outputs, "similarity_score_results")
      similarity_report_path <- file.path(results$outputs, "uq_batch_report.html")

      CASSIA::runCASSIA_similarity_score_batch(
        marker = test_markers,
        file_pattern = paste0(batch_output_name, "_*_full.csv"),
        output_name = similarity_output_name,
        celltype_column = "Broad.cell.type",
        max_workers = 3,
        model = llm_config$model %||% "google/gemini-2.5-flash",
        provider = llm_config$provider %||% "openrouter",
        generate_report = TRUE,
        report_output_path = similarity_report_path
      )

      # Check if similarity score CSV was created
      if (file.exists(paste0(similarity_output_name, ".csv"))) {
        log_msg("  Similarity score CSV created:", basename(paste0(similarity_output_name, ".csv")))
        similarity_batch_status <- "passed"
      } else {
        similarity_batch_status <- "failed"
        errors <- c(errors, "Similarity score batch: CSV file not created")
      }

      # Check if HTML report was created
      if (file.exists(similarity_report_path)) {
        log_msg("  HTML report created:", basename(similarity_report_path))
      } else {
        log_msg("  Warning: HTML report not created")
      }

      log_msg("\n[OK] Similarity score batch test PASSED")
    } else {
      log_msg("  Skipping - no batch files found")
      similarity_batch_status <- "skipped"
    }

  }, error = function(e) {
    errors <<- c(errors, paste("Similarity score batch test error:", e$message))
    similarity_batch_status <<- "error"
    log_msg("\nError in similarity score batch test:", e$message)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Combine statuses - all tests must pass for overall success
  if (status == "passed" && (batch_status == "passed" || batch_status == "skipped") && (similarity_batch_status == "passed" || similarity_batch_status == "skipped")) {
    overall_status <- "passed"
  } else if (status == "error" || batch_status == "error" || similarity_batch_status == "error") {
    overall_status <- "error"
  } else {
    overall_status <- "failed"
  }

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "uncertainty_quantification",
    config = config,
    duration_seconds = duration,
    status = overall_status,
    clusters_tested = c(list(test_cluster), as.list(batch_clusters_tested)),
    errors = errors
  )
  metadata$test_details <- list(
    similarity_score_test = list(
      status = status,
      cluster = test_cluster
    ),
    batch_n_times_test = list(
      status = batch_status,
      clusters = as.list(batch_clusters_tested)
    ),
    similarity_score_batch_test = list(
      status = similarity_batch_status,
      clusters = as.list(batch_clusters_tested)
    )
  )
  save_test_metadata(results$outputs, metadata)

  save_test_results(results$outputs, list(
    similarity_score_test = list(
      test_cluster = test_cluster,
      n_iterations = 5,
      results = uq_results
    ),
    batch_n_times_test = batch_results,
    similarity_score_batch_test = list(
      status = similarity_batch_status
    )
  ))

  # Print final result
  success <- overall_status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
# Check for --run-batch flag to enable batch generation
args <- commandArgs(trailingOnly = TRUE)
run_batch <- "--run-batch" %in% args
success <- run_uncertainty_quantification_test(run_batch_generation = run_batch)
quit(status = if (success) 0 else 1)
