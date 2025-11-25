# CASSIA Test 06: Annotation Boost (R)
# =====================================
# Tests the runCASSIA_annotationboost function via R package.
#
# Usage:
#     Rscript test_annotation_boost.R
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

run_annotation_boost_test <- function() {
  print_test_header("06 - Annotation Boost (R)")

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
  results_dir <- create_results_dir("06_annotation_boost")
  cat("Results will be saved to:", results_dir, "\n")

  # Check for existing batch results from Test 02
  batch_results_dir <- get_latest_results("02_batch_annotation")
  batch_results_file <- NULL

  if (!is.null(batch_results_dir)) {
    potential_file <- file.path(batch_results_dir, "batch_results_full.csv")
    if (file.exists(potential_file)) {
      batch_results_file <- potential_file
      cat("\nUsing existing batch results:", batch_results_file, "\n")
    }
  }

  # If no existing results, run batch annotation first
  if (is.null(batch_results_file)) {
    cat("\nNo existing batch results found. Running batch annotation first...\n")
    marker_df <- get_full_marker_dataframe()
    batch_output <- file.path(results_dir, "batch_for_boost")

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

  # Load batch results to get cluster names
  batch_df <- read.csv(batch_results_file)
  cluster_col <- if ("True.Cell.Type" %in% names(batch_df)) "True.Cell.Type" else names(batch_df)[1]
  available_clusters <- batch_df[[cluster_col]]

  # Test cluster for annotation boost
  test_cluster <- "plasma cell"
  if (!(test_cluster %in% available_clusters)) {
    test_cluster <- available_clusters[1]
  }

  cat("\nTesting annotation boost for:", test_cluster, "\n")

  # Get marker data path
  marker_path <- get_marker_file_path()

  # Run annotation boost
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  boost_results <- list()

  output_name <- file.path(results_dir, paste0("boost_", gsub(" ", "_", test_cluster)))

  tryCatch({
    cat("\nRunning annotation boost...\n")
    cat("  Cluster:", test_cluster, "\n")
    cat("  Model:", llm_config$model %||% "google/gemini-2.5-flash", "\n")
    cat("  Provider:", llm_config$provider %||% "openrouter", "\n")
    cat("  Search strategy: breadth\n")
    cat("  Max iterations: 3 (reduced for testing)\n")

    result <- CASSIA::runCASSIA_annotationboost(
      full_result_path = batch_results_file,
      marker = marker_path,
      cluster_name = test_cluster,
      major_cluster_info = paste(data_config$species %||% "human",
                                  data_config$tissue %||% "large intestine"),
      output_name = output_name,
      num_iterations = 3,  # Reduced for faster testing
      model = llm_config$model %||% "google/gemini-2.5-flash",
      provider = llm_config$provider %||% "openrouter",
      temperature = llm_config$temperature %||% 0.3,
      conversation_history_mode = "final",
      search_strategy = "breadth",
      report_style = "per_iteration"
    )

    # Check result
    if (is.list(result)) {
      boost_results <- result

      if (!is.null(result$status) && result$status == "success") {
        cat("\nAnnotation Boost Results:\n")
        cat("  Status:", result$status, "\n")
        cat("  Execution time:", round(result$execution_time %||% 0, 1), "s\n")

        if (!is.null(result$summary_report_path)) {
          cat("  Summary report:", basename(result$summary_report_path), "\n")
        }
        if (!is.null(result$raw_text_path)) {
          cat("  Raw conversation:", basename(result$raw_text_path), "\n")
        }

        # Check if analysis text is valid
        analysis_text <- result$analysis_text %||% ""
        if (nchar(analysis_text) > 100) {
          status <- "passed"
        } else {
          status <- "failed"
          errors <- list("Analysis text is empty or too short")
        }
      } else {
        status <- "failed"
        errors <- list(result$error_message %||% "Unknown error")
      }
    } else {
      status <- "failed"
      errors <- list("Unexpected result format from annotation boost")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    cat("\nError:", e$message, "\n")
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "annotation_boost",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(test_cluster),
    errors = errors
  )
  save_test_metadata(results_dir, metadata)

  save_test_results(results_dir, list(
    cluster = test_cluster,
    batch_results_file = batch_results_file,
    boost_results = list(
      status = boost_results$status,
      execution_time = boost_results$execution_time,
      summary_report = boost_results$summary_report_path,
      raw_text = boost_results$raw_text_path
    )
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  return(success)
}

# Run test
success <- run_annotation_boost_test()
quit(status = if (success) 0 else 1)
