# CASSIA Test 06: Annotation Boost (R) - INSTALL MODE
# ====================================================
# Tests the runCASSIA_annotationboost function via R package.
# Uses devtools::install_local() for full package installation testing.
#
# Usage:
#     Rscript test_annotation_boost_install.R
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

run_annotation_boost_test <- function() {
  print_test_header("06 - Annotation Boost (R) [INSTALL MODE]")

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
  results_dirs <- create_results_dir("06_annotation_boost", get_test_mode())

  start_logging(results_dirs$logs)
  log_msg("Results will be saved to:", results_dirs$base)

  # Check for existing batch results from Test 02
  batch_results_dir <- get_latest_results("02_batch_annotation")
  batch_results_file <- NULL

  if (!is.null(batch_results_dir)) {
    potential_file <- file.path(batch_results_dir, "batch_results_summary.csv")
    if (file.exists(potential_file)) {
      batch_results_file <- potential_file
      log_msg("\nUsing existing batch results:", batch_results_file)
    }
  }

  # If no existing results, run batch annotation first
  if (is.null(batch_results_file)) {
    log_msg("\nNo existing batch results found. Running batch annotation first...")
    marker_df <- get_full_marker_dataframe()
    batch_output <- file.path(results_dirs$outputs, "batch_for_boost")

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
    batch_results_file <- paste0(batch_output, "_summary.csv")
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

  log_msg("\nTesting annotation boost for:", test_cluster)

  # Get marker data path - annotation_boost requires RAW FindAllMarkers format
  # Use local findallmarkers_output.csv instead of processed.csv
  marker_path <- file.path(script_dir, "data", "findallmarkers_output.csv")

  # Run annotation boost
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  boost_results <- list()

  output_name <- file.path(results_dirs$outputs, paste0("boost_", gsub(" ", "_", test_cluster)))

  tryCatch({
    log_msg("\nRunning annotation boost (install mode)...")
    log_msg("  Cluster:", test_cluster)
    log_msg("  Model:", llm_config$model %||% "google/gemini-2.5-flash")
    log_msg("  Provider:", llm_config$provider %||% "openrouter")
    log_msg("  Search strategy: breadth")
    log_msg("  Max iterations: 3 (reduced for testing)")

    # Run annotation boost (returns NULL, outputs are files)
    CASSIA::runCASSIA_annotationboost(
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

    # Check output files exist (function outputs files, not return values)
    summary_report_path <- paste0(output_name, "_summary.html")
    raw_text_path <- paste0(output_name, "_raw_conversation.txt")

    log_msg("\nAnnotation Boost Results:")

    if (file.exists(summary_report_path)) {
      file_size <- file.info(summary_report_path)$size
      log_msg("  Summary report:", basename(summary_report_path), "(", round(file_size/1024, 1), "KB )")
      boost_results$summary_report_path <- summary_report_path
    }

    if (file.exists(raw_text_path)) {
      # Read raw conversation to check content
      raw_content <- readLines(raw_text_path, warn = FALSE)
      raw_text <- paste(raw_content, collapse = "\n")
      file_size <- file.info(raw_text_path)$size
      log_msg("  Raw conversation:", basename(raw_text_path), "(", round(file_size/1024, 1), "KB )")
      boost_results$raw_text_path <- raw_text_path

      # Check if raw conversation has meaningful content (>100 chars)
      if (nchar(raw_text) > 100) {
        status <- "passed"
        log_msg("  Content length:", nchar(raw_text), "characters")
      } else {
        status <- "failed"
        errors <- list("Raw conversation text is empty or too short")
      }
    } else {
      status <- "failed"
      errors <- list("Output files not generated")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    log_msg("\nError:", e$message)
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "annotation_boost_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(test_cluster),
    errors = errors
  )
  save_test_metadata(results_dirs$outputs, metadata)

  save_test_results(results_dirs$outputs, list(
    cluster = test_cluster,
    batch_results_file = batch_results_file,
    boost_results = list(
      status = boost_results$status,
      execution_time = boost_results$execution_time,
      summary_report = boost_results$summary_report_path,
      raw_text = boost_results$raw_text_path
    ),
    mode = "install"
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_annotation_boost_test()
quit(status = if (success) 0 else 1)
