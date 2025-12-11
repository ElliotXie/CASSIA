# CASSIA Test 15: CASSIA Pipeline (R)
# =====================================
# Tests the runCASSIA_pipeline function, the complete end-to-end cell type
# annotation orchestrator via R package.
#
# Usage:
#     Rscript test_cassia_pipeline.R

# Get script directory (works with Rscript and source())
get_script_dir <- function() {
  # Method 1: commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg), winslash = "/")))
  }
  # Method 2: sys.frame (works with source())
  for (i in sys.nframe():1) {
    if (!is.null(sys.frame(i)$ofile)) {
      return(dirname(normalizePath(sys.frame(i)$ofile, winslash = "/")))
    }
  }
  # Fallback: current working directory
  return(normalizePath(getwd(), winslash = "/"))
}
script_dir <- get_script_dir()

# Source shared utilities
source(file.path(script_dir, "..", "shared", "r", "test_utils.R"))
source(file.path(script_dir, "..", "shared", "r", "fixtures.R"))
source(file.path(script_dir, "..", "shared", "r", "result_manager.R"))
source(file.path(script_dir, "..", "shared", "r", "logging_manager.R"))

run_cassia_pipeline_test <- function() {
  print_test_header("15 - CASSIA Pipeline (R)")

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

  # Use only 2 clusters for faster testing
  test_clusters <- c("monocyte", "plasma cell")

  # Load raw FindAllMarkers data from data folder for annotation boost
  data_folder <- file.path(script_dir, "data")
  raw_markers_path <- file.path(data_folder, "findallmarkers_output.csv")
  
  if (!file.exists(raw_markers_path)) {
    stop("Raw marker file not found: ", raw_markers_path)
  }
  
  # Load and filter to only the 2 test clusters
  raw_markers <- read.csv(raw_markers_path)
  marker_df <- raw_markers[raw_markers$cluster %in% test_clusters, ]

  log_msg("\nTesting pipeline for", length(test_clusters), "clusters:")
  for (cluster in test_clusters) {
    log_msg("  -", cluster)
  }
  log_msg("\nLoaded raw marker data:", nrow(marker_df), "rows")
  log_msg("Data source:", raw_markers_path)

  # Create results directory
  results <- create_results_dir("15_cassia_pipeline", get_test_mode())

  start_logging(results$logs)
  output_name <- file.path(results$outputs, "pipeline_test")
  log_msg("Results will be saved to:", results$base)

  # Change to results directory so pipeline output goes there
  original_dir <- getwd()
  setwd(results$outputs)

  # Run the test
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  pipeline_output_dir <- NULL

  tryCatch({
    log_msg("\nRunning runCASSIA_pipeline via R package...")

    # Call CASSIA R function
    CASSIA::runCASSIA_pipeline(
      output_file_name = output_name,
      tissue = data_config$tissue %||% "large intestine",
      species = data_config$species %||% "human",
      marker = marker_df,
      max_workers = llm_config$max_workers %||% 3,
      annotation_model = llm_config$model %||% "google/gemini-2.5-flash",
      annotation_provider = llm_config$provider %||% "openrouter",
      score_model = llm_config$model %||% "google/gemini-2.5-flash",
      score_provider = llm_config$provider %||% "openrouter",
      annotationboost_model = llm_config$model %||% "google/gemini-2.5-flash",
      annotationboost_provider = llm_config$provider %||% "openrouter",
      score_threshold = 99,
      do_merge_annotations = TRUE,
      merge_model = llm_config$model %||% "google/gemini-2.5-flash",
      merge_provider = llm_config$provider %||% "openrouter",
      validator_involvement = config$validator$default %||% "v1"
    )

    # Find the pipeline output directory (starts with CASSIA_)
    items <- list.dirs(results$outputs, recursive = FALSE, full.names = FALSE)
    cassia_dirs <- items[grepl("^CASSIA_", items)]
    if (length(cassia_dirs) > 0) {
      pipeline_output_dir <- file.path(results$outputs, cassia_dirs[1])
    }

    if (!is.null(pipeline_output_dir) && dir.exists(pipeline_output_dir)) {
      log_msg("\nPipeline output directory:", basename(pipeline_output_dir))

      # Check expected subdirectories
      annotation_dir <- file.path(pipeline_output_dir, "01_annotation_results")
      reports_dir_path <- file.path(pipeline_output_dir, "02_reports")
      boost_dir <- file.path(pipeline_output_dir, "03_boost_analysis")

      checks_passed <- TRUE
      log_msg("\nValidating output structure:")

      # Check 01_annotation_results
      if (dir.exists(annotation_dir)) {
        log_msg("  [OK] 01_annotation_results exists")
        # Check for FINAL_RESULTS.csv
        final_results <- file.path(annotation_dir, "FINAL_RESULTS.csv")
        if (file.exists(final_results)) {
          log_msg("  [OK] FINAL_RESULTS.csv exists")
          results_df <- read.csv(final_results)
          log_msg("       - Contains", nrow(results_df), "rows")
        } else {
          log_msg("  [WARN] FINAL_RESULTS.csv not found")
        }
      } else {
        log_msg("  [FAIL] 01_annotation_results missing")
        checks_passed <- FALSE
        errors <- list("01_annotation_results directory missing")
      }

      # Check 02_reports
      if (dir.exists(reports_dir_path)) {
        log_msg("  [OK] 02_reports exists")
        html_files <- list.files(reports_dir_path, pattern = "\\.html$")
        log_msg("       - Contains", length(html_files), "HTML report(s)")
      } else {
        log_msg("  [WARN] 02_reports missing (may be expected if no reports generated)")
      }

      # Check 03_boost_analysis
      if (dir.exists(boost_dir)) {
        log_msg("  [OK] 03_boost_analysis exists")
      } else {
        log_msg("  [INFO] 03_boost_analysis missing (expected if all scores above threshold)")
      }

      if (checks_passed) {
        status <- "passed"
      } else {
        status <- "failed"
      }
    } else {
      status <- "failed"
      errors <- list("Pipeline output directory not created")
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    log_msg("\nError:", e$message)
  })

  # Change back to original directory
  setwd(original_dir)

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata
  metadata <- create_test_metadata(
    test_name = "cassia_pipeline",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(test_clusters),
    errors = errors
  )
  save_test_metadata(results$outputs, metadata)

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_cassia_pipeline_test()
quit(status = if (success) 0 else 1)
