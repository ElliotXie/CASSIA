# CASSIA Test 09: Subclustering (R)
# ==================================
# Tests the subclustering annotation functions via R package.
#
# Usage:
#     Rscript test_subclustering.R
#
# Functions tested:
# - runCASSIA_subclusters(): Single-run subcluster annotation

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

run_subclustering_test <- function() {
  print_test_header("09 - Subclustering (R)")

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
  results_dir <- create_results_dir("09_subclustering")
  cat("Results will be saved to:", results_dir, "\n")

  # Get marker data and create simulated subclusters
  marker_df <- get_full_marker_dataframe()

  # Use first 3 rows as simulated subclusters
  subcluster_df <- marker_df[1:3, ]

  # Define major cluster context
  major_cluster_info <- paste(
    data_config$species %||% "human",
    data_config$tissue %||% "large intestine",
    "immune cells"
  )

  # Run tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  subcluster_results <- list()

  tryCatch({
    cat("\n--- Test: runCASSIA_subclusters ---\n")
    cat("  Major cluster:", major_cluster_info, "\n")
    cat("  Subclusters to annotate:", nrow(subcluster_df), "\n")
    cat("  Model:", llm_config$model %||% "google/gemini-2.5-flash", "\n")
    cat("  Provider:", llm_config$provider %||% "openrouter", "\n")

    output_name <- file.path(results_dir, "subcluster_results")

    CASSIA::runCASSIA_subclusters(
      marker = subcluster_df,
      major_cluster_info = major_cluster_info,
      output_name = output_name,
      model = llm_config$model %||% "google/gemini-2.5-flash",
      temperature = llm_config$temperature %||% 0.3,
      provider = llm_config$provider %||% "openrouter",
      n_genes = data_config$n_genes %||% 30
    )

    # Check output files
    csv_file <- paste0(output_name, ".csv")
    if (file.exists(csv_file)) {
      result_df <- read.csv(csv_file)
      cat("\nSubclustering Results:\n")
      cat("  Output file:", basename(csv_file), "\n")
      cat("  Subclusters annotated:", nrow(result_df), "\n")

      if (nrow(result_df) > 0) {
        cat("\n  Sample annotations:\n")
        for (i in 1:min(3, nrow(result_df))) {
          cat("   ", result_df$Result.ID[i] %||% i, ":",
              result_df$main_cell_type[i] %||% "N/A", "/",
              result_df$sub_cell_type[i] %||% "N/A", "\n")
        }
      }

      subcluster_results <- list(
        output_file = csv_file,
        num_subclusters = nrow(result_df),
        columns = names(result_df)
      )

      # Check for HTML report
      html_file <- paste0(output_name, ".html")
      if (file.exists(html_file)) {
        cat("  HTML report:", basename(html_file), "\n")
        subcluster_results$html_report <- html_file
      }

      status <- "passed"
    } else {
      status <- "failed"
      errors <- list(paste("Output file not created:", csv_file))
    }

  }, error = function(e) {
    errors <<- list(e$message)
    status <<- "error"
    cat("\nError:", e$message, "\n")
  })

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "subclustering",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(major_cluster_info),
    errors = errors
  )
  save_test_metadata(results_dir, metadata)

  save_test_results(results_dir, list(
    major_cluster_info = major_cluster_info,
    num_input_subclusters = nrow(subcluster_df),
    results = subcluster_results
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  return(success)
}

# Run test
success <- run_subclustering_test()
quit(status = if (success) 0 else 1)
