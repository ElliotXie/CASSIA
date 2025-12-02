# CASSIA Test 11: Report Generation (R) - INSTALL MODE
# =====================================================
# Tests the generate_batch_html_report function for generating interactive HTML
# reports from batch analysis results WITHOUT re-running the analysis.
# Uses devtools::install_local() for full package installation testing.
#
# Usage:
#     Rscript test_report_generation_install.R
#
# Functions tested:
# - generate_batch_html_report(): Generate HTML report from CSV file (via Python)

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

find_existing_batch_results <- function() {
  # Find existing batch results from Test 02
  test_02_dir <- file.path(script_dir, "..", "02_batch_annotation", "results")

  if (!dir.exists(test_02_dir)) {
    return(NULL)
  }

  # Find most recent results folder
  result_folders <- list.dirs(test_02_dir, recursive = FALSE)
  result_folders <- sort(result_folders, decreasing = TRUE)

  for (folder in result_folders) {
    full_csv <- file.path(folder, "batch_results_full.csv")
    if (file.exists(full_csv)) {
      return(full_csv)
    }
  }

  return(NULL)
}

create_sample_batch_data <- function() {
  # Create sample batch data for testing if no existing results
  cluster_names <- get_all_clusters()
  sample_data <- list()

  for (i in 1:min(3, length(cluster_names))) {
    cluster_name <- cluster_names[i]
    markers <- get_cluster_markers(cluster_name)
    marker_list <- paste(markers[1:min(15, length(markers))], collapse = ", ")

    sample_data[[i]] <- list(
      `Cluster ID` = cluster_name,
      `Predicted General Cell Type` = tools::toTitleCase(cluster_name),
      `Predicted Detailed Cell Type` = paste0(cluster_name, " subtype 1, ", cluster_name, " subtype 2"),
      `Possible Mixed Cell Types` = "",
      `Marker Number` = as.character(min(15, length(markers))),
      `Marker List` = marker_list,
      `Iterations` = "1",
      `Model` = "google/gemini-2.5-flash",
      `Provider` = "openrouter",
      `Tissue` = "large intestine",
      `Species` = "human",
      `Additional Info` = "",
      `Conversation History` = paste0(
        "Final Annotation Agent: Analysis of ", cluster_name, " markers shows characteristic expression pattern. | ",
        "Coupling Validator: VALIDATION PASSED - Cell type identification is consistent. | ",
        'Formatting Agent: {"main_cell_type": "', tools::toTitleCase(cluster_name), '", "sub_cell_types": ["', cluster_name, ' subtype 1", "', cluster_name, ' subtype 2"], "possible_mixed_cell_types": []}'
      )
    )
  }

  # Convert to data frame
  df <- do.call(rbind, lapply(sample_data, as.data.frame, stringsAsFactors = FALSE))
  return(df)
}

run_report_generation_install_test <- function() {
  print_test_header("11 - Report Generation (R) [INSTALL MODE]")

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
  results_dirs <- create_results_dir("11_report_generation", get_test_mode())

  start_logging(results_dirs$logs)
  log_msg("Results will be saved to:", results_dirs$base)

  # Run tests
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  report_results <- list()

  tryCatch({
    # Import the reports module directly from CASSIA Python package
    # Note: py_cassia is internal to the CASSIA package, not in global env
    reports <- reticulate::import("CASSIA.reports.generate_batch_report")

    log_msg("\n--- Test 1: Find existing batch results ---")
    existing_csv <- find_existing_batch_results()

    if (!is.null(existing_csv)) {
      log_msg("  Found existing results:", basename(existing_csv))

      # Read CSV to count rows
      csv_data <- read.csv(existing_csv, stringsAsFactors = FALSE)
      log_msg("  CSV contains:", nrow(csv_data), "clusters")

      # Generate report from existing CSV
      log_msg("\n--- Test 2: Generate HTML report from CSV ---")
      output_html <- file.path(results_dirs$outputs, "test_report_from_csv.html")

      result_path <- reports$generate_batch_html_report(
        full_csv_path = existing_csv,
        output_path = output_html,
        report_title = "CASSIA Test Report (from CSV via R) [INSTALL MODE]"
      )

      if (file.exists(result_path)) {
        file_size <- file.info(result_path)$size
        log_msg("  Generated:", basename(result_path))
        log_msg("  File size:", round(file_size / 1024, 1), "KB")

        report_results$csv_report <- list(
          source = "existing_csv",
          input_file = existing_csv,
          output_file = result_path,
          file_size_bytes = file_size,
          clusters = nrow(csv_data)
        )
      } else {
        errors <- c(errors, "CSV report not created")
      }
    } else {
      log_msg("  No existing batch results found")
      log_msg("  Will create sample data for testing")
    }

    # Test 3: Generate report from sample data
    log_msg("\n--- Test 3: Generate HTML report from data ---")
    sample_data <- create_sample_batch_data()
    log_msg("  Created sample data:", nrow(sample_data), "clusters")

    # Save sample data to CSV first (required by Python function)
    sample_csv_path <- file.path(results_dirs$outputs, "sample_batch_results.csv")
    write.csv(sample_data, sample_csv_path, row.names = FALSE)

    output_html_data <- file.path(results_dirs$outputs, "test_report_from_data.html")

    result_path_data <- reports$generate_batch_html_report(
      full_csv_path = sample_csv_path,
      output_path = output_html_data,
      report_title = "CASSIA Test Report (from R Data) [INSTALL MODE]"
    )

    if (file.exists(result_path_data)) {
      file_size_data <- file.info(result_path_data)$size
      log_msg("  Generated:", basename(result_path_data))
      log_msg("  File size:", round(file_size_data / 1024, 1), "KB")

      report_results$data_report <- list(
        source = "sample_data",
        output_file = result_path_data,
        file_size_bytes = file_size_data,
        clusters = nrow(sample_data)
      )
    } else {
      errors <- c(errors, "Data report not created")
    }

    # Test 4: Validate HTML content
    log_msg("\n--- Test 4: Validate HTML content ---")
    report_to_check <- if (file.exists(result_path_data)) {
      result_path_data
    } else if (!is.null(report_results$csv_report)) {
      report_results$csv_report$output_file
    } else {
      NULL
    }

    if (!is.null(report_to_check) && file.exists(report_to_check)) {
      html_content <- readLines(report_to_check, warn = FALSE)
      html_content <- paste(html_content, collapse = "\n")

      # Check for key HTML elements
      checks <- list(
        DOCTYPE = grepl("<!DOCTYPE html>", html_content),
        report_header = grepl("report-header", html_content),
        cluster_cards = grepl("cluster-card", html_content),
        search_functionality = grepl("search-input", html_content),
        filter_dropdowns = grepl("filter-select", html_content),
        modal_popups = grepl("modal-overlay", html_content),
        javascript = grepl("<script>", html_content),
        css_styles = grepl("<style>", html_content)
      )

      log_msg("  Validation results:")
      all_passed <- TRUE
      for (check_name in names(checks)) {
        passed <- checks[[check_name]]
        status_str <- if (passed) "PASS" else "FAIL"
        log_msg("   ", check_name, ":", status_str)
        if (!passed) all_passed <- FALSE
      }

      report_results$validation <- list(
        file_checked = report_to_check,
        checks = checks,
        all_passed = all_passed
      )

      if (all_passed) {
        status <<- "passed"
      } else {
        status <<- "failed"
        errors <<- c(errors, "Some HTML validation checks failed")
      }
    } else {
      status <- "failed"
      errors <- c(errors, "No HTML file to validate")
    }

    # Summary
    log_msg("\n--- Report Generation Summary ---")
    if (!is.null(report_results$csv_report)) {
      log_msg("  CSV Report:", basename(report_results$csv_report$output_file))
    }
    if (!is.null(report_results$data_report)) {
      log_msg("  Data Report:", basename(report_results$data_report$output_file))
    }
    if (!is.null(report_results$validation)) {
      log_msg("  Validation:", if (report_results$validation$all_passed) "PASSED" else "FAILED")
    }

  }, error = function(e) {
    errors <<- c(errors, e$message)
    status <<- "error"
    log_msg("\nError:", e$message)
  })

  # Set final status based on validation results (workaround for tryCatch scoping)
  if (!is.null(report_results$validation) && report_results$validation$all_passed) {
    status <- "passed"
  } else if (status != "error") {
    status <- "failed"
  }

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Save metadata and results
  metadata <- create_test_metadata(
    test_name = "report_generation_install",
    config = config,
    duration_seconds = duration,
    status = status,
    clusters_tested = list(),
    errors = errors
  )
  save_test_metadata(results_dirs$outputs, metadata)

  save_test_results(results_dirs$outputs, list(
    report_results = report_results,
    mode = "install"
  ))

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))


  stop_logging()
  return(success)
}

# Run test
success <- run_report_generation_install_test()
quit(status = if (success) 0 else 1)
