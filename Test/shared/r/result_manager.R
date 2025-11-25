# CASSIA Test Suite - R Results Manager
# ======================================
# Functions for managing test results in R
#
# Note: This file requires test_utils.R to be sourced first (for get_test_root())

library(jsonlite)

#' Create a timestamped results directory
#'
#' @param test_folder Name of the test folder (e.g., '01_single_annotation')
#' @return Path to the created results directory
create_results_dir <- function(test_folder) {
  test_root <- get_test_root()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_dir <- file.path(test_root, test_folder, "results", timestamp)

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  return(results_dir)
}

#' Save test metadata to JSON
#'
#' @param results_dir Path to results directory
#' @param metadata List containing test metadata
save_test_metadata <- function(results_dir, metadata) {
  metadata_path <- file.path(results_dir, "test_metadata.json")
  jsonlite::write_json(metadata, metadata_path, pretty = TRUE, auto_unbox = TRUE)
}

#' Save test results to JSON
#'
#' @param results_dir Path to results directory
#' @param results List containing test results
#' @param filename Name of the output file
save_test_results <- function(results_dir, results, filename = "results.json") {
  results_path <- file.path(results_dir, filename)
  jsonlite::write_json(results, results_path, pretty = TRUE, auto_unbox = TRUE)
}

#' Create test metadata
#'
#' @param test_name Name of the test
#' @param config Test configuration
#' @param duration_seconds Test duration in seconds
#' @param status Test status ('passed', 'failed', 'error')
#' @param clusters_tested List of clusters tested
#' @param errors List of error messages
#' @return Test metadata list
create_test_metadata <- function(test_name, config, duration_seconds, status,
                                  clusters_tested = NULL, errors = NULL) {
  list(
    test_name = test_name,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    language = "R",
    config = list(
      model = config$llm$model,
      provider = config$llm$provider,
      tissue = config$data$tissue,
      species = config$data$species,
      n_genes = config$data$n_genes
    ),
    duration_seconds = round(duration_seconds, 2),
    status = status,
    clusters_tested = clusters_tested %||% list(),
    errors = errors %||% list()
  )
}

#' Cleanup old results directories
#'
#' @param test_folder Name of the test folder
#' @param keep_n Number of recent results to keep
cleanup_old_results <- function(test_folder, keep_n = 10) {
  test_root <- get_test_root()
  results_base <- file.path(test_root, test_folder, "results")

  if (!dir.exists(results_base)) {
    return(invisible(NULL))
  }

  # Get all timestamped directories
  result_dirs <- list.dirs(results_base, recursive = FALSE, full.names = TRUE)
  result_dirs <- result_dirs[order(basename(result_dirs), decreasing = TRUE)]

  # Remove old directories
  if (length(result_dirs) > keep_n) {
    for (old_dir in result_dirs[(keep_n + 1):length(result_dirs)]) {
      unlink(old_dir, recursive = TRUE)
      message("Removed old results: ", basename(old_dir))
    }
  }
}

#' Get the most recent results directory
#'
#' @param test_folder Name of the test folder
#' @return Path to most recent results, or NULL
get_latest_results <- function(test_folder) {
  test_root <- get_test_root()
  results_base <- file.path(test_root, test_folder, "results")

  if (!dir.exists(results_base)) {
    return(NULL)
  }

  result_dirs <- list.dirs(results_base, recursive = FALSE, full.names = TRUE)
  result_dirs <- result_dirs[order(basename(result_dirs), decreasing = TRUE)]

  if (length(result_dirs) > 0) {
    return(result_dirs[1])
  }
  return(NULL)
}

#' List all results directories
#'
#' @param test_folder Name of the test folder
#' @return Vector of results directory paths
list_all_results <- function(test_folder) {
  test_root <- get_test_root()
  results_base <- file.path(test_root, test_folder, "results")

  if (!dir.exists(results_base)) {
    return(character(0))
  }

  result_dirs <- list.dirs(results_base, recursive = FALSE, full.names = TRUE)
  result_dirs[order(basename(result_dirs), decreasing = TRUE)]
}

# Null coalesce operator
`%||%` <- function(x, y) if (is.null(x)) y else x
