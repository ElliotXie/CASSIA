# CASSIA Test Suite - R Test Utilities
# =====================================
# Common utilities for running R tests

library(yaml)

#' Get the test suite root directory
#'
#' @return Character path to test root
get_test_root <- function() {
  # Try multiple strategies to find the Test root directory

  # Strategy 1: Check if TEST_ROOT environment variable is set
  env_root <- Sys.getenv("CASSIA_TEST_ROOT", unset = "")
  if (env_root != "" && dir.exists(env_root)) {
    return(normalizePath(env_root, winslash = "/"))
  }

  # Strategy 2: Try sys.frame approach (works with Rscript)
  tryCatch({
    for (i in 1:sys.nframe()) {
      frame_file <- sys.frame(i)$ofile
      if (!is.null(frame_file) && file.exists(frame_file)) {
        # Check if this looks like a test file or test_utils.R
        if (grepl("test_utils\\.R$", frame_file)) {
          # This file is in Test/shared/r/
          candidate <- normalizePath(file.path(dirname(frame_file), "..", ".."), winslash = "/")
          if (file.exists(file.path(candidate, "config", "test_config.yaml"))) {
            return(candidate)
          }
        } else if (grepl("test_.*\\.R$", frame_file)) {
          # Test file is in Test/XX_testname/
          candidate <- normalizePath(file.path(dirname(frame_file), ".."), winslash = "/")
          if (file.exists(file.path(candidate, "config", "test_config.yaml"))) {
            return(candidate)
          }
        }
      }
    }
  }, error = function(e) NULL)

  # Strategy 3: Search upward from current working directory
  current <- getwd()
  for (i in 1:10) {
    # Check if current directory is the Test root
    if (file.exists(file.path(current, "config", "test_config.yaml"))) {
      return(normalizePath(current, winslash = "/"))
    }
    # Check if Test subdirectory exists
    test_dir <- file.path(current, "Test")
    if (file.exists(file.path(test_dir, "config", "test_config.yaml"))) {
      return(normalizePath(test_dir, winslash = "/"))
    }
    # Go up one level
    parent <- dirname(current)
    if (parent == current) break  # Reached root
    current <- parent
  }

  # Strategy 4: Check common locations relative to CASSIA package
  common_paths <- c(
    file.path(getwd(), "Test"),
    file.path(getwd(), "..", "Test"),
    file.path(getwd(), "..", "..", "Test"),
    "C:/Users/ellio/OneDrive - UW-Madison/CASSIA_enjoy/CASSIA/Test"
  )

  for (path in common_paths) {
    if (file.exists(file.path(path, "config", "test_config.yaml"))) {
      return(normalizePath(path, winslash = "/"))
    }
  }

  stop("Could not find Test root directory. Set CASSIA_TEST_ROOT environment variable or run from the Test directory.")
}

#' Get path to CASSIA R package
#'
#' @return Character path to CASSIA_R directory
get_cassia_r_path <- function() {
  # Test root is in CASSIA/Test, so CASSIA_R is at CASSIA/CASSIA_R
  test_root <- get_test_root()
  cassia_r <- file.path(dirname(test_root), "CASSIA_R")

  if (!dir.exists(cassia_r)) {
    # Try alternative: maybe we're already in CASSIA directory
    alt_path <- file.path(test_root, "..", "CASSIA_R")
    if (dir.exists(alt_path)) {
      cassia_r <- normalizePath(alt_path, winslash = "/")
    } else {
      stop(paste("CASSIA_R directory not found at:", cassia_r))
    }
  }

  normalizePath(cassia_r, winslash = "/")
}

#' Setup CASSIA R package using devtools::load_all() (Development Mode)
#'
#' Uses load_all() for fast development testing without package installation.
#' This is the DEFAULT and RECOMMENDED method for development.
#'
#' @return TRUE if successful
setup_cassia_dev <- function() {
  cassia_path <- get_cassia_r_path()

  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools package is required for development mode testing")
  }

  message("Loading CASSIA R package using devtools::load_all()...")

  # Unload CASSIA if previously loaded to avoid conflicts
  if ("CASSIA" %in% loadedNamespaces()) {
    try(unloadNamespace("CASSIA"), silent = TRUE)
  }

  devtools::load_all(cassia_path, export_all = FALSE)
  return(TRUE)
}

#' Setup CASSIA R package using devtools::install_local() (Install Mode)
#'
#' Installs CASSIA from local source. Use this for full package installation testing.
#'
#' @param force Force reinstall even if already installed (default: FALSE)
#' @return TRUE if successful
setup_cassia_install <- function(force = FALSE) {
  cassia_path <- get_cassia_r_path()

  # Check if CASSIA is already installed
  if (force || !requireNamespace("CASSIA", quietly = TRUE)) {
    message("Installing CASSIA R package from source...")
    devtools::install_local(cassia_path, force = TRUE, quiet = TRUE)
  }

  library(CASSIA)
  return(TRUE)
}

#' Setup CASSIA R package (default: development mode)
#'
#' Wrapper that defaults to setup_cassia_dev() for fast iteration.
#' Use setup_cassia_install() for full package installation testing.
#'
#' @return TRUE if successful
setup_cassia <- function() {
  setup_cassia_dev()
}

#' Load test configuration from YAML
#'
#' @return List containing configuration
load_config <- function() {
  test_root <- get_test_root()
  config_path <- file.path(test_root, "config", "test_config.yaml")

  if (!file.exists(config_path)) {
    stop(paste("Config file not found:", config_path))
  }

  yaml::read_yaml(config_path)
}

#' Setup API keys from environment file
#'
#' @return TRUE if file exists, FALSE otherwise
setup_api_keys <- function() {
  test_root <- get_test_root()
  env_path <- file.path(test_root, "config", "api_keys.env")

  if (file.exists(env_path)) {
    # Read and set environment variables
    lines <- readLines(env_path)
    for (line in lines) {
      if (grepl("=", line) && !grepl("^#", line)) {
        parts <- strsplit(line, "=")[[1]]
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = "="))
        Sys.setenv(key = value)
        do.call(Sys.setenv, setNames(list(value), key))
      }
    }
    return(TRUE)
  }
  return(FALSE)
}

#' Get LLM configuration
#'
#' @return List with provider, model, temperature, etc.
get_llm_config <- function() {
  config <- load_config()
  config$llm
}

#' Validate annotation result structure
#'
#' @param result The annotation result
#' @return List with 'valid' and 'errors'
validate_result <- function(result) {
  errors <- c()

  if (is.null(result$main_cell_type) || result$main_cell_type == "") {
    errors <- c(errors, "Missing or empty main_cell_type")
  }

  list(
    valid = length(errors) == 0,
    errors = errors,
    has_main_cell_type = !is.null(result$main_cell_type) && result$main_cell_type != "",
    has_sub_cell_types = !is.null(result$sub_cell_types) && length(result$sub_cell_types) > 0
  )
}

#' Print test header
#'
#' @param test_name Name of the test
print_test_header <- function(test_name) {
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  CASSIA TEST:", test_name, "\n")
  cat(strrep("=", 60), "\n")
}

#' Print test result
#'
#' @param success Boolean indicating test success
#' @param message Optional message to display
print_test_result <- function(success, message = "") {
  status <- if (success) "PASSED" else "FAILED"
  symbol <- if (success) "[OK]" else "[X]"
  cat("\n", symbol, "Test", status, "\n")
  if (message != "") {
    cat("    ", message, "\n")
  }
}

#' Print configuration summary
#'
#' @param config Configuration list
print_config_summary <- function(config) {
  llm <- config$llm
  data <- config$data
  cat("\nConfiguration:\n")
  cat("  Provider:", llm$provider %||% "N/A", "\n")
  cat("  Model:", llm$model %||% "N/A", "\n")
  cat("  Tissue:", data$tissue %||% "N/A", "\n")
  cat("  Species:", data$species %||% "N/A", "\n")
}

# Null coalesce operator
`%||%` <- function(x, y) if (is.null(x)) y else x
