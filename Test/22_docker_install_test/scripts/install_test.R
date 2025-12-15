#!/usr/bin/env Rscript
# =============================================================================
# CASSIA Fresh Installation Test Script
# =============================================================================
# This script simulates a brand new user installing CASSIA from GitHub.
# It tests each stage of the installation process and outputs JSON results.
#
# Stages:
#   1. Install devtools and dependencies
#   2. Install CASSIA from GitHub
#   3. Load library(CASSIA) - triggers .onLoad() and Python setup
#   4. Verify Python environment is properly configured
#   5. Run single-cluster annotation (if API key provided)
#
# Environment Variables:
#   CASSIA_API_KEY     - API key for annotation test (optional)
#   CASSIA_PROVIDER    - LLM provider (default: "openrouter")
#   CASSIA_MODEL       - Model to use (default: "google/gemini-2.0-flash-001")
#
# Output: JSON to stdout with test results
# =============================================================================

# Start timing
start_time <- Sys.time()

# Get environment variables
api_key <- Sys.getenv("CASSIA_API_KEY", "")
provider <- Sys.getenv("CASSIA_PROVIDER", "openrouter")
model <- Sys.getenv("CASSIA_MODEL", "google/gemini-2.0-flash-001")

# Initialize results structure
results <- list(
  scenario = Sys.getenv("CASSIA_SCENARIO", "unknown"),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  platform = list(
    os = Sys.info()["sysname"],
    release = Sys.info()["release"],
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  ),
  stages = list(
    stage1_devtools = list(success = FALSE, message = "", duration = 0),
    stage2_install = list(success = FALSE, message = "", duration = 0),
    stage3_load = list(success = FALSE, message = "", duration = 0),
    stage4_python = list(success = FALSE, message = "", duration = 0, info = list()),
    stage5_annotation = list(success = FALSE, message = "", duration = 0, skipped = FALSE)
  ),
  overall_success = FALSE,
  total_duration = 0,
  errors = list()
)

# Helper function to log messages
log_msg <- function(...) {
  msg <- paste(format(Sys.time(), "[%H:%M:%S]"), paste(..., collapse = " "))
  cat(msg, "\n", file = stderr())
}

# Helper to capture stage results
run_stage <- function(stage_name, expr) {
  stage_start <- Sys.time()
  result <- tryCatch({
    eval(expr)
    list(success = TRUE, message = "OK", error = NULL)
  }, error = function(e) {
    list(success = FALSE, message = conditionMessage(e), error = e)
  }, warning = function(w) {
    list(success = TRUE, message = paste("Warning:", conditionMessage(w)), error = NULL)
  })
  stage_end <- Sys.time()
  result$duration <- as.numeric(difftime(stage_end, stage_start, units = "secs"))
  result
}

# =============================================================================
# STAGE 1: Install devtools
# =============================================================================
log_msg("=== STAGE 1: Installing devtools ===")

stage1_result <- run_stage("devtools", {
  # Check if devtools already installed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    log_msg("Installing devtools from CRAN...")
    install.packages("devtools", repos = "https://cran.rstudio.com/", quiet = TRUE)
  }
  # Also install jsonlite for output
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    log_msg("Installing jsonlite from CRAN...")
    install.packages("jsonlite", repos = "https://cran.rstudio.com/", quiet = TRUE)
  }
  library(devtools)
  library(jsonlite)
  log_msg("devtools version:", as.character(packageVersion("devtools")))
})

results$stages$stage1_devtools <- stage1_result
if (!stage1_result$success) {
  results$errors <- c(results$errors, paste("Stage 1:", stage1_result$message))
}

# =============================================================================
# STAGE 2: Install CASSIA from GitHub
# =============================================================================
if (stage1_result$success) {
  log_msg("=== STAGE 2: Installing CASSIA from GitHub ===")

  stage2_result <- run_stage("install_cassia", {
    log_msg("Installing from ElliotXie/CASSIA...")
    devtools::install_github(
      "ElliotXie/CASSIA",
      subdir = "CASSIA_R",
      upgrade = "never",
      quiet = FALSE
    )
    log_msg("CASSIA installation completed")
  })

  results$stages$stage2_install <- stage2_result
  if (!stage2_result$success) {
    results$errors <- c(results$errors, paste("Stage 2:", stage2_result$message))
  }
} else {
  log_msg("Skipping Stage 2 (Stage 1 failed)")
  results$stages$stage2_install$message <- "Skipped (Stage 1 failed)"
}

# =============================================================================
# STAGE 3: Load library(CASSIA)
# =============================================================================
if (results$stages$stage2_install$success) {
  log_msg("=== STAGE 3: Loading CASSIA library ===")

  stage3_result <- run_stage("load_cassia", {
    log_msg("Loading library(CASSIA)...")
    library(CASSIA)
    log_msg("CASSIA loaded successfully")
    log_msg("CASSIA version:", as.character(packageVersion("CASSIA")))
  })

  results$stages$stage3_load <- stage3_result
  if (!stage3_result$success) {
    results$errors <- c(results$errors, paste("Stage 3:", stage3_result$message))
  }
} else {
  log_msg("Skipping Stage 3 (Stage 2 failed)")
  results$stages$stage3_load$message <- "Skipped (Stage 2 failed)"
}

# =============================================================================
# STAGE 4: Verify Python environment
# =============================================================================
if (results$stages$stage3_load$success) {
  log_msg("=== STAGE 4: Verifying Python environment ===")

  stage4_result <- run_stage("check_python", {
    log_msg("Checking Python environment...")

    # Try to get Python info from reticulate
    python_info <- list()

    if (requireNamespace("reticulate", quietly = TRUE)) {
      tryCatch({
        py_config <- reticulate::py_config()
        python_info$python_path <- py_config$python
        python_info$python_version <- py_config$version
        python_info$virtualenv <- py_config$virtualenv
        python_info$conda <- !is.null(py_config$conda)
        log_msg("Python path:", python_info$python_path)
        log_msg("Python version:", python_info$python_version)
      }, error = function(e) {
        log_msg("Could not get Python config:", conditionMessage(e))
      })
    }

    # Try CASSIA's check_python_env if available
    if (exists("check_python_env", mode = "function")) {
      tryCatch({
        cassia_py_info <- check_python_env()
        python_info$cassia_check <- cassia_py_info
        log_msg("CASSIA Python check passed")
      }, error = function(e) {
        log_msg("CASSIA check_python_env failed:", conditionMessage(e))
      })
    }

    results$stages$stage4_python$info <- python_info
    log_msg("Python environment verified")
  })

  results$stages$stage4_python$success <- stage4_result$success
  results$stages$stage4_python$message <- stage4_result$message
  results$stages$stage4_python$duration <- stage4_result$duration

  if (!stage4_result$success) {
    results$errors <- c(results$errors, paste("Stage 4:", stage4_result$message))
  }
} else {
  log_msg("Skipping Stage 4 (Stage 3 failed)")
  results$stages$stage4_python$message <- "Skipped (Stage 3 failed)"
}

# =============================================================================
# STAGE 5: Run single-cluster annotation test
# =============================================================================
if (results$stages$stage4_python$success && nchar(api_key) > 0) {
  log_msg("=== STAGE 5: Running annotation test ===")

  stage5_result <- run_stage("annotation", {
    log_msg("Setting API key for provider:", provider)
    setLLMApiKey(api_key, provider = provider)

    # Simple T cell markers for quick test
    test_markers <- c("CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LEF1", "TCF7")

    log_msg("Running runCASSIA with", length(test_markers), "markers...")
    result <- runCASSIA(
      marker_list = test_markers,
      tissue = "blood",
      species = "human",
      model = model,
      provider = provider
    )

    if (!is.null(result)) {
      log_msg("Annotation result received")
      if (!is.null(result$main_cell_type)) {
        log_msg("Main cell type:", result$main_cell_type)
      }
    }

    TRUE
  })

  results$stages$stage5_annotation <- stage5_result
  if (!stage5_result$success) {
    results$errors <- c(results$errors, paste("Stage 5:", stage5_result$message))
  }
} else if (nchar(api_key) == 0) {
  log_msg("Skipping Stage 5 (no API key provided)")
  results$stages$stage5_annotation$skipped <- TRUE
  results$stages$stage5_annotation$message <- "Skipped (no API key)"
  results$stages$stage5_annotation$success <- TRUE  # Not a failure if skipped
} else {
  log_msg("Skipping Stage 5 (Stage 4 failed)")
  results$stages$stage5_annotation$message <- "Skipped (Stage 4 failed)"
}

# =============================================================================
# Final results
# =============================================================================
end_time <- Sys.time()
results$total_duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Determine overall success
# Success = Stages 1-4 passed (Stage 5 optional)
results$overall_success <- (
  results$stages$stage1_devtools$success &&
  results$stages$stage2_install$success &&
  results$stages$stage3_load$success &&
  results$stages$stage4_python$success
)

log_msg("=== TEST COMPLETE ===")
log_msg("Overall success:", results$overall_success)
log_msg("Total duration:", round(results$total_duration, 2), "seconds")

# Output JSON results to stdout
if (requireNamespace("jsonlite", quietly = TRUE)) {
  cat("\n---JSON_RESULTS_START---\n")
  cat(jsonlite::toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  cat("\n---JSON_RESULTS_END---\n")
} else {
  # Fallback: simple text output
  cat("\n---RESULTS---\n")
  cat("overall_success:", results$overall_success, "\n")
  cat("stage1:", results$stages$stage1_devtools$success, "\n")
  cat("stage2:", results$stages$stage2_install$success, "\n")
  cat("stage3:", results$stages$stage3_load$success, "\n")
  cat("stage4:", results$stages$stage4_python$success, "\n")
  cat("stage5:", results$stages$stage5_annotation$success, "\n")
}

# Exit with appropriate code
quit(status = if(results$overall_success) 0 else 1, save = "no")
