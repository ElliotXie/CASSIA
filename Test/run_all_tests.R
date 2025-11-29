# CASSIA Test Runner - Run All Tests (R)
# =======================================
# Run all CASSIA R tests sequentially and generate a summary report.
#
# Usage:
#     Rscript run_all_tests.R                  # Run with load_all() (fast, default)
#     Rscript run_all_tests.R --install        # Run with install_local() (full package)
#     Rscript run_all_tests.R --skip 03,04
#     Rscript run_all_tests.R --only 01,02

library(jsonlite)

# Get script directory (works with Rscript)
args <- commandArgs(trailingOnly = TRUE)
get_script_dir <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", all_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg), winslash = "/")))
  }
  return(normalizePath(getwd(), winslash = "/"))
}
script_dir <- get_script_dir()

get_test_folders <- function() {
  all_items <- list.dirs(script_dir, recursive = FALSE, full.names = TRUE)
  test_folders <- all_items[grepl("^[0-9]", basename(all_items))]
  sort(test_folders)
}

run_r_test <- function(test_folder, install_mode = FALSE) {
  # Find the R test script based on mode
  if (install_mode) {
    # Look for *_install.R files
    test_scripts <- list.files(test_folder, pattern = "^test_.*_install\\.R$", full.names = TRUE)
    mode_label <- "(install)"
  } else {
    # Look for regular test_*.R files (excluding *_install.R)
    all_scripts <- list.files(test_folder, pattern = "^test_.*\\.R$", full.names = TRUE)
    test_scripts <- all_scripts[!grepl("_install\\.R$", all_scripts)]
    mode_label <- "(load_all)"
  }

  if (length(test_scripts) == 0) {
    return(list(
      name = basename(test_folder),
      status = "skipped",
      reason = paste("No R test script found", mode_label),
      duration = 0,
      mode = if (install_mode) "install" else "load_all"
    ))
  }

  test_script <- test_scripts[1]
  start_time <- Sys.time()

  result <- tryCatch({
    # Run the test in a separate R process
    # Use shQuote to handle paths with spaces
    exit_code <- system2("Rscript", shQuote(test_script), wait = TRUE)
    duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    list(
      name = basename(test_folder),
      status = if (exit_code == 0) "passed" else "failed",
      returncode = exit_code,
      duration = round(duration, 2),
      mode = if (install_mode) "install" else "load_all",
      script = basename(test_script)
    )
  }, error = function(e) {
    list(
      name = basename(test_folder),
      status = "error",
      duration = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      reason = e$message,
      mode = if (install_mode) "install" else "load_all"
    )
  })

  return(result)
}

print_summary <- function(results, total_duration) {
  cat("\n", strrep("=", 60), "\n")
  cat("CASSIA TEST SUITE (R) - SUMMARY\n")
  cat(strrep("=", 60), "\n")

  passed <- sum(sapply(results, function(r) r$status == "passed"))
  failed <- sum(sapply(results, function(r) r$status == "failed"))
  errors <- sum(sapply(results, function(r) r$status %in% c("error", "timeout")))
  skipped <- sum(sapply(results, function(r) r$status == "skipped"))

  cat("\nResults:", passed, "passed,", failed, "failed,", errors, "errors,", skipped, "skipped\n")
  cat("Total Duration:", round(total_duration, 1), "s\n")

  cat("\n", strrep("-", 60), "\n")
  cat(sprintf("%-40s %-10s %-10s\n", "Test", "Status", "Duration"))
  cat(strrep("-", 60), "\n")

  for (result in results) {
    status_symbol <- switch(result$status,
      "passed" = "[OK]",
      "failed" = "[X]",
      "error" = "[ERR]",
      "timeout" = "[TO]",
      "skipped" = "[--]",
      "[?]"
    )
    duration_str <- if (result$duration > 0) paste0(round(result$duration, 1), "s") else "N/A"
    cat(sprintf("%-40s %-10s %-10s\n", result$name, status_symbol, duration_str))
  }

  cat(strrep("-", 60), "\n")

  if (passed == length(results)) {
    cat("\nAll tests PASSED!\n")
  } else if (failed + errors > 0) {
    cat("\n", failed + errors, "test(s) FAILED or had ERRORS\n")
  }

  return(passed == length(results[sapply(results, function(r) r$status != "skipped")]))
}

save_report <- function(results, total_duration) {
  report_path <- file.path(script_dir, "test_report_r.json")

  report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    total_duration = round(total_duration, 2),
    summary = list(
      passed = sum(sapply(results, function(r) r$status == "passed")),
      failed = sum(sapply(results, function(r) r$status == "failed")),
      errors = sum(sapply(results, function(r) r$status %in% c("error", "timeout"))),
      skipped = sum(sapply(results, function(r) r$status == "skipped"))
    ),
    results = results
  )

  jsonlite::write_json(report, report_path, pretty = TRUE, auto_unbox = TRUE)
  cat("\nReport saved to:", report_path, "\n")
}

parse_test_list <- function(test_string) {
  if (is.null(test_string) || test_string == "") {
    return(character(0))
  }
  trimws(strsplit(test_string, ",")[[1]])
}

# Parse arguments
skip_tests <- character(0)
only_tests <- character(0)
no_report <- FALSE
use_install_mode <- FALSE

i <- 1
while (i <= length(args)) {
  if (args[i] == "--skip" && i < length(args)) {
    skip_tests <- parse_test_list(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--only" && i < length(args)) {
    only_tests <- parse_test_list(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--no-report") {
    no_report <- TRUE
    i <- i + 1
  } else if (args[i] == "--install") {
    use_install_mode <- TRUE
    i <- i + 1
  } else {
    i <- i + 1
  }
}

# Main
cat(strrep("=", 60), "\n")
cat("CASSIA TEST SUITE (R)\n")
cat(strrep("=", 60), "\n")
mode_str <- if (use_install_mode) "install_local (full package)" else "load_all (fast dev)"
cat("Mode:", mode_str, "\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Get test folders
test_folders <- get_test_folders()

# Apply filters
if (length(only_tests) > 0) {
  test_folders <- test_folders[sapply(test_folders, function(f) {
    any(sapply(only_tests, function(t) startsWith(basename(f), t)))
  })]
} else if (length(skip_tests) > 0) {
  test_folders <- test_folders[!sapply(test_folders, function(f) {
    any(sapply(skip_tests, function(t) startsWith(basename(f), t)))
  })]
}

cat("\nRunning", length(test_folders), "tests:\n")
for (folder in test_folders) {
  cat("  -", basename(folder), "\n")
}

# Run tests
results <- list()
total_start <- Sys.time()

for (i in seq_along(test_folders)) {
  test_folder <- test_folders[i]

  cat("\n", strrep("=", 60), "\n")
  cat("[", i, "/", length(test_folders), "] Running:", basename(test_folder), "\n")
  cat(strrep("=", 60), "\n")

  result <- run_r_test(test_folder, install_mode = use_install_mode)
  results <- c(results, list(result))

  # Print immediate status
  if (result$status == "passed") {
    cat("\n[OK]", basename(test_folder), "PASSED (", round(result$duration, 1), "s)\n")
  } else {
    cat("\n[X]", basename(test_folder), toupper(result$status), "\n")
    if (!is.null(result$reason)) {
      cat("    Reason:", result$reason, "\n")
    }
  }
}

total_duration <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))

# Print summary
all_passed <- print_summary(results, total_duration)

# Save report
if (!no_report) {
  save_report(results, total_duration)
}

quit(status = if (all_passed) 0 else 1)
