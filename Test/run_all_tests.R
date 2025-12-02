# CASSIA Test Runner - Run All Tests (R)
# =======================================
# Run all CASSIA R tests sequentially and generate a summary report.
#
# Usage:
#     Rscript run_all_tests.R                  # Run with load_all() (fast, default)
#     Rscript run_all_tests.R --install        # Run with install_local() (full package)
#     Rscript run_all_tests.R --skip 03,04
#     Rscript run_all_tests.R --only 01,02
#     & "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" "c:\Users\ellio\OneDrive - UW-Madison\CASSIA_enjoy\CASSIA\Test\run_all_tests.R" --install

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

  # Determine test mode from install_mode parameter
  current_test_mode <- if (install_mode) "installed" else "development"

  # Check if CASSIA was pre-installed by the test runner
  preinstalled <- Sys.getenv("CASSIA_PREINSTALLED", unset = "")

  result <- tryCatch({
    # Run the test in a separate R process and capture output
    # Use shQuote to handle paths with spaces
    # Set CASSIA_TEST_MODE and CASSIA_PREINSTALLED env vars so child process inherits them

    # Create temp files to capture stdout/stderr
    stdout_file <- tempfile("stdout_", fileext = ".txt")
    stderr_file <- tempfile("stderr_", fileext = ".txt")

    if (.Platform$OS.type == "windows") {
      # On Windows, use cmd /c to run the command with the env vars set
      # Note: No space before && to avoid trailing space in env var value
      env_vars <- paste0('set CASSIA_TEST_MODE=', current_test_mode)
      if (preinstalled == "TRUE") {
        env_vars <- paste0(env_vars, '&& set CASSIA_PREINSTALLED=TRUE')
      }
      cmd <- paste0('cmd /c "', env_vars, '&& Rscript ', shQuote(test_script), '"')
      exit_code <- system2("cmd", args = c("/c", paste0('"', env_vars, '&& Rscript ', shQuote(test_script), '"')),
                           stdout = stdout_file, stderr = stderr_file, wait = TRUE)
    } else {
      # On Unix, set env vars inline
      Sys.setenv(CASSIA_TEST_MODE = current_test_mode)
      if (preinstalled == "TRUE") {
        Sys.setenv(CASSIA_PREINSTALLED = "TRUE")
      }
      exit_code <- system2("Rscript", args = shQuote(test_script),
                           stdout = stdout_file, stderr = stderr_file, wait = TRUE)
    }
    duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # Read captured output (last 2000 chars for stderr)
    stdout_content <- ""
    stderr_content <- ""
    if (file.exists(stdout_file)) {
      stdout_content <- paste(readLines(stdout_file, warn = FALSE), collapse = "\n")
      unlink(stdout_file)
    }
    if (file.exists(stderr_file)) {
      stderr_content <- paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
      # Keep last 2000 chars
      if (nchar(stderr_content) > 2000) {
        stderr_content <- substr(stderr_content, nchar(stderr_content) - 1999, nchar(stderr_content))
      }
      unlink(stderr_file)
    }

    res <- list(
      name = basename(test_folder),
      status = if (exit_code == 0) "passed" else "failed",
      returncode = exit_code,
      duration = round(duration, 2),
      mode = if (install_mode) "install" else "load_all",
      script = basename(test_script)
    )

    # Add error info for failed tests
    if (exit_code != 0) {
      res$stderr <- stderr_content
      res$reason <- if (nchar(stderr_content) > 0) {
        # Extract last meaningful line as reason
        lines <- strsplit(stderr_content, "\n")[[1]]
        lines <- lines[nchar(trimws(lines)) > 0]
        if (length(lines) > 0) tail(lines, 1) else paste("Exit code:", exit_code)
      } else {
        paste("Exit code:", exit_code)
      }
    }

    res
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

save_report <- function(results, total_duration, mode = "development") {
  # Create organized folder structure: reports/r/{mode}/{timestamp}/
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  reports_dir <- file.path(script_dir, "reports", "r", mode, timestamp)
  if (!dir.exists(reports_dir)) {
    dir.create(reports_dir, recursive = TRUE)
  }

  report_path <- file.path(reports_dir, "report.json")

  report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    mode = mode,
    language = "R",
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
  cat("\nJSON report saved to:", report_path, "\n")

  # Return the reports_dir so HTML report can use the same folder
  return(reports_dir)
}

save_html_report <- function(results, total_duration, mode = "development", reports_dir = NULL) {
  # If reports_dir not provided, create organized folder structure
  if (is.null(reports_dir)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    reports_dir <- file.path(script_dir, "reports", "r", mode, timestamp)
    if (!dir.exists(reports_dir)) {
      dir.create(reports_dir, recursive = TRUE)
    }
  }

  html_path <- file.path(reports_dir, "report.html")

  passed <- sum(sapply(results, function(r) r$status == "passed"))
  failed <- sum(sapply(results, function(r) r$status == "failed"))
  errors <- sum(sapply(results, function(r) r$status %in% c("error", "timeout")))
  skipped <- sum(sapply(results, function(r) r$status == "skipped"))
  total <- length(results)

  # Determine overall status color
  if (failed + errors == 0) {
    status_color <- "#28a745"  # green
    status_text <- "ALL TESTS PASSED"
  } else {
    status_color <- "#dc3545"  # red
    status_text <- paste(failed + errors, "TEST(S) FAILED")
  }

  # Build test rows with error details
  test_rows <- sapply(results, function(r) {
    status_class <- switch(r$status,
      "passed" = "passed",
      "failed" = "failed",
      "error" = "error",
      "timeout" = "error",
      "skipped" = "skipped",
      "unknown"
    )
    status_icon <- switch(r$status,
      "passed" = "&#10004;",
      "failed" = "&#10008;",
      "error" = "&#9888;",
      "timeout" = "&#8986;",
      "skipped" = "&#8212;",
      "?"
    )
    duration_str <- if (r$duration > 0) paste0(round(r$duration, 1), "s") else "N/A"

    # Build the main row
    row_html <- sprintf('<tr class="%s"><td>%s</td><td>%s %s</td><td>%s</td></tr>',
            status_class, r$name, status_icon, toupper(r$status), duration_str)

    # Add error details row for failed/error tests
    if (r$status %in% c("failed", "error", "timeout")) {
      error_content <- ""
      if (!is.null(r$stderr) && nchar(r$stderr) > 0) {
        # Escape HTML characters
        escaped_stderr <- gsub("&", "&amp;", r$stderr)
        escaped_stderr <- gsub("<", "&lt;", escaped_stderr)
        escaped_stderr <- gsub(">", "&gt;", escaped_stderr)
        error_content <- escaped_stderr
      } else if (!is.null(r$reason) && nchar(r$reason) > 0) {
        escaped_reason <- gsub("&", "&amp;", r$reason)
        escaped_reason <- gsub("<", "&lt;", escaped_reason)
        escaped_reason <- gsub(">", "&gt;", escaped_reason)
        error_content <- escaped_reason
      }

      if (nchar(error_content) > 0) {
        row_html <- paste0(row_html, sprintf(
          '<tr class="error-details-row"><td colspan="3"><pre class="error-details">%s</pre></td></tr>',
          error_content
        ))
      }
    }

    row_html
  })

  html_content <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CASSIA Test Report</title>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f5f5f5; padding: 20px; }
    .container { max-width: 900px; margin: 0 auto; background: white; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); overflow: hidden; }
    .header { background: %s; color: white; padding: 30px; text-align: center; }
    .header h1 { font-size: 24px; margin-bottom: 10px; }
    .header .status { font-size: 18px; font-weight: bold; }
    .meta { padding: 20px; background: #f8f9fa; border-bottom: 1px solid #dee2e6; display: flex; justify-content: space-around; flex-wrap: wrap; }
    .meta-item { text-align: center; padding: 10px; }
    .meta-item .label { font-size: 12px; color: #6c757d; text-transform: uppercase; }
    .meta-item .value { font-size: 20px; font-weight: bold; color: #333; }
    .summary { display: flex; justify-content: center; gap: 20px; padding: 20px; flex-wrap: wrap; }
    .summary-item { padding: 15px 25px; border-radius: 6px; text-align: center; min-width: 100px; }
    .summary-item.passed { background: #d4edda; color: #155724; }
    .summary-item.failed { background: #f8d7da; color: #721c24; }
    .summary-item.errors { background: #fff3cd; color: #856404; }
    .summary-item.skipped { background: #e2e3e5; color: #383d41; }
    .summary-item .count { font-size: 28px; font-weight: bold; }
    .summary-item .label { font-size: 12px; text-transform: uppercase; }
    table { width: 100%%; border-collapse: collapse; }
    th { background: #343a40; color: white; padding: 12px 15px; text-align: left; }
    td { padding: 12px 15px; border-bottom: 1px solid #dee2e6; }
    tr:hover { background: #f8f9fa; }
    tr.passed td:nth-child(2) { color: #28a745; }
    tr.failed td:nth-child(2) { color: #dc3545; }
    tr.error td:nth-child(2) { color: #ffc107; }
    tr.skipped td:nth-child(2) { color: #6c757d; }
    tr.error-details-row { background: #fff8f8; }
    tr.error-details-row:hover { background: #fff8f8; }
    .error-details { background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; font-size: 12px; white-space: pre-wrap; word-wrap: break-word; color: #721c24; max-height: 300px; overflow-y: auto; margin: 0; border-left: 3px solid #dc3545; }
    .footer { padding: 15px; text-align: center; color: #6c757d; font-size: 12px; background: #f8f9fa; }
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <h1>CASSIA Test Report</h1>
      <div class="status">%s</div>
    </div>
    <div class="meta">
      <div class="meta-item">
        <div class="label">Mode</div>
        <div class="value">%s</div>
      </div>
      <div class="meta-item">
        <div class="label">Duration</div>
        <div class="value">%.1fs</div>
      </div>
      <div class="meta-item">
        <div class="label">Timestamp</div>
        <div class="value">%s</div>
      </div>
    </div>
    <div class="summary">
      <div class="summary-item passed"><div class="count">%d</div><div class="label">Passed</div></div>
      <div class="summary-item failed"><div class="count">%d</div><div class="label">Failed</div></div>
      <div class="summary-item errors"><div class="count">%d</div><div class="label">Errors</div></div>
      <div class="summary-item skipped"><div class="count">%d</div><div class="label">Skipped</div></div>
    </div>
    <table>
      <thead><tr><th>Test Name</th><th>Status</th><th>Duration</th></tr></thead>
      <tbody>%s</tbody>
    </table>
    <div class="footer">Generated by CASSIA Test Suite</div>
  </div>
</body>
</html>',
    status_color, status_text, mode, total_duration,
    format(Sys.time(), "%%Y-%%m-%%d %%H:%%M:%%S"),
    passed, failed, errors, skipped,
    paste(test_rows, collapse = "\n")
  )

  writeLines(html_content, html_path)
  cat("HTML report saved to:", html_path, "\n")

  # Try to open in browser
  tryCatch({
    browseURL(paste0("file://", normalizePath(html_path, winslash = "/")))
  }, error = function(e) {
    # Silently ignore if browser can't be opened
  })
}

parse_test_list <- function(test_string) {
  if (is.null(test_string) || test_string == "") {
    return(character(0))
  }
  trimws(strsplit(test_string, ",")[[1]])
}

# Parse arguments
skip_tests <- c("13", "14", "16")  # Skip tests 13, 14, 16 by default
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

# Set environment variable for test mode (so individual tests can read it)
test_mode <- if (use_install_mode) "installed" else "development"
Sys.setenv(CASSIA_TEST_MODE = test_mode)

# Main
cat(strrep("=", 60), "\n")
cat("CASSIA TEST SUITE (R)\n")
cat(strrep("=", 60), "\n")
mode_str <- if (use_install_mode) "install_github (full package)" else "load_all (fast dev)"
cat("Mode:", mode_str, "\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Pre-install CASSIA once if in install mode (so individual tests don't reinstall)
if (use_install_mode) {
  cat("\nPre-installing CASSIA package (one-time)...\n")
  source(file.path(script_dir, "shared", "r", "test_utils.R"))
  tryCatch({
    setup_cassia_install(force = TRUE)
    cat("CASSIA package installed successfully.\n")
    # Set env var so individual tests know CASSIA is already installed
    Sys.setenv(CASSIA_PREINSTALLED = "TRUE")
  }, error = function(e) {
    cat("Warning: Pre-installation failed:", e$message, "\n")
    cat("Individual tests will attempt installation.\n")
  })
}

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
  reports_dir <- save_report(results, total_duration, test_mode)
  save_html_report(results, total_duration, test_mode, reports_dir)
}

quit(status = if (all_passed) 0 else 1)
