# CASSIA Test Framework - Logging Manager
# ========================================
# Provides timestamped console logging with file output.
# All log messages have format: [YYYY-MM-DD HH:MM:SS] message

# Global state for logging (environment to avoid global variable issues)
.log_state <- new.env(parent = emptyenv())
.log_state$log_file <- NULL
.log_state$log_con <- NULL

#' Timestamped log function (replaces cat())
#'
#' Writes timestamped message to both console and log file (if active).
#'
#' @param ... Message components to paste together
#' @param sep Separator between components (default: " ")
#' @examples
#' log_msg("Running annotation for cluster:", cluster_name)
#' log_msg("Result:", result$main_cell_type)
log_msg <- function(..., sep = " ") {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  msg <- paste(..., sep = sep)
  full_msg <- paste(timestamp, msg)

  # Write to console
  cat(full_msg, "\n")

  # Write to log file if active
  if (!is.null(.log_state$log_con)) {
    tryCatch({
      writeLines(full_msg, .log_state$log_con)
      flush(.log_state$log_con)
    }, error = function(e) {
      # Silently ignore write errors to avoid disrupting test flow
    })
  }
}

#' Start logging to file
#'
#' Opens a log file in the logs directory and begins capturing output.
#'
#' @param logs_dir Path to the logs directory (e.g., results$logs)
#' @return The path to the log file (invisible)
start_logging <- function(logs_dir) {
  log_path <- file.path(logs_dir, "console_output.log")
  .log_state$log_file <- log_path

  tryCatch({
    .log_state$log_con <- file(log_path, open = "wt")
    log_msg("Logging started")
    log_msg("Log file:", log_path)
  }, error = function(e) {
    warning("Could not open log file: ", e$message)
    .log_state$log_con <- NULL
  })

  invisible(log_path)
}

#' Stop logging
#'
#' Closes the log file connection and stops file logging.
stop_logging <- function() {
  if (!is.null(.log_state$log_con)) {
    tryCatch({
      log_msg("Logging stopped")
      close(.log_state$log_con)
    }, error = function(e) {
      # Silently ignore close errors
    })
    .log_state$log_con <- NULL
  }
}

#' Log error with traceback
#'
#' Logs an error message along with the R traceback for debugging.
#'
#' @param e An error condition object (from tryCatch)
log_error <- function(e) {
  log_msg("ERROR:", e$message)

  # Capture and log traceback if available
  tb <- tryCatch({
    capture.output(traceback())
  }, error = function(err) {
    character(0)
  })

  if (length(tb) > 0 && !all(tb == "No traceback available")) {
    log_msg("Traceback:")
    for (line in tb) {
      if (nchar(trimws(line)) > 0) {
        log_msg("  ", line)
      }
    }
  }
}

#' Log a separator line
#'
#' Logs a line of repeated characters for visual separation.
#'
#' @param char Character to repeat (default: "=")
#' @param width Number of repetitions (default: 50)
log_separator <- function(char = "=", width = 50) {
  log_msg(strrep(char, width))
}

#' Check if logging is active
#'
#' @return TRUE if file logging is active, FALSE otherwise
is_logging_active <- function() {
  !is.null(.log_state$log_con)
}

#' Get current log file path
#'
#' @return Path to current log file, or NULL if not logging
get_log_file <- function() {
  .log_state$log_file
}
