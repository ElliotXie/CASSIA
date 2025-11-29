# CASSIA Test Runner - Run Individual Test (R)
# =============================================
# Run a specific R test by number or name.
#
# Usage:
#     Rscript run_test.R 01              # Run with load_all() (default, fast)
#     Rscript run_test.R 01 --install    # Run with install_local() (full package)
#     Rscript run_test.R batch
#     Rscript run_test.R 03_validator

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

# Parse --install flag
use_install_mode <- "--install" %in% args
args <- args[args != "--install"]

find_test_folder <- function(test_id) {
  # Get all test folders (those starting with digits)
  all_items <- list.dirs(script_dir, recursive = FALSE, full.names = TRUE)
  test_folders <- all_items[grepl("^[0-9]", basename(all_items))]
  test_folders <- sort(test_folders)

  # Try exact match first
  for (folder in test_folders) {
    if (basename(folder) == test_id) {
      return(folder)
    }
  }

  # Try matching by number prefix
  for (folder in test_folders) {
    if (startsWith(basename(folder), test_id) ||
        startsWith(basename(folder), paste0(test_id, "_"))) {
      return(folder)
    }
  }

  # Try matching by name substring
  for (folder in test_folders) {
    if (grepl(tolower(test_id), tolower(basename(folder)))) {
      return(folder)
    }
  }

  return(NULL)
}

run_r_test <- function(test_folder, install_mode = FALSE) {
  # Find the R test script based on mode
  if (install_mode) {
    # Look for *_install.R files
    test_scripts <- list.files(test_folder, pattern = "^test_.*_install\\.R$", full.names = TRUE)
    mode_label <- "(install mode)"
  } else {
    # Look for regular test_*.R files (excluding *_install.R)
    all_scripts <- list.files(test_folder, pattern = "^test_.*\\.R$", full.names = TRUE)
    test_scripts <- all_scripts[!grepl("_install\\.R$", all_scripts)]
    mode_label <- "(load_all mode)"
  }

  if (length(test_scripts) == 0) {
    cat("No R test script found in", test_folder, mode_label, "\n")
    return(FALSE)
  }

  test_script <- test_scripts[1]
  cat("\nRunning:", basename(test_script), mode_label, "\n")
  cat(strrep("-", 50), "\n")

  # Run the test in a separate R process (allows proper path detection)
  # Use stdout="" and stderr="" to show output in real-time
  result <- tryCatch({
    exit_code <- system2("Rscript", test_script, wait = TRUE,
                         stdout = "", stderr = "")
    exit_code == 0
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    FALSE
  })

  return(result)
}

list_available_tests <- function() {
  all_items <- list.dirs(script_dir, recursive = FALSE, full.names = TRUE)
  test_folders <- all_items[grepl("^[0-9]", basename(all_items))]
  test_folders <- sort(test_folders)

  cat("\nAvailable Tests:\n")
  cat(strrep("-", 40), "\n")
  for (folder in test_folders) {
    cat("  ", basename(folder), "\n")
  }
  cat("\n")
}

# Main
if (length(args) == 0) {
  cat("CASSIA Test Runner - Individual Test (R)\n")
  cat(strrep("=", 40), "\n")
  cat("\nUsage:\n")
  cat("  Rscript run_test.R <test_id>             # Run with load_all() (fast, default)\n")
  cat("  Rscript run_test.R <test_id> --install   # Run with install_local() (full package)\n")
  cat("\nExamples:\n")
  cat("  Rscript run_test.R 01\n")
  cat("  Rscript run_test.R 01 --install\n")
  cat("  Rscript run_test.R batch\n")
  cat("  Rscript run_test.R 03_validator\n")
  list_available_tests()
  quit(status = 1)
}

test_id <- args[1]

# Handle special cases
if (test_id %in% c("--list", "-l", "list")) {
  list_available_tests()
  quit(status = 0)
}

# Find and run the test
test_folder <- find_test_folder(test_id)

if (is.null(test_folder)) {
  cat("Error: Test '", test_id, "' not found\n", sep = "")
  list_available_tests()
  quit(status = 1)
}

cat("\n", strrep("=", 50), "\n")
mode_str <- if (use_install_mode) "(install mode)" else "(load_all mode)"
cat("Running Test:", basename(test_folder), mode_str, "\n")
cat(strrep("=", 50), "\n")

success <- run_r_test(test_folder, install_mode = use_install_mode)

cat("\n", strrep("=", 50), "\n")
if (success) {
  cat("Test PASSED\n")
} else {
  cat("Test FAILED\n")
}
cat(strrep("=", 50), "\n")

quit(status = if (success) 0 else 1)
