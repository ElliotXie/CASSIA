# install_cassia_local.R
# Script to install the CASSIA R package locally for testing purposes

cat("=== CASSIA Local Package Installation ===\n\n")

# Check and install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  cat("Installing devtools...\n")
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# Get the directory where this script is located
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  getwd()
})

# Path to CASSIA_R package (relative to Test folder)
cassia_r_path <- file.path(script_dir, "..", "CASSIA_R")

# Normalize the path
cassia_r_path <- normalizePath(cassia_r_path, mustWork = FALSE)

cat("Package source path:", cassia_r_path, "\n\n")

# Check if the package directory exists
if (!dir.exists(cassia_r_path)) {
  stop("CASSIA_R package directory not found at: ", cassia_r_path)
}

# Check for DESCRIPTION file
if (!file.exists(file.path(cassia_r_path, "DESCRIPTION"))) {
  stop("DESCRIPTION file not found in CASSIA_R directory")
}

cat("Installing CASSIA package from local source...\n")
cat("This may take a few minutes...\n\n")

# Install the package locally with force to allow reinstallation
tryCatch({
  devtools::install_local(
    cassia_r_path,
    force = TRUE,
    upgrade = "never",
    quiet = FALSE
  )

  cat("\n=== Installation Complete ===\n\n")

  # Verify installation by loading the package
  cat("Verifying installation...\n")
  library(CASSIA)
  cat("CASSIA package loaded successfully!\n")
  cat("Package version:", as.character(packageVersion("CASSIA")), "\n")

}, error = function(e) {
  cat("\n=== Installation Failed ===\n")
  cat("Error:", conditionMessage(e), "\n")
  stop(e)
})
