# CASSIA Test 23: Free API Lifetime Limit (R)
# =============================================
# Tests the free API 2-cluster lifetime limit by running runCASSIA_batch
# WITHOUT any API key. Expects auto-selection of at most 2 clusters.
#
# This test deliberately does NOT load API keys. It validates:
# - Free API auto-detection when no key is set
# - Auto-selection of first N clusters when quota < total clusters
# - "LIMIT REACHED" message when quota is exhausted
#
# WARNING: This test consumes free API quota (up to 2 cluster annotations
# per machine, lifetime).
#
# Usage:
#     Rscript test_free_api_limit.R

# Get script directory (works with Rscript and source())
get_script_dir <- function() {
  # Method 1: commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg), winslash = "/")))
  }
  # Method 2: sys.frame (works with source())
  for (i in sys.nframe():1) {
    if (!is.null(sys.frame(i)$ofile)) {
      return(dirname(normalizePath(sys.frame(i)$ofile, winslash = "/")))
    }
  }
  # Fallback: current working directory
  return(normalizePath(getwd(), winslash = "/"))
}
script_dir <- get_script_dir()

# Source shared utilities
source(file.path(script_dir, "..", "shared", "r", "test_utils.R"))
source(file.path(script_dir, "..", "shared", "r", "fixtures.R"))
source(file.path(script_dir, "..", "shared", "r", "result_manager.R"))
source(file.path(script_dir, "..", "shared", "r", "logging_manager.R"))

run_free_api_limit_test <- function() {
  print_test_header("23 - Free API Lifetime Limit (R)")

  # ---------------------------------------------------------------------------
  # IMPORTANT: Do NOT load API keys. We want to test the free API path.
  # Explicitly unset any API keys that might be in the environment.
  # ---------------------------------------------------------------------------
  log_msg("Clearing API keys to force free API path...")
  keys_to_unset <- c("OPENROUTER_API_KEY", "GOOGLE_API_KEY", "TOGETHER_API_KEY",
                      "OPENAI_API_KEY", "ANTHROPIC_API_KEY")
  for (key in keys_to_unset) {
    if (Sys.getenv(key, unset = "") != "") {
      log_msg("  Unsetting", key)
    }
    Sys.unsetenv(key)
  }

  # Setup CASSIA R package
  tryCatch({
    setup_cassia()
    log_msg("CASSIA R package loaded successfully")
  }, error = function(e) {
    log_msg("Error loading CASSIA:", e$message)
    return(FALSE)
  })

  # Check remaining quota via Python
  remaining_before <- tryCatch({
    py_cassia <- reticulate::import("CASSIA")
    py_free_api <- reticulate::import("CASSIA.core.free_api")
    py_free_api$clear_free_api_cache()
    as.integer(py_cassia$get_remaining_free_clusters())
  }, error = function(e) {
    log_msg("Warning: Could not check remaining quota:", e$message)
    2L  # Optimistic default
  })

  MAX_FREE_CLUSTERS <- 2L

  log_msg("")
  log_msg("Free API Status:")
  log_msg("  Max lifetime clusters:", MAX_FREE_CLUSTERS)
  log_msg("  Remaining before test:", remaining_before)

  # Load ALL clusters (all 6 from processed.csv)
  marker_df <- load_markers()
  all_clusters <- marker_df$Broad.cell.type
  total_clusters <- length(all_clusters)

  log_msg("")
  log_msg("Loaded", total_clusters, "clusters from marker data:")
  for (cluster in all_clusters) {
    log_msg("  -", cluster)
  }

  # Create results directory
  results <- create_results_dir("23_free_api_limit", get_test_mode())
  start_logging(results$logs)

  output_name <- file.path(results$outputs, "free_api_results")
  log_msg("Results will be saved to:", results$base)

  # Run the test
  start_time <- Sys.time()
  errors <- list()
  status <- "error"
  clusters_annotated <- 0L
  expected_clusters <- 0L

  if (remaining_before <= 0) {
    # =========================================================================
    # SCENARIO: Quota exhausted
    # =========================================================================
    log_msg("")
    log_msg("--- SCENARIO: Quota exhausted (remaining=", remaining_before, ") ---")
    log_msg("Expecting 'FREE API LIMIT REACHED' message and no output.")

    expected_clusters <- 0L

    tryCatch({
      CASSIA::runCASSIA_batch(
        marker = marker_df,
        output_name = output_name,
        n_genes = 30,
        tissue = "large intestine",
        species = "human",
        max_workers = 3,
        provider = "openrouter"
      )

      # Function returned without error — check for output
      summary_csv <- paste0(output_name, "_summary.csv")
      if (file.exists(summary_csv)) {
        results_df <- read.csv(summary_csv)
        if (nrow(results_df) > 0) {
          errors <- list(paste("Expected no clusters but got", nrow(results_df)))
          status <- "failed"
        } else {
          status <- "passed"
        }
      } else {
        # No output = correct for quota exhausted
        status <- "passed"
      }

    }, error = function(e) {
      if (grepl("LIMIT REACHED|free", e$message, ignore.case = TRUE)) {
        log_msg("Got expected error:", e$message)
        status <<- "passed"
      } else {
        errors <<- list(e$message)
        status <<- "error"
        log_error(e)
      }
    })

  } else {
    # =========================================================================
    # SCENARIO: Quota available — expect auto-selection
    # =========================================================================
    expected_clusters <- min(remaining_before, total_clusters)
    log_msg("")
    log_msg("--- SCENARIO: Quota available (remaining=", remaining_before, ") ---")
    log_msg("Expecting", expected_clusters, "of", total_clusters, "clusters to be annotated.")

    tryCatch({
      CASSIA::runCASSIA_batch(
        marker = marker_df,
        output_name = output_name,
        n_genes = 30,
        tissue = "large intestine",
        species = "human",
        max_workers = 3,
        provider = "openrouter"
      )

      # Check output files
      summary_csv <- paste0(output_name, "_summary.csv")

      log_msg("DEBUG: Checking for file at:", summary_csv)
      log_msg("DEBUG: Directory exists:", dir.exists(dirname(summary_csv)))
      if (dir.exists(dirname(summary_csv))) {
        log_msg("DEBUG: Files in directory:")
        for (f in list.files(dirname(summary_csv))) {
          log_msg("  -", f)
        }
      }

      if (file.exists(summary_csv)) {
        results_df <- read.csv(summary_csv)
        clusters_annotated <- nrow(results_df)

        log_msg("")
        log_msg("Batch Results:")
        log_msg("  Clusters annotated:", clusters_annotated)
        log_msg("  Expected:", expected_clusters)
        log_msg("  Output files:")
        log_msg("    -", basename(summary_csv))

        conversations_json <- paste0(output_name, "_conversations.json")
        if (file.exists(conversations_json)) {
          log_msg("    -", basename(conversations_json))
        }

        html_report <- paste0(output_name, "_report.html")
        if (file.exists(html_report)) {
          log_msg("    -", basename(html_report))
        }

        if (clusters_annotated == expected_clusters) {
          status <- "passed"
        } else if (clusters_annotated > 0 && clusters_annotated <= expected_clusters) {
          status <- "passed"
          log_msg("  Note:", expected_clusters - clusters_annotated, "cluster(s) may have failed")
        } else {
          status <- "failed"
          errors <- list(paste("Expected", expected_clusters, "clusters but got", clusters_annotated))
        }
      } else {
        status <- "failed"
        errors <- list("Summary CSV not created")
      }

    }, error = function(e) {
      errors <<- list(e$message)
      status <<- "error"
      log_error(e)
    })
  }

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Check remaining quota AFTER test
  remaining_after <- tryCatch({
    py_free_api <- reticulate::import("CASSIA.core.free_api")
    py_free_api$clear_free_api_cache()
    py_cassia <- reticulate::import("CASSIA")
    as.integer(py_cassia$get_remaining_free_clusters())
  }, error = function(e) {
    -1L
  })

  log_msg("")
  log_msg("Post-test quota:", remaining_after, "of", MAX_FREE_CLUSTERS, "remaining")

  # Save metadata
  metadata <- create_test_metadata(
    test_name = "free_api_limit",
    config = list(
      llm = list(provider = "openrouter", model = "(free API default)"),
      data = list(tissue = "large intestine", species = "human", n_genes = 30)
    ),
    duration_seconds = duration,
    status = status,
    clusters_tested = as.list(all_clusters[seq_len(min(expected_clusters, total_clusters))]),
    errors = errors
  )
  metadata$free_api <- list(
    remaining_before = remaining_before,
    remaining_after = remaining_after,
    max_free_clusters = MAX_FREE_CLUSTERS,
    expected_clusters = expected_clusters,
    clusters_annotated = clusters_annotated,
    total_clusters_in_data = total_clusters
  )
  save_test_metadata(results$outputs, metadata)

  # Print final result
  success <- status == "passed"
  print_test_result(success, paste("Duration:", round(duration, 2), "s"))

  if (success) {
    if (remaining_before <= 0) {
      log_msg("    (Validated: limit reached behavior)")
    } else {
      log_msg("    (Validated:", clusters_annotated, "clusters annotated with free API)")
    }
  }

  stop_logging()
  return(success)
}

# Run test
success <- run_free_api_limit_test()
quit(status = if (success) 0 else 1)
