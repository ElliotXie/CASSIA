#' Merge Cell Type Annotations using LLM
#'
#' This function uses large language models to merge and group single-cell RNA-seq cluster annotations
#' at different levels of granularity. Results are saved directly to a CSV file.
#'
#' @param csv_path Path to the CSV file containing cluster annotations
#' @param output_path Path to save the results (if NULL, saves back to input CSV file)
#' @param provider LLM provider to use ("openai", "anthropic", or "openrouter")
#' @param model Specific model to use (if NULL, uses default for provider)
#' @param api_key API key for the provider (if NULL, gets from environment)
#' @param additional_context Optional domain-specific context to help with annotation
#' @param batch_size Number of clusters to process in each LLM call
#' @param detail_level Level of detail for the groupings: "broad", "detailed", or "very_detailed"
#' @param process_all If TRUE, processes all detail levels sequentially
#' @param debug If TRUE, print additional debugging information
#'
#' @return Invisibly returns TRUE if successful
#'
#' @importFrom reticulate import source_python py_set_seed use_python py_last_error py_has_attr
#' @export
runCASSIA_merge_annotations <- function(csv_path,
                                        output_path = NULL,
                                        provider = "openrouter",
                                        model = "deepseek/deepseek-chat-v3-0324",
                                        additional_context = NULL,
                                        batch_size = 20,
                                        detail_level = "broad",
                                        process_all = FALSE, # Deprecated
                                        debug = FALSE) {

  # Handle deprecated 'process_all' argument
  if (process_all) {
    warning("'process_all' is deprecated. Please use detail_level = 'all' instead.")
    detail_level <- "all"
  }

  tryCatch({
    # Verify that the input file exists
    if (!file.exists(csv_path)) {
      stop("Input CSV file does not exist: ", csv_path)
    }

    # Use the same output path as input if not provided
    if (is.null(output_path)) {
      output_path <- csv_path
    }

    # Check if the output file is accessible (basic check on Windows)
    tryCatch({
      test_con <- file(output_path, "a") # Open for appending
      close(test_con)
    }, error = function(e) {
      stop("The output file appears to be open or inaccessible: ", output_path,
           ". Please close it in other programs or check permissions.")
    })


    # Use the already imported Python module
    if (is.null(py_merging)) {
      stop("Python merging module not loaded. Please restart R and load the CASSIA package.")
    }
    merging_module <- py_merging
    message("Using pre-loaded merging_annotation module.")

    # Read the data once
    df <- read.csv(csv_path)

    # Determine which levels to process
    levels_to_process <- if (detail_level == "all") {
      c("broad", "detailed", "very_detailed")
    } else {
      c(detail_level)
    }

    message(paste("Starting annotation merging for detail level(s):", paste(levels_to_process, collapse=", ")))

    # Process each level sequentially
    for (level in levels_to_process) {
      message(paste0("Processing '", level, "' level groupings..."))

      # Call the Python function, which returns a pandas DataFrame
      result_df_py <- merging_module$merge_annotations(
        csv_path = csv_path, # Pass path for initial read inside python
        output_path = NULL, # Do not save from Python
        provider = provider,
        model = model,
        api_key = NULL, # Handled by llm_utils
        additional_context = additional_context,
        batch_size = as.integer(batch_size),
        detail_level = level
      )

      # Convert pandas DataFrame to R data.frame
      result_df_r <- reticulate::py_to_r(result_df_py)

      # Determine the column name for the current level
      result_column <- switch(level,
                              "broad" = "Merged_Grouping_1",
                              "detailed" = "Merged_Grouping_2",
                              "very_detailed" = "Merged_Grouping_3")

      # Add or update the column in the R dataframe
      if (!is.null(result_column) && result_column %in% names(result_df_r)) {
        df[[result_column]] <- result_df_r[[result_column]]
        message(paste0("Completed and updated '", level, "' level groupings."))
      } else {
        warning(paste("Could not find expected result column for level:", level))
      }
    }

    # Save the final, combined data.frame to the output file
    tryCatch({
      write.csv(df, file = output_path, row.names = FALSE)
      message("Annotation merging completed. Results saved to: ", output_path)
    }, error = function(e) {
      stop("Failed to save the final CSV file: ", e$message)
    })

    # Return TRUE invisibly on success
    invisible(TRUE)

  }, error = function(e) {
    message("Error in merging annotations: ", e$message)
    if (reticulate::py_has_attr(reticulate::py, "last_error") && !is.null(reticulate::py_last_error())) {
      traceback <- reticulate::py_last_error()
      message("Python traceback: \n", traceback$traceback)
    }
    stop("Failed to run annotation merging.")
  })
}

