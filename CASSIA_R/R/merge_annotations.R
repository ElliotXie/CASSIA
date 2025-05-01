#' Merge Cell Type Annotations using LLM
#'
#' This function uses large language models to merge and group single-cell RNA-seq cluster annotations
#' at different levels of granularity.
#'
#' @param csv_path Path to the CSV file containing cluster annotations
#' @param output_path Path to save the results (if NULL, returns data frame without saving)
#' @param provider LLM provider to use ("openai", "anthropic", or "openrouter")
#' @param model Specific model to use (if NULL, uses default for provider)
#' @param api_key API key for the provider (if NULL, gets from environment)
#' @param additional_context Optional domain-specific context to help with annotation
#' @param batch_size Number of clusters to process in each LLM call
#' @param detail_level Level of detail for the groupings: "broad", "detailed", or "very_detailed"
#' @param process_all If TRUE, processes all detail levels sequentially
#'
#' @return A data frame with the original annotations and suggested cell groupings (if output_path is NULL),
#'         otherwise saves to file and returns invisibly
#'
#' @importFrom reticulate import source_python py_set_seed use_python
#' @export
runCASSIA_merge_annotations <- function(csv_path,
                             output_path = NULL,
                             provider = "openrouter",
                             model = NULL,
                             additional_context = NULL,
                             batch_size = 100,
                             detail_level = "broad",
                             process_all = FALSE) {

  tryCatch({
    # Check if py_merging is loaded
    if (is.null(py_merging)) {
      # Try to load it if not already loaded
      py_merging <<- reticulate::import_from_path("merging_annotation", path = system.file("python", package = "CASSIA"))
    }

    # Prepare arguments
    args <- list(
      csv_path = csv_path,
      provider = provider,
      model = model,
      additional_context = additional_context,
      batch_size = as.integer(batch_size)
    )
    
    # If output path is provided, add it to args
    if (!is.null(output_path)) {
      args$output_path <- output_path
    }
    
    # Choose which Python function to call based on process_all parameter
    if (process_all) {
      # Process each detail level sequentially
      
      # First, process broad level (Merged_Grouping_1) with a temporary output
      temp_output1 <- tempfile(fileext = ".csv")
      args_temp <- args
      args_temp$output_path <- temp_output1
      args_temp$detail_level <- "broad"
      message("Processing broad level groupings...")
      do.call(py_merging$merge_annotations, args_temp)
      
      # Now process detailed level (Merged_Grouping_2) with a temporary output
      temp_output2 <- tempfile(fileext = ".csv")
      args_temp$output_path <- temp_output2
      args_temp$detail_level <- "detailed"
      message("Processing detailed level groupings...")
      do.call(py_merging$merge_annotations, args_temp)
      
      # Finally process very_detailed level (Merged_Grouping_3) with a temporary output
      temp_output3 <- tempfile(fileext = ".csv")
      args_temp$output_path <- temp_output3
      args_temp$detail_level <- "very_detailed"
      message("Processing very detailed level groupings...")
      do.call(py_merging$merge_annotations, args_temp)
      
      # Read all the temporary files
      df1 <- read.csv(temp_output1)
      df2 <- read.csv(temp_output2)
      df3 <- read.csv(temp_output3)
      
      # Combine the results into a single dataframe
      # Start with df3 (very_detailed) as the base
      combined_df <- df3
      
      # Add columns from df1 and df2
      if ("Merged_Grouping_1" %in% colnames(df1)) {
        combined_df$Merged_Grouping_1 <- df1$Merged_Grouping_1
      }
      if ("Merged_Grouping_2" %in% colnames(df2)) {
        combined_df$Merged_Grouping_2 <- df2$Merged_Grouping_2
      }
      
      # Save the combined dataframe to the output path
      if (!is.null(output_path)) {
        write.csv(combined_df, output_path, row.names = FALSE)
        message("Saved combined results to ", output_path)
      }
      
      # Clean up temporary files
      file.remove(temp_output1, temp_output2, temp_output3)
      
      message("All grouping levels processed and combined successfully")
      
    } else {
      # Call the function that processes a single detail level
      args$detail_level <- detail_level
      message(paste0("Processing ", detail_level, " level groupings..."))
      result <- do.call(py_merging$merge_annotations, args)
      message("Processing completed successfully")
    }
    
  }, error = function(e) {
    message("Error importing Python module: ", e$message)
    message("Python traceback: ", reticulate::py_last_error())
    stop("Failed to load Python module. Make sure Python dependencies are installed.")
  })
}

