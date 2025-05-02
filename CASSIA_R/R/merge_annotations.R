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
                             model = "deepseek/deepseek-chat-v3-0324",
                             additional_context = NULL,
                             batch_size = 20,
                             detail_level = "broad",
                             process_all = FALSE) {

  tryCatch({
    # Check if reticulate is available
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("Package 'reticulate' is required. Please install it using install.packages('reticulate')")
    }
    
    # Get the conda environment name from the package configuration
    conda_env <- getOption("CASSIA.conda_env", default = "cassia_env")
    
    # Use the package's conda environment
    tryCatch({
      # Check if the environment exists, if not create it
      if (!conda_env %in% reticulate::conda_list()$name) {
        message("Creating conda environment: ", conda_env)
        setup_cassia_env(conda_env)
      }
      
      reticulate::use_condaenv(conda_env, required = TRUE)
    }, error = function(e) {
      message("Could not load conda environment '", conda_env, "'. Using default Python.")
    })
    
    # Make sure the Python module is loaded
    if (is.null(py_merging)) {
      message("Loading Python module...")
      py_merging <<- reticulate::import_from_path("merging_annotation", path = system.file("python", package = "CASSIA"))
    }

    # Prepare arguments
    args <- list(
      csv_path = csv_path,
      provider = provider,
      model = model,
      additional_context = additional_context,
      batch_size = as.integer(batch_size)  # Ensure this is an integer
    )
    
    # If output path is provided, add it to args
    if (!is.null(output_path)) {
      args$output_path <- output_path
    }
    
    # Choose which Python function to call based on process_all parameter
    if (process_all) {
      # Process all levels sequentially using the Python function
      message("Processing all detail levels...")
      result <- do.call(py_merging$merge_annotations_all, args)
      message("All detail levels processed successfully")
    } else {
      # Call the function that processes a single detail level
      args$detail_level <- detail_level
      message(paste0("Processing ", detail_level, " level groupings..."))
      result <- do.call(py_merging$merge_annotations, args)
      message("Processing completed successfully")
    }
    
    # Convert result to R dataframe if it's a Python object
    if (inherits(result, "python.builtin.object")) {
      result <- reticulate::py_to_r(result)
    }
    
    # Return the result invisibly if output_path was provided
    if (is.null(output_path)) {
      return(result)
    } else {
      invisible(result)
    }
    
  }, error = function(e) {
    message("Error in merging annotations: ", e$message)
    if (reticulate::py_has_attr(reticulate::py, "last_error")) {
      message("Python traceback: ", reticulate::py_last_error())
    }
    stop("Failed to run annotation merging. Make sure Python dependencies are installed.")
  })
}

