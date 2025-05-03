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
#' @importFrom reticulate import source_python py_set_seed use_python
#' @export
runCASSIA_merge_annotations <- function(csv_path,
                             output_path = NULL,
                             provider = "openrouter",
                             model = "deepseek/deepseek-chat-v3-0324",
                             additional_context = NULL,
                             batch_size = 20,
                             detail_level = "broad",
                             process_all = FALSE,
                             debug = FALSE) {

  tryCatch({
    # Verify that the input file exists
    if (!file.exists(csv_path)) {
      stop("Input CSV file does not exist: ", csv_path)
    }
    
    # Check if the input file is already open (basic check on Windows)
    if (is.null(output_path)) {
      tryCatch({
        test_con <- file(csv_path, "r+")
        close(test_con)
      }, error = function(e) {
        stop("The input CSV file appears to be open in another program. ",
             "Please close the file or specify a different output_path.")
      })
    }
    
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
    
    # Get the Python module path
    py_module_path <- system.file("python", package = "CASSIA")
    if (py_module_path == "") {
      # If system.file returns empty string, the directory doesn't exist in the package
      # Try to find it relative to the current script
      package_root <- find.package("CASSIA")
      py_module_path <- file.path(package_root, "inst", "python")
      if (!dir.exists(py_module_path)) {
        stop("Cannot find the Python module directory. Looked in: ", 
             system.file("python", package = "CASSIA"), " and ", py_module_path)
      }
    }
    
    if (debug) {
      message("Python module path: ", py_module_path)
      message("Files in the module directory: ", paste(list.files(py_module_path), collapse = ", "))
    }
    
    # Source the Python file directly instead of importing it
    py_file_path <- file.path(py_module_path, "merging_annotation.py")
    if (!file.exists(py_file_path)) {
      stop("Python module file not found: ", py_file_path)
    }
    
    # Import the module or source the file
    tryCatch({
      # Try direct import first
      message("Importing Python module from: ", py_module_path)
      reticulate::source_python(py_file_path)
      # Create references to the functions
      merge_annotations <- reticulate::py_get_attr(reticulate::py, "merge_annotations")
      merge_annotations_all <- reticulate::py_get_attr(reticulate::py, "merge_annotations_all")
      
      message("Successfully loaded Python functions")
    }, error = function(e) {
      stop("Failed to import Python module: ", e$message)
    })
    
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
      message("Results will be saved to: ", output_path)
    }
    # Otherwise, Python will save back to input file
    else {
      message("No output_path specified, results will be saved back to the input file: ", csv_path)
    }
    
    # Print debug info if requested
    if (debug) {
      message("Function arguments: ")
      for (name in names(args)) {
        message("  ", name, ": ", args[[name]])
      }
      if (process_all) {
        message("Will call merge_annotations_all")
      } else {
        message("Will call merge_annotations with detail_level: ", detail_level)
        args$detail_level <- detail_level
      }
    }
    
    # Choose which Python function to call based on process_all parameter
    if (process_all) {
      # Process all detail levels sequentially by calling merge_annotations three times
      message("Processing all detail levels sequentially...")
      
      # Define all three detail levels
      all_detail_levels <- c("broad", "detailed", "very_detailed")
      
      # Call merge_annotations for each detail level
      for (level in all_detail_levels) {
        message(paste0("Processing ", level, " level groupings..."))
        merge_annotations(
          csv_path = args$csv_path,
          output_path = args$output_path,
          provider = args$provider,
          model = args$model,
          additional_context = args$additional_context,
          batch_size = args$batch_size,
          detail_level = level
        )
        message(paste0("Completed ", level, " level groupings"))
      }
      
      message("All detail levels processed successfully")
    } else {
      # Call the function that processes a single detail level
      args$detail_level <- detail_level
      message(paste0("Processing ", detail_level, " level groupings..."))
      merge_annotations(
        csv_path = args$csv_path,
        output_path = args$output_path,
        provider = args$provider,
        model = args$model,
        additional_context = args$additional_context,
        batch_size = args$batch_size,
        detail_level = args$detail_level
      )
      message("Processing completed successfully")
    }
    
    # Return TRUE invisibly on success
    invisible(TRUE)
    
  }, error = function(e) {
    # Check for permission error
    if (grepl("Permission denied", e$message)) {
      target_file <- if (is.null(output_path)) csv_path else output_path
      stop("Permission denied when saving to ", target_file, 
           ". Please close the file if it's open in another program, or specify a different output_path.")
    }
    
    # Other errors
    message("Error in merging annotations: ", e$message)
    if (reticulate::py_has_attr(reticulate::py, "last_error")) {
      message("Python traceback: ", reticulate::py_last_error())
    }
    stop("Failed to run annotation merging. Make sure Python dependencies are installed.")
  })
}

