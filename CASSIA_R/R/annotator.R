#' @import reticulate

# Use reticulate to import Python modules
py_main <- NULL
py_tools <- NULL
py_determinator <- NULL  # Add this line

.onLoad <- function(libname, pkgname) {
  # Get the conda environment name from the package configuration
  conda_env <- getOption("CASSIA.conda_env", default = "cassia_env")
  
  # Set up the Python environment
  tryCatch({
    # Check if the environment exists, if not create it
    if (!conda_env %in% reticulate::conda_list()$name) {
      setup_cassia_env(conda_env)
    }
    
    reticulate::use_condaenv(conda_env, required = TRUE)
    
    # Import Python modules
    py_main <<- reticulate::import_from_path("main_function_code", path = system.file("python", package = "CASSIA"))
    py_tools <<- reticulate::import_from_path("tools_function", path = system.file("python", package = "CASSIA"))
    py_merging <<- reticulate::import_from_path("merging_annotation_code", path = system.file("python", package = "CASSIA"))
  }, error = function(e) {
    warning("Failed to set up Python environment. Please run setup_cassia_env() manually to set up the required environment.")
  })
}



#' Set up CASSIA Python Environment
#'
#' This function sets up the required Python environment for CASSIA.
#' It can be used to create a new environment or update an existing one.
#' 
#' @param conda_env The name of the conda environment to use. If NULL, uses the default from package configuration.
#' @param python_version The Python version to use. Default is "3.10".
#' @param pip_packages A character vector of pip packages to install.
#'
#' @return Invisible NULL. Called for side effects.
#' @export
setup_cassia_env <- function(conda_env = NULL, python_version = "3.10", 
                           pip_packages = c("openai", "pandas", "numpy", "scikit-learn", 
                                          "requests", "anthropic", "charset-normalizer")) {
  if (is.null(conda_env)) {
    conda_env <- getOption("CASSIA.conda_env", default = "cassia_env")
  }
  
  # Check if the environment exists
  if (!conda_env %in% reticulate::conda_list()$name) {
    # Create conda environment
    reticulate::conda_create(envname = conda_env, python_version = python_version)
  }
  
  # Use the environment
  reticulate::use_condaenv(conda_env, required = TRUE)
  
  # Install or update required packages
  reticulate::py_install(pip_packages, pip = TRUE)
  
  # Set the CASSIA.conda_env option
  options(CASSIA.conda_env = conda_env)
  
  invisible(NULL)
}


#' Set API Key for LLM Provider
#'
#' @param api_key Character string containing the API key
#' @param provider Character string specifying the provider ('openai', 'anthropic', or 'openrouter')
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
setLLMApiKey <- function(api_key, provider = "anthropic", persist = FALSE) {
  # Set environment variable based on provider
  if (provider == "openai") {
    Sys.setenv(OPENAI_API_KEY = api_key)
    env_var_name <- "OPENAI_API_KEY"
  } else if (provider == "anthropic") {
    Sys.setenv(ANTHROPIC_API_KEY = api_key)
    env_var_name <- "ANTHROPIC_API_KEY"
  } else if (provider == "openrouter") {
    Sys.setenv(OPENROUTER_API_KEY = api_key)
    env_var_name <- "OPENROUTER_API_KEY"
  } else {
    stop("Unsupported provider. Use 'openai' or 'anthropic' or 'openrouter'")
  }
  
  # Set API key in Python tools if available
  if (!is.null(py_tools)) {
    py_tools$set_api_key(api_key, provider)
  }
  
  # Persist to .Renviron if requested
  if (persist) {
    renviron_path <- path.expand("~/.Renviron")
    
    # Read existing content
    if (file.exists(renviron_path)) {
      lines <- readLines(renviron_path)
      # Remove existing entry if present
      lines <- lines[!grepl(paste0("^", env_var_name, "="), lines)]
    } else {
      lines <- character()
    }
    
    # Add new entry
    lines <- c(lines, paste0(env_var_name, "='", api_key, "'"))
    
    # Write back to .Renviron
    writeLines(lines, renviron_path)
    
    # Reload .Renviron
    readRenviron("~/.Renviron")
    
    message("API key has been saved to .Renviron and will be loaded in future sessions")
  }
  
  invisible(NULL)
}


#' Set OpenAI API Key
#'
#' @param api_key Character string containing the OpenAI API key
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
set_openai_api_key <- function(api_key, persist = FALSE) {
  setLLMApiKey(api_key, provider = "openai", persist = persist)
}


#' Set Anthropic API Key
#'
#' @param api_key Character string containing the Anthropic API key
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
setAnthropicApiKey <- function(api_key, persist = FALSE) {
  setLLMApiKey(api_key, provider = "anthropic", persist = persist)
}


#' Set OpenRouter API Key
#'
#' @param api_key Character string containing the OpenRouter API key
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
setOpenRouterApiKey <- function(api_key, persist = FALSE) {
  setLLMApiKey(api_key, provider = "openrouter", persist = persist)
}



#' Run Cell Type Analysis
#'
#' @param model Character string specifying the model to use.
#' @param temperature Numeric value for temperature parameter.
#' @param marker_list List of marker genes.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter', default='openai')
#'
#' @return A list containing two elements: structured output and conversation history.
#' @export
runCASSIA <- function(model = "google/gemini-2.5-flash-preview", temperature, marker_list, tissue, species, additional_info = NULL, provider = "openrouter") {
  tryCatch({
    result <- py_tools$run_cell_type_analysis_wrapper(
      model = model,
      temperature = temperature,
      marker_list = marker_list,
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      provider = provider
    )
    
    # Convert structured_output (result[[1]])
    structured_output <- as.list(result[[1]])
    
    # Convert conversation_history (result[[2]])
    conversation_history <- lapply(result[[2]], function(entry) {
      if (is.list(entry) && length(entry) == 2) {
        list(agent = as.character(entry[[1]]), message = as.character(entry[[2]]))
      } else {
        warning("Unexpected entry structure in conversation_history")
        list(agent = "Unknown", message = "Conversion failed")
      }
    })
    
    return(list(structured_output = structured_output, conversation_history = conversation_history))
  }, error = function(e) {
    error_msg <- paste("Error in run_cell_type_analysis:", e$message, "\n",
                       "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}

#' Run Cell Type Analysis Multiple Times
#'
#' @param n Number of times to run the analysis.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param temperature Numeric value for temperature parameter.
#' @param marker_list List of marker genes.
#' @param model Character string specifying the model to use.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#'
#' @return A list containing results from multiple runs.
#' @export
runCASSIA_n_times <- function(n, tissue, species, additional_info, temperature, marker_list, 
                           model = "google/gemini-2.5-flash-preview", max_workers = 10, provider = "openrouter") {
  tryCatch({
    result <- py_tools$run_analysis_n_times(
      n = as.integer(n),
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      temperature = as.numeric(temperature),
      marker_list = marker_list,
      model = model,
      max_workers = as.integer(max_workers),
      provider = provider
    )
    
    # Convert the result to an R list
    converted_result <- lapply(seq_len(n), function(i) {
      run_result <- result[[as.character(i-1)]]  # Python uses 0-based indexing
      list(
        analysis_result = list(
          main_cell_type = run_result[[1]]$main_cell_type,
          sub_cell_types = run_result[[1]]$sub_cell_types,
          possible_mixed_cell_types = run_result[[1]]$possible_mixed_cell_types,
          num_markers = run_result[[1]]$num_markers,
          iterations = run_result[[1]]$iterations
        ),
        conversation_history = lapply(run_result[[2]], function(entry) {
          list(agent = entry[[1]], message = entry[[2]])
        })
      )
    })
    
    names(converted_result) <- as.character(seq_len(n) - 1)  # Match Python's 0-based indexing
    return(converted_result)
  }, error = function(e) {
    error_msg <- paste("Error in run_analysis_n_times:", e$message, "\n",
                       "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}




#' Process Cell Type Analysis Results
#'
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param temperature Numeric value for temperature parameter.
#' @param marker_list List of marker genes.
#' @param model Character string specifying the model to use.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param n Number of times to run the analysis.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#'
#' @return A list containing processed results including variance analysis.
#' @export
runCASSIA_n_times_similarity_score <- function(tissue, species, additional_info, temperature, marker_list, model = "google/gemini-2.5-flash-preview", max_workers, n, provider = "openrouter") {
  tryCatch({
    # Call the Python function with the new parameter structure
    processed_results <- py_tools$process_cell_type_analysis_single_wrapper(
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      temperature = as.numeric(temperature),
      marker_list = marker_list,
      model = model,
      max_workers = as.integer(max_workers),
      n = as.integer(n),
      provider=provider
    )
    
    # Convert the processed results back to R
    r_results <- list(
      unified_results = as.character(processed_results$unified_results),
      consensus_types = list(
        general = as.character(processed_results$consensus_types[[1]]),
        sub = as.character(processed_results$consensus_types[[2]])
      ),
      general_celltype_llm = as.character(processed_results$general_celltype_llm),
      sub_celltype_llm = as.character(processed_results$sub_celltype_llm),
      Possible_mixed_celltypes_llm = as.list(processed_results$Possible_mixed_celltypes_llm),
      similarity_score = as.numeric(processed_results$similarity_score),
      original_results = lapply(processed_results$original_results, function(x) as.character(x)),
      llm_response = as.character(processed_results$llm_response)
    )
    
    return(r_results)
  }, error = function(e) {
    error_msg <- paste("Error in process_cell_type_analysis_single:", e$message, "\n",
                       "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}







#' Run Cell Type Analysis Batch
#'
#' @param df_input A data frame or path to the CSV file containing marker data.
#' @param output_name Name of the output JSON file.
#' @param model Character string specifying the model to use.
#' @param temperature Numeric value for temperature parameter.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param celltype_column Name of the column containing cell types.
#' @param gene_column_name Name of the column containing gene names.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param n_genes Number of top genes to use (default: 50)
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#'
#' @return None. This function creates output files and prints execution time.
#' @export
runCASSIA_batch <- function(marker, output_name = "cell_type_analysis_results.json", 
                          model = "google/gemini-2.5-flash-preview", temperature = 0, tissue = "lung", 
                          species = "human", additional_info = NULL, 
                          celltype_column = NULL, gene_column_name = NULL, 
                          max_workers = 10, provider = "openrouter", n_genes = 50,
                          max_retries = 1) {
  execution_time <- system.time({
    # Convert R dataframe to Python if df_input is a dataframe
if (is.data.frame(marker)) {
  pd <- reticulate::import("pandas")
  # Direct conversion with convert=TRUE
  marker <- reticulate::r_to_py(marker, convert = TRUE)
} else if (!is.character(marker)) {
  stop("marker must be either a data frame or a character vector")
}
    
    py_tools$run_cell_type_analysis_batchrun(
      marker = marker,  # Changed parameter name to match Python function
      output_name = output_name,
      model = model,
      temperature = temperature,
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      celltype_column = celltype_column,
      gene_column_name = gene_column_name,
      max_workers = as.integer(max_workers),
      provider = provider,
      n_genes = as.integer(n_genes),
      max_retries = as.integer(max_retries)
    )
  })
  
  print(paste("Execution time for runCASSIA_batch:"))
  print(execution_time)
}





#' Run Batch Analysis Multiple Times
#'
#' @param n Number of times to run the batch analysis.
#' @param marker Path to the CSV file containing marker data.
#' @param output_name Prefix for output JSON files.
#' @param model Character string specifying the model to use.
#' @param temperature Numeric value for temperature parameter.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param celltype_column Name of the column containing cell types.
#' @param gene_column_name Name of the column containing gene names.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param batch_max_workers Maximum number of workers for batch processing.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#'
#' @return None. This function creates output files and prints execution time.
#' @export
runCASSIA_batch_n_times <- function(n, marker, output_name = "cell_type_analysis_results", 
                                  model = "google/gemini-2.5-flash-preview", temperature = 0, tissue = "lung", 
                                  species = "human", additional_info = NULL, 
                                  celltype_column = NULL, gene_column_name = NULL, 
                                  max_workers = 10, batch_max_workers = 5, 
                                  provider = "openrouter", max_retries = 1) {

  if (is.data.frame(marker)) {
  pd <- reticulate::import("pandas")
  # Direct conversion with convert=TRUE
  marker <- reticulate::r_to_py(marker, convert = TRUE)
} else if (!is.character(marker)) {
  stop("marker must be either a data frame or a character vector")
}

  execution_time <- system.time({
    tryCatch({
      py_tools$run_batch_analysis_n_times(
        as.integer(n), marker, output_name, model, temperature, tissue, 
        species, additional_info, celltype_column, gene_column_name, 
        as.integer(max_workers), as.integer(batch_max_workers), provider,
        as.integer(max_retries)
      )
    }, error = function(e) {
      stop(paste("Error in run_batch_analysis_n_times:", e$message))
    })
  })
  
  print(paste("Execution time for runCASSIA_batch_n_times:"))
  print(execution_time)
}




#' Process and Save Batch Results
#'
#' @param marker Path to the marker file.
#' @param file_pattern Pattern to match result files.
#' @param output_name Name of the output CSV file.
#' @param celltype_column Name of the column containing cell types.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param model Model to use for processing (default: "gpt-4o")
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param main_weight Weight for the main cell type.
#' @param sub_weight Weight for the sub cell type.
#'
#' @return None. This function processes and saves results to a CSV file and prints execution time.
#' @export
runCASSIA_similarity_score_batch <- function(marker, file_pattern, output_name, 
                                               celltype_column = NULL, max_workers = 10, model = "google/gemini-2.5-flash-preview", provider = "openrouter", main_weight=0.5, sub_weight=0.5) {


  if (is.data.frame(marker)) {
  pd <- reticulate::import("pandas")
  # Direct conversion with convert=TRUE
  marker <- reticulate::r_to_py(marker, convert = TRUE)
} else if (!is.character(marker)) {
  stop("marker must be either a data frame or a character vector")
}

  execution_time <- system.time({
    tryCatch({
      py_tools$process_and_save_batch_results(marker, file_pattern, output_name, 
                                              celltype_column, max_workers, model, provider, main_weight, sub_weight)
    }, error = function(e) {
      stop(paste("Error in process_and_save_batch_results:", e$message))
    })
  })
  
  print(paste("Execution time for runCASSIA_batch_get_similarity_score:"))
  print(execution_time)
}







#' Generate Cell Type Analysis Report
#'
#' @param full_result_path Path to the full results CSV file
#' @param marker_path Path to the marker genes CSV file
#' @param cluster_name Name of the cluster to analyze
#' @param major_cluster_info General information about the dataset (e.g., "Human PBMC")
#' @param output_name Name of the output HTML file
#' @param num_iterations Number of iterations for marker analysis (default=5)
#' @param model Model to use for analysis (default="gpt-4o")
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#'
#' @return None
#' @export
runCASSIA_annotationboost<- function(full_result_path, 
                                                     marker, 
                                                     cluster_name, 
                                                     major_cluster_info, 
                                                     output_name, 
                                                     num_iterations = 5,
                                                     model = "google/gemini-2.5-flash-preview",
                                                     provider = "openrouter") {

  if (is.data.frame(marker)) {
  pd <- reticulate::import("pandas")
  # Direct conversion with convert=TRUE
  marker <- reticulate::r_to_py(marker, convert = TRUE)
} else if (!is.character(marker)) {
  stop("marker must be either a data frame or a character vector")
}

                                                      
  tryCatch({
    result <- py_tools$generate_cell_type_analysis_report_wrapper(
      full_result_path = full_result_path,
      marker = marker,
      cluster_name = cluster_name,
      major_cluster_info = major_cluster_info,
      output_name = output_name,
      num_iterations = as.integer(num_iterations),
      model = model,
      provider = provider
    )
    
    # Convert the result to an R list
    analysis_result <- result[[1]]
    messages_history <- lapply(result[[2]], function(msg) {
      list(
        role = msg$role,
        content = msg$content
      )
    })
    
    
  }, error = function(e) {
    error_msg <- paste("Error in generate_cell_type_analysis_report_wrapper:", e$message)
    stop(error_msg)
  })
}




#' Generate Cell Type Analysis Report
#'
#' @param full_result_path Path to the full results CSV file
#' @param marker_path Path to the marker genes CSV file
#' @param cluster_name Name of the cluster to analyze
#' @param major_cluster_info General information about the dataset (e.g., "Human PBMC")
#' @param output_name Name of the output HTML file
#' @param num_iterations Number of iterations for marker analysis (default=5)
#' @param model Model to use for analysis (default="gpt-4o")
#' @param additional_task Additional task to perform
#' @return None
#' @export
runCASSIA_annottaionboost_additional_task<- function(full_result_path, 
                                                     marker, 
                                                     cluster_name, 
                                                     major_cluster_info, 
                                                     output_name, 
                                                     num_iterations = 5,
                                                     model = "google/gemini-2.5-flash-preview",
                                                     additional_task = "") {

  if (is.data.frame(marker)) {
  pd <- reticulate::import("pandas")
  # Direct conversion with convert=TRUE
  marker <- reticulate::r_to_py(marker, convert = TRUE)
} else if (!is.character(marker)) {
  stop("marker must be either a data frame or a character vector")
}

                                                      
  tryCatch({
    result <- py_tools$generate_cell_type_analysis_report_openrouter_additional_task(
      full_result_path = full_result_path,
      marker = marker,
      cluster_name = cluster_name,
      major_cluster_info = major_cluster_info,
      output_name = output_name,
      num_iterations = as.integer(num_iterations),
      model = model,
      additional_task = additional_task
    )
    
    # Convert the result to an R list
    analysis_result <- result[[1]]
    messages_history <- lapply(result[[2]], function(msg) {
      list(
        role = msg$role,
        content = msg$content
      )
    })
    
    
  }, error = function(e) {
    error_msg <- paste("Error in generate_cell_type_analysis_report_wrapper:", e$message)
    stop(error_msg)
  })
}


#' Run Scoring with Progress Updates
#'
#' @param input_file Path to input CSV file
#' @param output_file Path to output CSV file (optional)
#' @param max_workers Maximum number of parallel workers
#' @param model Model to use
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#'
#' @return None
#' @export
runCASSIA_score_batch <- function(input_file, 
                                    output_file = NULL, 
                                    max_workers = 4, 
                                    model = "deepseek/deepseek-chat-v3-0324",
                                    provider = "openrouter",
                                    max_retries = 1) {
  tryCatch({
    results <- py_tools$run_scoring_with_progress(
      input_file = input_file,
      output_file = output_file,
      max_workers = as.integer(max_workers),
      model = model,
      provider = provider,
      max_retries = as.integer(max_retries)
    )
    
    # Convert Python DataFrame to R data.frame if results is not NULL
    if (!is.null(results)) {
      results <- reticulate::py_to_r(results)
    }
    
    
  }, error = function(e) {
    error_msg <- paste("Error in run_scoring_with_progress:", e$message)
    stop(error_msg)
  })
}

#' Generate HTML Reports from Scored Results
#'
#' This function processes a CSV file containing scored annotation results and generates 
#' individual HTML reports for each row, along with an index page.
#'
#' @param csv_path Character string. Path to the CSV file containing scored results.
#' @param index_name Character string. Base name for the index file (default: "CASSIA_reports_summary").
#'
#' @details 
#' The function generates:
#' 1. Individual HTML reports for each annotation result
#' 2. An index page linking to all generated reports
#' 
#' The CSV file should contain columns:
#' - Conversation History
#' - Scoring_Reasoning
#' - Score
#'
#' @return None. Files are written to the current directory.
#' @export
#'
#' @examples
#' \dontrun{
#' runCASSIA_generate_score_report("path/to/scored_results.csv")
#' }
runCASSIA_generate_score_report <- function(csv_path, output_name = "CASSIA_reports_summary") {
  tryCatch({
    py_tools$process_all_reports(
      csv_path = csv_path,
      index_name = output_name
    )
    message("Reports generated successfully. Check the current directory for the generated HTML files.")
  }, error = function(e) {
    error_msg <- paste("Error in generating reports:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}



#' Run Complete Cell Analysis Pipeline
#'
#' @param output_file_name Base name for output files
#' @param tissue Tissue type being analyzed
#' @param species Species being analyzed
#' @param marker Marker data (data frame or file path)
#' @param max_workers Maximum number of concurrent workers (default: 4)
#' @param annotation_model Model to use for initial annotation (default: "meta-llama/llama-4-maverick")
#' @param annotation_provider Provider for initial annotation (default: "openrouter")
#' @param score_model Model to use for scoring (default: "google/gemini-2.5-pro-preview-03-25")
#' @param score_provider Provider for scoring (default: "openrouter")
#' @param annotationboost_model Model to use for boosting low-scoring annotations (default: "google/gemini-2.5-flash-preview")
#' @param annotationboost_provider Provider for boosting low-scoring annotations (default: "openrouter")
#' @param score_threshold Threshold for identifying low-scoring clusters (default: 75)
#' @param additional_info Additional information for analysis (default: NULL)
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#' @param do_merge_annotations Whether to run the merging annotations step (default: TRUE)
#' @param merge_model Model to use for merging annotations (default: "deepseek/deepseek-chat-v3-0324")
#'
#' @return None. Creates output files and generates reports.
#' @export
runCASSIA_pipeline <- function(
    output_file_name,
    tissue,
    species,
    marker,
    max_workers = 4,
    annotation_model = "google/gemini-2.5-flash-preview",
    annotation_provider = "openrouter",
    score_model = "deepseek/deepseek-chat-v3-0324",
    score_provider = "openrouter",
    annotationboost_model = "google/gemini-2.5-flash-preview",
    annotationboost_provider = "openrouter",
    score_threshold = 75,
    additional_info = NULL,
    max_retries = 1,
    do_merge_annotations = TRUE,
    merge_model = "deepseek/deepseek-chat-v3-0324"
) {
  # Convert marker data frame if necessary
  if (is.data.frame(marker)) {
    pd <- reticulate::import("pandas")
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }

  tryCatch({
    py_tools$run_cell_analysis_pipeline(
      output_file_name = output_file_name,
      tissue = tissue,
      species = species,
      marker_path = marker,
      max_workers = as.integer(max_workers),
      annotation_model = annotation_model,
      annotation_provider = annotation_provider,
      score_model = score_model,
      score_provider = score_provider,
      annotationboost_model = annotationboost_model,
      annotationboost_provider = annotationboost_provider,
      score_threshold = as.numeric(score_threshold),
      additional_info = if(is.null(additional_info)) "None" else additional_info,
      max_retries = as.integer(max_retries),
      do_merge_annotations = do_merge_annotations,
      merge_model = merge_model
    )
  }, error = function(e) {
    error_msg <- paste("Error in run_cell_analysis_pipeline:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}



#' Compare Cell Types Using Different Models
#'
#' @param tissue Character string specifying the tissue type
#' @param celltypes Character vector of 2-4 cell types to compare
#' @param marker Character vector or string of marker genes
#' @param species Character string specifying the species (default: "human")
#' @param model_list Character vector of model names to use (optional)
#'
#' @return A list containing responses from different models
#' @export
compareCelltypes <- function(tissue, celltypes, marker, species = "human", model_list = NULL,output_file = NULL) {
  tryCatch({
    # Input validation
    if (length(celltypes) < 2 || length(celltypes) > 5) {
      stop("Please provide 2-5 cell types to compare")
    }
    
    # Convert marker_set to string if it's a vector
    if (is.vector(marker)) {
      marker <- paste(marker, collapse = ", ")
    }
    
    # Call the Python function
    responses <- py_tools$compare_celltypes(
      tissue = tissue,
      celltypes = celltypes,
      marker_set = marker,
      species = species,
      model_list = model_list,
      output_file = output_file
    )

  }, error = function(e) {
    error_msg <- paste("Error in compareModels:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}






#' Process Subclusters
#'
#' @param marker Marker data (data frame or file path)
#' @param major_cluster_info Description of the major cluster type
#' @param output_name Base name for the output file (will add .csv if not present)
#' @param model Model name for Claude API (default: "claude-3-5-sonnet-20241022")
#' @param temperature Temperature parameter for API calls (default: 0)
#' @param provider AI provider to use (default: "anthropic")
#' @param n_genes Number of top genes to use (default: 50)
#'
#' @return None. This function processes subclusters and saves results to a CSV file.
#' @export
runCASSIA_subclusters <- function(marker, major_cluster_info, output_name, 
                               model = "google/gemini-2.5-flash-preview", temperature = 0, 
                               provider = "openrouter", n_genes = 50L) {
  py_tools$runCASSIA_subclusters(
    marker = marker,
    major_cluster_info = major_cluster_info,
    output_name = output_name,
    model = model,
    temperature = as.numeric(temperature),
    provider = provider,
    n_genes = as.integer(n_genes)
  )
}

#' Run Analysis Multiple Times for Subclusters
#'
#' @param n Number of times to run the analysis
#' @param marker Marker data (data frame or file path)
#' @param major_cluster_info Description of the major cluster type
#' @param base_output_name Base name for output CSV files
#' @param model Model name for Claude API (default: "claude-3-5-sonnet-20241022")
#' @param temperature Temperature parameter for API calls (default: 0)
#' @param provider AI provider to use (default: "anthropic")
#' @param max_workers Maximum number of workers for parallel processing (default: 5)
#' @param n_genes Number of top genes to use (default: 50)
#'
#' @return None. This function runs the analysis multiple times and saves results to CSV files.
#' @export
runCASSIA_n_subcluster <- function(n, marker, major_cluster_info, base_output_name, 
                                               model = "google/gemini-2.5-flash-preview", temperature = 0, 
                                               provider = "openrouter", max_workers = 5,n_genes=50L) {
  py_tools$runCASSIA_n_subcluster(
    n = as.integer(n),
    marker = marker,
    major_cluster_info = major_cluster_info,
    base_output_name = base_output_name,
    model = model,
    temperature = temperature,
    provider = provider,
    max_workers = as.integer(max_workers),
    n_genes = as.integer(n_genes)
  )
}