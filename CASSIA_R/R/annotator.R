#' @import reticulate

# Use reticulate to import Python modules
py_main <- NULL
py_tools <- NULL
py_merging <- NULL
py_annotation_boost <- NULL
py_cell_comparison <- NULL
py_subclustering <- NULL
py_uncertainty <- NULL
py_generate_reports <- NULL
py_llm_utils <- NULL
py_generate_hypothesis_report <- NULL
py_hypothesis_geneartion <- NULL
py_summarize_hypothesis_runs <- NULL
py_debug_genes <- NULL
py_super_annottaion_boost <- NULL
py_symphony_compare <- NULL
py_model_settings <- NULL

.onLoad <- function(libname, pkgname) {
  # Get the environment name from the package configuration
  # Support both new and legacy option names
  env_name <- getOption("CASSIA.env_name", default = NULL)
  if (is.null(env_name)) {
    env_name <- getOption("CASSIA.conda_env", default = "cassia_env")
  }
  
  # Set up the Python environment
  tryCatch({
    # Check if environment exists in either virtualenv or conda
    env_exists <- .check_env_exists(env_name)
    
    if (!env_exists$virtualenv && !env_exists$conda) {
      # Environment doesn't exist, create it using the new setup function
      message("CASSIA Python environment not found. Setting up automatically...")
      setup_cassia_env(conda_env = env_name)
    } else {
      # Environment exists, activate it using the appropriate method
      if (env_exists$virtualenv) {
        reticulate::use_virtualenv(env_name, required = TRUE)
        options(CASSIA.env_name = env_name)
        options(CASSIA.env_method = "virtualenv")
      } else if (env_exists$conda) {
        reticulate::use_condaenv(env_name, required = TRUE)
        options(CASSIA.env_name = env_name)
        options(CASSIA.env_method = "conda")
        # Maintain backward compatibility
        options(CASSIA.conda_env = env_name)
      }
    }
    
    # Import Python modules
    py_main <<- reticulate::import_from_path("main_function_code", path = system.file("python", package = "CASSIA"))
    py_tools <<- reticulate::import_from_path("tools_function", path = system.file("python", package = "CASSIA"))
    py_merging <<- reticulate::import_from_path("merging_annotation", path = system.file("python", package = "CASSIA"))
    py_annotation_boost <<- reticulate::import_from_path("annotation_boost", path = system.file("python", package = "CASSIA"))
    py_cell_comparison <<- reticulate::import_from_path("cell_type_comparison", path = system.file("python", package = "CASSIA"))
    py_subclustering <<- reticulate::import_from_path("subclustering", path = system.file("python", package = "CASSIA"))
    py_uncertainty <<- reticulate::import_from_path("Uncertainty_quantification", path = system.file("python", package = "CASSIA"))
    py_generate_reports <<- reticulate::import_from_path("generate_reports", path = system.file("python", package = "CASSIA"))
    py_llm_utils <<- reticulate::import_from_path("llm_utils", path = system.file("python", package = "CASSIA"))
    py_generate_hypothesis_report <<- reticulate::import_from_path("generate_hypothesis_report", path = system.file("python", package = "CASSIA"))
    py_hypothesis_geneartion <<- reticulate::import_from_path("hypothesis_geneartion", path = system.file("python", package = "CASSIA"))
    py_summarize_hypothesis_runs <<- reticulate::import_from_path("summarize_hypothesis_runs", path = system.file("python", package = "CASSIA"))
    py_debug_genes <<- reticulate::import_from_path("debug_genes", path = system.file("python", package = "CASSIA"))
    py_super_annottaion_boost <<- reticulate::import_from_path("super_annottaion_boost", path = system.file("python", package = "CASSIA"))
    py_symphony_compare <<- reticulate::import_from_path("symphony_compare", path = system.file("python", package = "CASSIA"))
    py_model_settings <<- reticulate::import_from_path("model_settings", path = system.file("python", package = "CASSIA"))
    
    message("CASSIA loaded successfully!")
    
  }, error = function(e) {
    # If setup fails, try to run setup_cassia_env() automatically
    message("Initial setup failed, attempting automatic environment setup...")
    tryCatch({
      setup_cassia_env(conda_env = env_name)
      
      # Try to import Python modules again after successful setup
      py_main <<- reticulate::import_from_path("main_function_code", path = system.file("python", package = "CASSIA"))
      py_tools <<- reticulate::import_from_path("tools_function", path = system.file("python", package = "CASSIA"))
      py_merging <<- reticulate::import_from_path("merging_annotation", path = system.file("python", package = "CASSIA"))
      py_annotation_boost <<- reticulate::import_from_path("annotation_boost", path = system.file("python", package = "CASSIA"))
      py_cell_comparison <<- reticulate::import_from_path("cell_type_comparison", path = system.file("python", package = "CASSIA"))
      py_subclustering <<- reticulate::import_from_path("subclustering", path = system.file("python", package = "CASSIA"))
      py_uncertainty <<- reticulate::import_from_path("Uncertainty_quantification", path = system.file("python", package = "CASSIA"))
      py_generate_reports <<- reticulate::import_from_path("generate_reports", path = system.file("python", package = "CASSIA"))
      py_llm_utils <<- reticulate::import_from_path("llm_utils", path = system.file("python", package = "CASSIA"))
      py_generate_hypothesis_report <<- reticulate::import_from_path("generate_hypothesis_report", path = system.file("python", package = "CASSIA"))
      py_hypothesis_geneartion <<- reticulate::import_from_path("hypothesis_geneartion", path = system.file("python", package = "CASSIA"))
      py_summarize_hypothesis_runs <<- reticulate::import_from_path("summarize_hypothesis_runs", path = system.file("python", package = "CASSIA"))
      py_debug_genes <<- reticulate::import_from_path("debug_genes", path = system.file("python", package = "CASSIA"))
      py_super_annottaion_boost <<- reticulate::import_from_path("super_annottaion_boost", path = system.file("python", package = "CASSIA"))
      py_symphony_compare <<- reticulate::import_from_path("symphony_compare", path = system.file("python", package = "CASSIA"))
      py_model_settings <<- reticulate::import_from_path("model_settings", path = system.file("python", package = "CASSIA"))
      
      message("CASSIA environment setup completed successfully!")
      
    }, error = function(e2) {
      warning("Failed to set up Python environment automatically. Error: ", e2$message, 
              "\nPlease run setup_cassia_env() manually to set up the required environment.")
    })
  })
}



# Helper function to try virtualenv setup
.try_virtualenv_setup <- function(env_name, python_version, pip_packages) {
  tryCatch({
    # Check if virtualenv already exists
    existing_envs <- tryCatch(reticulate::virtualenv_list(), error = function(e) character(0))
    
    if (!env_name %in% existing_envs) {
      # Create virtualenv
      reticulate::virtualenv_create(envname = env_name, python_version = python_version)
    }
    
    # Use the virtualenv
    reticulate::use_virtualenv(env_name, required = TRUE)
    
    # Install packages
    reticulate::virtualenv_install(envname = env_name, packages = pip_packages)
    
    # Set options to track the environment and method
    options(CASSIA.env_name = env_name)
    options(CASSIA.env_method = "virtualenv")
    
    message("Successfully set up CASSIA environment using virtualenv: ", env_name)
    return(TRUE)
    
  }, error = function(e) {
    warning("Virtualenv setup failed: ", e$message)
    return(FALSE)
  })
}

# Helper function to try conda setup (current logic)
.try_conda_setup <- function(conda_env, python_version, pip_packages) {
  tryCatch({
    # Check if conda environment exists
    existing_envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))
    
    if (!conda_env %in% existing_envs) {
      # Create conda environment
      reticulate::conda_create(envname = conda_env, python_version = python_version)
    }
    
    # Use the conda environment
    reticulate::use_condaenv(conda_env, required = TRUE)
    
    # Install packages
    reticulate::py_install(pip_packages, pip = TRUE)
    
    # Set options to track the environment and method
    options(CASSIA.env_name = conda_env)
    options(CASSIA.env_method = "conda")
    
    message("Successfully set up CASSIA environment using conda: ", conda_env)
    return(TRUE)
    
  }, error = function(e) {
    warning("Conda setup failed: ", e$message)
    return(FALSE)
  })
}

# Helper function to get environment name
.get_env_name <- function(env_name) {
  if (is.null(env_name)) {
    # Check for existing environment name from options
    existing_env <- getOption("CASSIA.env_name", default = NULL)
    if (!is.null(existing_env)) {
      return(existing_env)
    }
    # Fall back to legacy conda env option
    return(getOption("CASSIA.conda_env", default = "cassia_env"))
  }
  return(env_name)
}

# Helper function to check if environment exists
.check_env_exists <- function(env_name) {
  # Check virtualenv
  virtualenv_exists <- tryCatch({
    existing_virtualenvs <- reticulate::virtualenv_list()
    env_name %in% existing_virtualenvs
  }, error = function(e) FALSE)
  
  # Check conda
  conda_exists <- tryCatch({
    existing_condas <- reticulate::conda_list()$name
    env_name %in% existing_condas
  }, error = function(e) FALSE)
  
  return(list(virtualenv = virtualenv_exists, conda = conda_exists))
}

#' Set up CASSIA Python Environment
#'
#' This function sets up the required Python environment for CASSIA.
#' It can be used to create a new environment or update an existing one.
#' By default, it tries virtualenv first (simpler, more reliable), then falls back to conda.
#' 
#' @param conda_env The name of the environment to use. If NULL, uses the default from package configuration.
#' @param python_version The Python version to use. Default is "3.10".
#' @param pip_packages A character vector of pip packages to install.
#' @param method The method to use for environment setup. Options: "auto" (try virtualenv first, then conda), "virtualenv", "conda". Default is "auto".
#'
#' @return Invisible NULL. Called for side effects.
#' @export
setup_cassia_env <- function(conda_env = NULL, python_version = "3.10", 
                           pip_packages = c("openai", "pandas", "numpy", "scikit-learn", 
                                          "requests", "anthropic", "charset-normalizer"),
                           method = "auto") {
  
  # Get environment name
  env_name <- .get_env_name(conda_env)
  
  # Validate method parameter
  if (!method %in% c("auto", "virtualenv", "conda")) {
    stop("Method must be one of: 'auto', 'virtualenv', 'conda'")
  }
  
  # Check if environment already exists
  env_exists <- .check_env_exists(env_name)
  
  success <- FALSE
  
  if (method == "auto") {
    # Try virtualenv first, then conda
    if (env_exists$virtualenv) {
      message("Using existing virtualenv: ", env_name)
      success <- .try_virtualenv_setup(env_name, python_version, pip_packages)
    } else if (env_exists$conda) {
      message("Using existing conda environment: ", env_name)
      success <- .try_conda_setup(env_name, python_version, pip_packages)
    } else {
      # Neither exists, try to create virtualenv first
      message("Creating new environment. Trying virtualenv first...")
      success <- .try_virtualenv_setup(env_name, python_version, pip_packages)
      
      if (!success) {
        message("Virtualenv failed. Trying conda as fallback...")
        success <- .try_conda_setup(env_name, python_version, pip_packages)
      }
    }
  } else if (method == "virtualenv") {
    success <- .try_virtualenv_setup(env_name, python_version, pip_packages)
  } else if (method == "conda") {
    success <- .try_conda_setup(env_name, python_version, pip_packages)
  }
  
  if (!success) {
    stop("Failed to set up CASSIA Python environment with method: ", method)
  }
  
  # Maintain backward compatibility - keep the old conda_env option
  if (getOption("CASSIA.env_method", default = "conda") == "conda") {
    options(CASSIA.conda_env = env_name)
  }
  
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
  } else if (startsWith(provider, "http")) {
    # Handle custom HTTP endpoints
    Sys.setenv(CUSTERMIZED_API_KEY = api_key)
    env_var_name <- "CUSTERMIZED_API_KEY"
  } else {
    stop("Unsupported provider. Use 'openai', 'anthropic', 'openrouter', or an HTTP URL for custom endpoints")
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
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @return A list containing two elements: structured output and conversation history.
#' @export
runCASSIA <- function(model = "google/gemini-2.5-flash-preview", temperature, marker_list, tissue, species, additional_info = NULL, provider = "openrouter", validator_involvement = "v1") {
  tryCatch({
    result <- py_tools$runCASSIA(
      model = model,
      temperature = temperature,
      marker_list = marker_list,
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      provider = provider,
      validator_involvement = validator_involvement
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
    error_msg <- paste("Error in runCASSIA:", e$message, "\n",
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
#' @param temperature Numeric value for temperature parameter. Default: 0.3 (updated to match Python).
#' @param marker_list List of marker genes.
#' @param model Character string specifying the model to use.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#'
#' @return A list containing results from multiple runs.
#' @export
runCASSIA_n_times <- function(n, tissue, species, additional_info, temperature = 0.3, marker_list, 
                           model = "google/gemini-2.5-flash-preview", max_workers = 10, provider = "openrouter", validator_involvement = "v1") {
  tryCatch({
    result <- py_tools$runCASSIA_n_times(
      n = as.integer(n),
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      temperature = as.numeric(temperature),
      marker_list = marker_list,
      model = model,
      max_workers = as.integer(max_workers),
      provider = provider,
      validator_involvement = validator_involvement
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
    error_msg <- paste("Error in runCASSIA_n_times:", e$message, "\n",
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
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#'
#' @return A list containing processed results including variance analysis.
#' @export
runCASSIA_n_times_similarity_score <- function(tissue, species, additional_info, temperature, marker_list, model = "google/gemini-2.5-flash-preview", max_workers, n, provider = "openrouter", validator_involvement = "v1") {
  tryCatch({
    # Call the Python function with the new parameter structure
    processed_results <- py_uncertainty$runCASSIA_n_times_similarity_score(
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      temperature = as.numeric(temperature),
      marker_list = marker_list,
      model = model,
      max_workers = as.integer(max_workers),
      n = as.integer(n),
      provider=provider,
      validator_involvement = validator_involvement
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
    error_msg <- paste("Error in runCASSIA_n_times_similarity_score:", e$message, "\n",
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
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#'
#' @return None. This function creates output files and prints execution time.
#' @export
runCASSIA_batch <- function(marker, output_name = "cell_type_analysis_results.json", 
                          model = "google/gemini-2.5-flash-preview", temperature = 0, tissue = "lung", 
                          species = "human", additional_info = NULL, 
                          celltype_column = NULL, gene_column_name = NULL, 
                          max_workers = 10, provider = "openrouter", n_genes = 50,
                          max_retries = 1, validator_involvement = "v1") {
  execution_time <- system.time({
    # Convert R dataframe to Python if df_input is a dataframe
if (is.data.frame(marker)) {
  # Determine the cluster column to check
  cluster_col <- if (!is.null(celltype_column)) celltype_column else "cluster"

  # If the cluster column exists and is a factor, convert it to character
  if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
    marker[[cluster_col]] <- as.character(marker[[cluster_col]])
  }
  pd <- reticulate::import("pandas")
  # Direct conversion with convert=TRUE
  marker <- reticulate::r_to_py(marker, convert = TRUE)
} else if (!is.character(marker)) {
  stop("marker must be either a data frame or a character vector")
}
    
    py_tools$runCASSIA_batch(
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
      max_retries = as.integer(max_retries),
      validator_involvement = validator_involvement
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
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#'
#' @return None. This function creates output files and prints execution time.
#' @export
runCASSIA_batch_n_times <- function(n, marker, output_name = "cell_type_analysis_results", 
                                  model = "google/gemini-2.5-flash-preview", temperature = 0, tissue = "lung", 
                                  species = "human", additional_info = NULL, 
                                  celltype_column = NULL, gene_column_name = NULL, 
                                  max_workers = 10, batch_max_workers = 5, 
                                  provider = "openrouter", max_retries = 1, validator_involvement = "v1") {

  if (is.data.frame(marker)) {
    # Determine the cluster column to check
    cluster_col <- if (!is.null(celltype_column)) celltype_column else "cluster"

    # If the cluster column exists and is a factor, convert it to character
    if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
      marker[[cluster_col]] <- as.character(marker[[cluster_col]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }

  execution_time <- system.time({
    tryCatch({
      py_uncertainty$runCASSIA_batch_n_times(
        as.integer(n), marker, output_name, model, temperature, tissue, 
        species, additional_info, celltype_column, gene_column_name, 
        as.integer(max_workers), as.integer(batch_max_workers), provider,
        as.integer(max_retries), validator_involvement
      )
    }, error = function(e) {
      stop(paste("Error in runCASSIA_batch_n_times:", e$message))
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
    # Determine the cluster column to check
    cluster_col <- if (!is.null(celltype_column)) celltype_column else "cluster"

    # If the cluster column exists and is a factor, convert it to character
    if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
      marker[[cluster_col]] <- as.character(marker[[cluster_col]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }

  execution_time <- system.time({
    tryCatch({
      py_uncertainty$runCASSIA_similarity_score_batch(marker, file_pattern, output_name, 
                                              celltype_column, max_workers, model, provider, main_weight, sub_weight)
    }, error = function(e) {
      stop(paste("Error in runCASSIA_similarity_score_batch:", e$message))
    })
  })
  
  print(paste("Execution time for runCASSIA_batch_get_similarity_score:"))
  print(execution_time)
}







#' Generate Cell Type Analysis Report
#'
#' @param full_result_path Path to the full results CSV file
#' @param marker Path to the marker genes CSV file or data frame with marker data
#' @param cluster_name Name of the cluster to analyze
#' @param major_cluster_info General information about the dataset (e.g., "Human PBMC")
#' @param output_name Name of the output HTML file
#' @param num_iterations Number of iterations for marker analysis (default=5)
#' @param model Model to use for analysis (default="google/gemini-2.5-flash-preview")
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param temperature Sampling temperature (0-1)
#' @param conversation_history_mode Mode for extracting conversation history ("full", "final", or "none")
#' @param search_strategy Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time) (default: "breadth")
#' @param report_style Style of report ("per_iteration" or "total_summary") (default: "per_iteration")
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#'
#' @return None. This function generates output files.
#' @export
runCASSIA_annotationboost <- function(full_result_path, 
                                     marker, 
                                     cluster_name, 
                                     major_cluster_info, 
                                     output_name, 
                                     num_iterations = 5,
                                     model = "google/gemini-2.5-flash-preview",
                                     provider = "openrouter",
                                     temperature = 0,
                                     conversation_history_mode = "final",
                                     search_strategy = "breadth",
                                     report_style = "per_iteration",
                                     validator_involvement = "v1",
                                     ...) {

  if (is.data.frame(marker)) {
    # Factors can cause issues with reticulate, convert to character
    if ("cluster" %in% names(marker) && is.factor(marker[["cluster"]])) {
      marker[["cluster"]] <- as.character(marker[["cluster"]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }
                                                      
  tryCatch({
    py_annotation_boost$runCASSIA_annotationboost(
      full_result_path = full_result_path,
      marker = marker,
      cluster_name = cluster_name,
      major_cluster_info = major_cluster_info,
      output_name = output_name,
      num_iterations = as.integer(num_iterations),
      model = model,
      provider = provider,
      temperature = as.numeric(temperature),
      conversation_history_mode = conversation_history_mode,
      search_strategy = search_strategy,
      report_style = report_style,
      validator_involvement = validator_involvement
    )
    
    invisible(NULL)
    
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_annotationboost:", e$message)
    stop(error_msg)
  })
}




#' Generate Cell Type Analysis Report with Additional Task
#'
#' @param full_result_path Path to the full results CSV file
#' @param marker Path to the marker genes CSV file or data frame with marker data
#' @param cluster_name Name of the cluster to analyze
#' @param major_cluster_info General information about the dataset (e.g., "Human PBMC")
#' @param output_name Name of the output HTML file
#' @param num_iterations Number of iterations for marker analysis (default=5)
#' @param model Model to use for analysis (default="google/gemini-2.5-flash-preview")
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param additional_task Additional task to perform during analysis
#' @param temperature Sampling temperature (0-1)
#' @param conversation_history_mode Mode for extracting conversation history ("full", "final", or "none")
#' @param search_strategy Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time) (default: "breadth")
#' @param report_style Style of report ("per_iteration" or "total_summary") (default: "per_iteration")
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param ... Additional parameters for future compatibility
#' 
#' @return None. This function generates output files.
#' @export
runCASSIA_annotationboost_additional_task <- function(full_result_path, 
                                                     marker, 
                                                     cluster_name, 
                                                     major_cluster_info, 
                                                     output_name, 
                                                     num_iterations = 5,
                                                     model = "google/gemini-2.5-flash-preview",
                                                     provider = "openrouter",
                                                     additional_task = "check if this is a cancer cluster",
                                                     temperature = 0,
                                                     conversation_history_mode = "final",
                                                     search_strategy = "breadth",
                                                     report_style = "per_iteration",
                                                     validator_involvement = "v1",
                                                     ...) {

  if (is.data.frame(marker)) {
    # Factors can cause issues with reticulate, convert to character
    if ("cluster" %in% names(marker) && is.factor(marker[["cluster"]])) {
      marker[["cluster"]] <- as.character(marker[["cluster"]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }
                                                      
  tryCatch({
    py_annotation_boost$runCASSIA_annotationboost_additional_task(
      full_result_path = full_result_path,
      marker = marker,
      cluster_name = cluster_name,
      major_cluster_info = major_cluster_info,
      output_name = output_name,
      num_iterations = as.integer(num_iterations),
      model = model,
      provider = provider,
      additional_task = additional_task,
      temperature = as.numeric(temperature),
      conversation_history_mode = conversation_history_mode,
      search_strategy = search_strategy,
      report_style = report_style,
      validator_involvement = validator_involvement,
      ...
    )
    
    invisible(NULL)
    
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_annotationboost_additional_task:", e$message)
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
    results <- py_tools$runCASSIA_score_batch(
      input_file = input_file,
      output_file = output_file,
      max_workers = as.integer(max_workers),
      model = model,
      provider = provider,
      max_retries = as.integer(max_retries)
    )
    
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
    py_tools$runCASSIA_generate_score_report(
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
#' @param annotation_model Model to use for initial annotation (default: "google/gemini-2.5-flash-preview")
#' @param annotation_provider Provider for initial annotation (default: "openrouter")
#' @param score_model Model to use for scoring (default: "deepseek/deepseek-chat-v3-0324")
#' @param score_provider Provider for scoring (default: "openrouter")
#' @param annotationboost_model Model to use for boosting low-scoring annotations (default: "google/gemini-2.5-flash-preview")
#' @param annotationboost_provider Provider for boosting low-scoring annotations (default: "openrouter")
#' @param score_threshold Threshold for identifying low-scoring clusters (default: 75)
#' @param additional_info Additional information for analysis (default: NULL)
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#' @param do_merge_annotations Whether to run the merging annotations step (default: TRUE)
#' @param merge_model Model to use for merging annotations (default: "deepseek/deepseek-chat-v3-0324")
#' @param merge_provider Provider to use for merging annotations (default: "openrouter")
#' @param conversation_history_mode Mode for extracting conversation history ("full", "final", or "none") (default: "final")
#' @param search_strategy Search strategy for annotation boost - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time) (default: "breadth")
#' @param report_style Style of report for annotation boost ("per_iteration" or "total_summary") (default: "per_iteration")
#' @param ranking_method Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score') (default: "avg_log2FC")
#' @param ascending Sort direction (NULL uses default for each method) (default: NULL)
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#'
#' @return None. Creates output files and generates reports.
#' @export
runCASSIA_pipeline <- function(
    output_file_name,
    tissue,
    species,
    marker,
    max_workers = 4,
    annotation_model = "google/gemini-2.5-flash",
    annotation_provider = "openrouter",
    score_model = "google/gemini-2.5-flash",
    score_provider = "openrouter",
    annotationboost_model = "google/gemini-2.5-flash",
    annotationboost_provider = "openrouter",
    score_threshold = 75,
    additional_info = NULL,
    max_retries = 1,
    do_merge_annotations = TRUE,
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter",
    conversation_history_mode = "final",
    search_strategy = "breadth",
    report_style = "per_iteration",
    ranking_method = "avg_log2FC",
    ascending = NULL,
    validator_involvement = "v1"
) {
  # Convert R dataframe to Python if marker is a dataframe
  if (is.data.frame(marker)) {
    # Determine the cluster column to check
    cluster_col <- "cluster"
    if ("cluster" %in% names(marker)) {
      cluster_col <- "cluster"
    } else if (length(names(marker)) > 0) {
      cluster_col <- names(marker)[1]  # Use first column if no cluster column
    }

    # If the cluster column exists and is a factor, convert it to character
    if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
      marker[[cluster_col]] <- as.character(marker[[cluster_col]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }

  tryCatch({
    py_tools$runCASSIA_pipeline(
      output_file_name = output_file_name,
      tissue = tissue,
      species = species,
      marker = marker,
      max_workers = as.integer(max_workers),
      annotation_model = annotation_model,
      annotation_provider = annotation_provider,
      score_model = score_model,
      score_provider = score_provider,
      annotationboost_model = annotationboost_model,
      annotationboost_provider = annotationboost_provider,
      score_threshold = as.numeric(score_threshold),
      additional_info = if (is.null(additional_info)) "None" else additional_info,
      max_retries = as.integer(max_retries),
      merge_annotations = do_merge_annotations,
      merge_model = merge_model,
      merge_provider = merge_provider,
      conversation_history_mode = conversation_history_mode,
      ranking_method = ranking_method,
      ascending = ascending,
      report_style = report_style,
      validator_involvement = validator_involvement
    )
    
    invisible(NULL)
    
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_pipeline:", e$message)
    if (!is.null(reticulate::py_last_error())) {
      error_msg <- paste(error_msg, "\nPython traceback:", reticulate::py_last_error())
    }
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
    responses <- py_cell_comparison$compareCelltypes(
      tissue = tissue,
      celltypes = celltypes,
      marker_set = marker,
      species = species,
      model_list = model_list,
      output_file = output_file
    )

  }, error = function(e) {
    error_msg <- paste("Error in compareCelltypes:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}


#' Symphony Compare - Advanced Multi-Model Cell Type Comparison with Consensus Building
#'
#' Orchestrate multiple AI models to compare cell types with automatic consensus building.
#' This function conducts a comprehensive cell type comparison using multiple AI models in parallel,
#' automatically triggering discussion rounds when models disagree on the best matching cell type.
#' Think of it as a virtual panel of expert biologists debating and reaching consensus.
#'
#' @param tissue Character string specifying the tissue type (e.g., "blood", "brain", "liver")
#' @param celltypes Character vector of 2-4 cell types to compare
#' @param marker_set Character vector or string of gene markers to analyze
#' @param species Character string specifying the species (default: "human")
#' @param model_preset Character string specifying preset model configuration. Options:
#'   \itemize{
#'     \item "symphony": High-performance ensemble (Claude, GPT-4, Gemini Pro)
#'     \item "quartet": Balanced 4-model ensemble
#'     \item "budget": Cost-effective models
#'     \item "custom": Use custom_models list
#'   }
#' @param custom_models Character vector of custom models to use (when model_preset="custom")
#' @param output_dir Character string specifying directory to save results (default: current directory)
#' @param output_basename Character string for base name of output files (auto-generated if NULL)
#' @param enable_discussion Logical indicating whether to enable automatic discussion rounds when no consensus (default: TRUE)
#' @param max_discussion_rounds Integer specifying maximum discussion rounds to perform (default: 2)
#' @param consensus_threshold Numeric value (0-1) specifying fraction of models that must agree for consensus (default: 0.8)
#' @param generate_report Logical indicating whether to generate interactive HTML report (default: TRUE)
#' @param api_key Character string for OpenRouter API key (uses environment variable if NULL)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item results: List of all model responses and scores
#'     \item consensus: The consensus cell type (if reached)
#'     \item confidence: Confidence level of the consensus (0-1)
#'     \item csv_file: Path to the generated CSV file
#'     \item html_file: Path to the generated HTML report (if enabled)
#'     \item summary: Summary statistics of the comparison
#'     \item dataframe: R data frame with structured results
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage - let Symphony Compare handle everything
#' results <- symphonyCompare(
#'   tissue = "peripheral blood",
#'   celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
#'   marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
#'   species = "human"
#' )
#' 
#' # Access the results
#' cat("Consensus:", results$consensus, "\n")
#' cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\n")
#' 
#' # Advanced usage with custom settings
#' results <- symphonyCompare(
#'   tissue = "brain",
#'   celltypes = c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"),
#'   marker_set = c("RBFOX3", "GFAP", "IBA1", "OLIG2", "MAP2", "S100B", "CD11B", "MBP"),
#'   species = "mouse",
#'   model_preset = "quartet",  # Use 4 models instead of 3
#'   enable_discussion = TRUE,  # Enable automatic discussion rounds
#'   max_discussion_rounds = 3,  # Allow up to 3 discussion rounds
#'   consensus_threshold = 0.75,  # 75% of models must agree
#'   output_dir = "./symphony_results",
#'   verbose = TRUE
#' )
#' }
#'
#' @export
symphonyCompare <- function(tissue, celltypes, marker_set, species = "human", 
                          model_preset = "symphony", custom_models = NULL,
                          output_dir = NULL, output_basename = NULL,
                          enable_discussion = TRUE, max_discussion_rounds = 2L,
                          consensus_threshold = 0.8, generate_report = TRUE,
                          api_key = NULL, verbose = TRUE) {
  tryCatch({
    # Input validation
    if (length(celltypes) < 2 || length(celltypes) > 4) {
      stop("Please provide 2-4 cell types to compare")
    }
    
    # Convert marker_set to string if it's a vector
    if (is.vector(marker_set)) {
      marker_set <- paste(marker_set, collapse = ", ")
    }
    
    # Convert R types to Python-compatible types
    max_discussion_rounds <- as.integer(max_discussion_rounds)
    consensus_threshold <- as.numeric(consensus_threshold)
    enable_discussion <- as.logical(enable_discussion)
    generate_report <- as.logical(generate_report)
    verbose <- as.logical(verbose)
    
    # Call the Python function
    results <- py_symphony_compare$symphonyCompare(
      tissue = tissue,
      celltypes = celltypes,
      marker_set = marker_set,
      species = species,
      model_preset = model_preset,
      custom_models = custom_models,
      output_dir = output_dir,
      output_basename = output_basename,
      enable_discussion = enable_discussion,
      max_discussion_rounds = max_discussion_rounds,
      consensus_threshold = consensus_threshold,
      generate_report = generate_report,
      api_key = api_key,
      verbose = verbose
    )
    
    # Convert Python dataframe to R dataframe if it exists
    if (!is.null(results$dataframe)) {
      results$dataframe <- reticulate::py_to_r(results$dataframe)
    }
    
    # Return the results
    return(results)
    
  }, error = function(e) {
    error_msg <- paste("Error in symphonyCompare:", e$message, "\n",
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
  py_subclustering$runCASSIA_subclusters(
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
  py_subclustering$runCASSIA_n_subcluster(
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

#' Generate Subclustering HTML Report
#' 
#' Generate a beautiful HTML report for subclustering batch results, showing annotation, 
#' reasoning, and top markers in an interactive format.
#' 
#' @param csv_path Character string specifying the path to the subclustering CSV results file
#' @param html_report_path Character string specifying the output path for the HTML report. 
#'   If NULL, will use the same name as csv_path but with .html extension
#' @param model_name Character string specifying the model name to display in the report. 
#'   If NULL, will extract from the CSV filename
#' 
#' @return None. This function generates an HTML report file.
#' 
#' @examples
#' \dontrun{
#' # Generate a subclustering report from CSV results
#' generate_subclustering_report(
#'   csv_path = "subclustering_results.csv",
#'   html_report_path = "subclustering_report.html",
#'   model_name = "Gemini-2.5-Flash"
#' )
#' 
#' # Auto-generate HTML filename and model name
#' generate_subclustering_report("subclustering_results.csv")
#' }
#' 
#' @export
generate_subclustering_report <- function(csv_path, html_report_path = NULL, model_name = NULL) {
  tryCatch({
    py_generate_reports$generate_subclustering_report(
      csv_path = csv_path,
      html_report_path = html_report_path,
      model_name = model_name
    )
  }, error = function(e) {
    error_msg <- paste("Error in generate_subclustering_report:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}

#' Generate HTML Evaluation Report
#' 
#' Generate a comprehensive HTML report for evaluation results with visualizations,
#' metrics, and sample analyses.
#' 
#' @param result_df Data frame containing the evaluation results
#' @param gold_col Character string specifying the column name for ground truth/gold standard values
#' @param pred_col Character string specifying the column name for predicted values
#' @param score_col Character string specifying the column name for evaluation scores (default: "score")
#' @param reasoning_col Character string specifying the column name for reasoning/explanation text (default: "reasoning")
#' @param metrics Named list of pre-calculated metrics (optional, will be calculated if NULL)
#' @param html_report_path Character string specifying the output path for HTML report (default: "report.html")
#' @param model_name Character string specifying the model name to display in the report (optional)
#' 
#' @return None. This function generates an HTML report file.
#' 
#' @examples
#' \dontrun{
#' # Generate evaluation report
#' generate_html_report(
#'   result_df = evaluation_results,
#'   gold_col = "True Cell Type",
#'   pred_col = "Predicted Cell Type",
#'   score_col = "evaluation_score",
#'   reasoning_col = "explanation",
#'   html_report_path = "evaluation_report.html",
#'   model_name = "Claude-3.5-Sonnet"
#' )
#' }
#' 
#' @export
generate_html_report <- function(result_df, gold_col, pred_col, score_col = "score", 
                                reasoning_col = "reasoning", metrics = NULL, 
                                html_report_path = "report.html", model_name = NULL) {
  tryCatch({
    py_generate_reports$generate_html_report(
      result_df = result_df,
      gold_col = gold_col,
      pred_col = pred_col,
      score_col = score_col,
      reasoning_col = reasoning_col,
      metrics = metrics,
      html_report_path = html_report_path,
      model_name = model_name
    )
  }, error = function(e) {
    error_msg <- paste("Error in generate_html_report:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}

#' Calculate Evaluation Metrics
#' 
#' Calculate comprehensive metrics from evaluation results including mean scores,
#' distributions, and performance ratios.
#' 
#' @param eval_df Data frame containing evaluation results
#' @param score_col Character string specifying the column name for evaluation scores (default: "score")
#' 
#' @return A named list containing calculated metrics:
#'   \itemize{
#'     \item mean_score: Average score
#'     \item median_score: Median score  
#'     \item min_score: Minimum score
#'     \item max_score: Maximum score
#'     \item std_score: Standard deviation of scores
#'     \item count: Number of evaluations
#'     \item For 0-5 scale: perfect_ratio, very_good_ratio, good_ratio, partial_ratio, poor_ratio, nonsensical_ratio
#'   }
#' 
#' @examples
#' \dontrun{
#' # Calculate metrics from evaluation results
#' metrics <- calculate_evaluation_metrics(
#'   eval_df = evaluation_results,
#'   score_col = "evaluation_score"
#' )
#' 
#' # Print key metrics
#' cat("Mean Score:", metrics$mean_score, "\n")
#' cat("Perfect Predictions:", sprintf("%.1f%%", metrics$perfect_ratio * 100), "\n")
#' }
#' 
#' @export
calculate_evaluation_metrics <- function(eval_df, score_col = "score") {
  tryCatch({
    result <- py_generate_reports$calculate_evaluation_metrics(
      eval_df = eval_df,
      score_col = score_col
    )
    # Convert Python dict to R list
    return(reticulate::py_to_r(result))
  }, error = function(e) {
    error_msg <- paste("Error in calculate_evaluation_metrics:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}