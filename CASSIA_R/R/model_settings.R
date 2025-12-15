# Model Settings Functions for CASSIA R Package

#' Get Python CASSIA Module
#'
#' @return Python CASSIA module (all functions available at package level)
#' @keywords internal
.get_py_model_settings <- function() {
  # Use the py_cassia from package namespace if available
  # All model settings functions are exported at the CASSIA package level

  # First check if py_cassia exists in the package namespace
  pkg_env <- tryCatch(asNamespace("CASSIA"), error = function(e) NULL)
  if (!is.null(pkg_env) && exists("py_cassia", envir = pkg_env)) {
    py_obj <- get("py_cassia", envir = pkg_env)
    if (!is.null(py_obj)) {
      return(py_obj)
    }
  }

  # Also check parent frames (for interactive use / devtools::load_all)
  py_obj <- get0("py_cassia", envir = parent.frame(), ifnotfound = NULL)
  if (!is.null(py_obj)) {
    return(py_obj)
  }

  # Check in global environment as fallback
  if (exists("py_cassia", envir = .GlobalEnv)) {
    py_obj <- get("py_cassia", envir = .GlobalEnv)
    if (!is.null(py_obj)) {
      return(py_obj)
    }
  }

  # Otherwise, try to import directly
  tryCatch({
    py_cassia_local <- reticulate::import_from_path("CASSIA", path = system.file("python", package = "CASSIA"))
    return(py_cassia_local)
  }, error = function(e) {
    warning("Could not load CASSIA Python package: ", e$message)
    return(NULL)
  })
}

#' Resolve Model Name
#'
#' Resolve a user-provided model name to the actual provider model name.
#' This function automatically maps simplified names (like "gpt4", "claude", "gemini") 
#' to their actual provider-specific names.
#'
#' @param model_name Character string of the model name (can be alias or actual name)
#' @param provider Character string of the provider name (if NULL, will try to infer from model_name)
#'
#' @return A list with two elements: 'model' (actual model name) and 'provider' (provider name)
#' @export
#'
#' @examples
#' \dontrun{
#' # Using simple names
#' resolve_model_name("gpt4")  # Returns gpt-4o with openai provider
#' resolve_model_name("claude")  # Returns claude-3-5-sonnet-latest with anthropic provider
#' resolve_model_name("gemini")  # Returns google/gemini-2.5-flash with openrouter provider
#' 
#' # Using aliases
#' resolve_model_name("sonnet")  # Returns claude-3-5-sonnet-latest
#' resolve_model_name("deepseek")  # Returns deepseek/deepseek-chat-v3-0324
#' 
#' # With specific provider
#' resolve_model_name("gpt-4o", "openai")
#' }
resolve_model_name <- function(model_name, provider = NULL) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    result <- py_settings$resolve_model_name(model_name, provider)
    
    # Convert Python tuple to R list
    return(list(
      model = result[[1]],
      provider = result[[2]]
    ))
  }, error = function(e) {
    warning("Error resolving model name: ", e$message)
    return(list(model = model_name, provider = provider))
  })
}

#' Get Recommended Model
#'
#' Get the recommended model for a provider or use case.
#'
#' @param provider Character string of the provider name (if NULL, returns overall best model)
#' @param use_case Character string of the use case (annotation, scoring, annotation_boost, merging)
#'
#' @return A list with two elements: 'model' (recommended model name) and 'provider' (provider name)
#' @export
#'
#' @examples
#' \dontrun{
#' # Get overall best model
#' get_recommended_model()
#' 
#' # Get best model for specific use case
#' get_recommended_model(use_case = "annotation")
#' get_recommended_model(use_case = "scoring")
#' get_recommended_model(use_case = "annotation_boost")
#' 
#' # Get recommended model for specific provider
#' get_recommended_model(provider = "openai")
#' get_recommended_model(provider = "anthropic")
#' }
get_recommended_model <- function(provider = NULL, use_case = NULL) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    result <- py_settings$get_recommended_model(provider, use_case)
    
    # Convert Python tuple to R list
    return(list(
      model = result[[1]],
      provider = result[[2]]
    ))
  }, error = function(e) {
    warning("Error getting recommended model: ", e$message)
    return(list(model = "google/gemini-2.5-flash", provider = "openrouter"))
  })
}

#' Get Model Information
#'
#' Get detailed information about a model.
#'
#' @param model_name Character string of the model name
#' @param provider Character string of the provider name (if NULL, will try to infer)
#'
#' @return A list with model information or NULL if not found
#' @export
#'
#' @examples
#' \dontrun{
#' # Get info for specific models
#' get_model_info("gpt-4o")
#' get_model_info("claude-3-5-sonnet-latest")
#' get_model_info("google/gemini-2.5-flash")
#' 
#' # Using aliases
#' get_model_info("gemini")
#' get_model_info("deepseek")
#' }
get_model_info <- function(model_name, provider = NULL) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    result <- py_settings$get_model_info(model_name, provider)
    
    if (is.null(result)) {
      return(NULL)
    }
    
    # Convert Python dict to R list
    return(reticulate::py_to_r(result))
  }, error = function(e) {
    warning("Error getting model info: ", e$message)
    return(NULL)
  })
}

#' List Available Models
#'
#' List available models with optional filtering.
#'
#' @param provider Character string to filter by provider
#' @param cost_tier Character string to filter by cost tier (very_low, low, medium, high)
#' @param use_case Character string to filter by use case
#'
#' @return A data frame with model information
#' @export
#'
#' @examples
#' \dontrun{
#' # List all models
#' list_models()
#' 
#' # Filter by provider
#' list_models(provider = "openai")
#' list_models(provider = "anthropic")
#' list_models(provider = "openrouter")
#' 
#' # Filter by cost tier
#' list_models(cost_tier = "low")
#' list_models(cost_tier = "high")
#' 
#' # Filter by use case
#' list_models(use_case = "annotation")
#' list_models(use_case = "scoring")
#' 
#' # Combine filters
#' list_models(provider = "openrouter", cost_tier = "low")
#' }
list_models <- function(provider = NULL, cost_tier = NULL, use_case = NULL) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    result <- py_settings$list_models(provider, cost_tier, use_case)
    
    # Convert Python list of dicts to R data frame
    if (length(result) == 0) {
      return(data.frame(
        name = character(0),
        provider = character(0),
        description = character(0),
        cost_tier = character(0),
        recommended = logical(0),
        stringsAsFactors = FALSE
      ))
    }
    
    # Extract information from Python result
    df <- data.frame(
      name = sapply(result, function(x) x$name),
      provider = sapply(result, function(x) x$provider),
      description = sapply(result, function(x) x$info$description),
      cost_tier = sapply(result, function(x) x$info$cost_tier),
      recommended = sapply(result, function(x) x$info$recommended),
      stringsAsFactors = FALSE
    )
    
    return(df)
  }, error = function(e) {
    warning("Error listing models: ", e$message)
    return(data.frame())
  })
}

#' Get Use Case Recommendations
#'
#' Get recommendations for a specific use case.
#'
#' @param use_case Character string of the use case name
#'
#' @return A list with recommendations
#' @export
#'
#' @examples
#' \dontrun{
#' # Get recommendations for different use cases
#' get_use_case_recommendations("annotation")
#' get_use_case_recommendations("scoring")
#' get_use_case_recommendations("annotation_boost")
#' get_use_case_recommendations("merging")
#' }
get_use_case_recommendations <- function(use_case) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    result <- py_settings$get_use_case_recommendations(use_case)
    
    # Convert Python dict to R list
    return(reticulate::py_to_r(result))
  }, error = function(e) {
    warning("Error getting use case recommendations: ", e$message)
    return(list())
  })
}

#' Print Model Recommendations
#'
#' Print model recommendations in a user-friendly format.
#'
#' @param use_case Character string of the specific use case to show recommendations for (optional)
#'
#' @return None (prints to console)
#' @export
#'
#' @examples
#' \dontrun{
#' # Print all recommendations
#' print_model_recommendations()
#' 
#' # Print recommendations for specific use case
#' print_model_recommendations("annotation")
#' print_model_recommendations("scoring")
#' }
print_model_recommendations <- function(use_case = NULL) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    py_settings$print_model_recommendations(use_case)
  }, error = function(e) {
    warning("Error printing model recommendations: ", e$message)
  })
}

#' Get Model Aliases
#'
#' Get all aliases for a model.
#'
#' @param model_name Character string of the model name
#' @param provider Character string of the provider name (if NULL, will try to infer)
#'
#' @return A character vector of aliases
#' @export
#'
#' @examples
#' \dontrun{
#' # Get aliases for models
#' get_model_aliases("gpt-4o")
#' get_model_aliases("claude-3-5-sonnet-latest")
#' get_model_aliases("google/gemini-2.5-flash")
#' }
get_model_aliases <- function(model_name, provider = NULL) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    result <- py_settings$get_model_aliases(model_name, provider)
    
    # Convert Python list to R vector
    return(as.character(result))
  }, error = function(e) {
    warning("Error getting model aliases: ", e$message)
    return(character(0))
  })
}

#' Get Cost Tier Models
#'
#' Get models in a specific cost tier.
#'
#' @param cost_tier Character string of the cost tier name
#'
#' @return A character vector of model names
#' @export
#'
#' @examples
#' \dontrun{
#' # Get models by cost tier
#' get_cost_tier_models("very_low")
#' get_cost_tier_models("low")
#' get_cost_tier_models("medium")
#' get_cost_tier_models("high")
#' }
get_cost_tier_models <- function(cost_tier) {
  tryCatch({
    py_settings <- .get_py_model_settings()
    result <- py_settings$get_cost_tier_models(cost_tier)
    
    # Convert Python list to R vector
    return(as.character(result))
  }, error = function(e) {
    warning("Error getting cost tier models: ", e$message)
    return(character(0))
  })
}