# Add any additional utility functions here
#' Check Python Environment
#'
#' @return TRUE if the Python environment is set up correctly, FALSE otherwise.
#' @export
check_python_env <- function() {
  tryCatch({
    py_main <- reticulate::import_from_path("main_function_code", path = system.file("python", package = "CASSIA"))
    py_tools <- reticulate::import_from_path("tools_function", path = system.file("python", package = "CASSIA"))
    return(TRUE)
  }, error = function(e) {
    warning(paste("Python environment not set up correctly:", e$message))
    return(FALSE)
  })
}

#' Set Python Environment
#'
#' @param conda_env The name of the conda environment to use.
#' @return TRUE if the Python environment is set successfully, FALSE otherwise.
#' @export
set_python_env <- function(conda_env) {
  tryCatch({
    reticulate::use_condaenv(conda_env, required = TRUE)
    options(CASSIA.conda_env = conda_env)
    return(TRUE)
  }, error = function(e) {
    warning(paste("Failed to set Python environment:", e$message))
    return(FALSE)
  })
}
