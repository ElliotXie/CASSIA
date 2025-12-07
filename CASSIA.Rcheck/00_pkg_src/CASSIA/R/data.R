#' Load Example Marker Data
#'
#' This function loads an example marker dataset that can be used to test CASSIA's
#' cell type annotation functionality. The dataset contains marker genes for
#' different immune cell types from a PBMC dataset.
#'
#' @param processed Logical. If TRUE, loads preprocessed RDS data. If FALSE, loads raw CSV data.
#' @export
#'
#' @examples
#' markers <- loadExampleMarkers()
#' head(markers)
loadExampleMarkers <- function(processed = FALSE) {
  file_path <- system.file(
    "extdata", 
    if (processed) "processed.csv" else "unprocessed.RDS",
    package = "CASSIA"
  )
  
  if (file_path == "") {
    stop("Example data files not found. Please ensure CASSIA is properly installed.")
  }
  
  markers <- if (processed) {
    read.csv(file_path)
  } else {
    readRDS(file_path)
  }
  
  return(markers)
}


#' Load Example Subcluster Results
#'
#' This function loads example subcluster results data that can be used to test CASSIA's
#' functionality.
#'
#' @return A data frame containing subcluster results
#' @export
#'
#' @examples
#' subcluster_data <- loadExampleMarkers_subcluster()
#' head(subcluster_data)
loadExampleMarkers_subcluster <- function() {
  file_path <- system.file("extdata", "subcluster_results.csv", package = "CASSIA")
  
  if (file_path == "") {
    stop("Example data files not found. Please ensure CASSIA is properly installed.")
  }
  
  return(read.csv(file_path))
}

#' Load Built-in Marker Data
#'
#' This function loads built-in marker files from the CASSIA package.
#'
#' @param marker_type Type of markers to load. Options:
#'   - "processed": For processed marker data (default)
#'   - "unprocessed": For raw unprocessed marker data
#'   - "subcluster_results": For subcluster analysis results
#' @return A data frame containing marker data
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- loadBuiltinMarkers()
#' head(markers)
#'
#' subcluster_results <- loadBuiltinMarkers("subcluster_results")
#' head(subcluster_results)
#' }
loadBuiltinMarkers <- function(marker_type = "processed") {
  # Import the Python function
  # tools_function <- reticulate::import_from_path("tools_function", 
  #                                                path = system.file("python", 
  #                                                                   package = "CASSIA"))
  
  # Call the Python function
  markers <- py_cassia$loadmarker(marker_type = marker_type)
  
  # Convert to R data frame if necessary
  if (!is.data.frame(markers)) {
    markers <- as.data.frame(markers)
  }
  
  return(markers)
}

#' List Available Built-in Marker Sets
#'
#' This function lists all available built-in marker sets in the CASSIA package.
#'
#' @return A character vector containing names of available marker files
#' @export
#'
#' @examples
#' \dontrun{
#' available_markers <- listAvailableMarkers()
#' print(available_markers)
#' }
listAvailableMarkers <- function() {
  # Import the Python function
  # tools_function <- reticulate::import_from_path("tools_function", 
  #                                                path = system.file("python", 
  #                                                                   package = "CASSIA"))
  
  # Call the Python function and convert to R character vector
  marker_names <- py_cassia$list_available_markers()
  
  return(marker_names)
}
