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


loadExampleMarkers_subcluster <- function(){
  return(read.csv(system.file("extdata","subcluster_results.csv",package = "CASSIA")))
}
