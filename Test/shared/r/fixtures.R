# CASSIA Test Suite - R Data Fixtures
# ====================================
# Functions for loading and preparing test data in R
#
# Note: This file requires test_utils.R to be sourced first (for get_test_root())

#' Get path to marker data file
#'
#' @return Character path to processed.csv
get_marker_file_path <- function() {
  file.path(get_test_root(), "data", "markers", "processed.csv")
}

#' Load the processed marker data
#'
#' @return Data frame with columns Broad.cell.type and Top.Markers
load_markers <- function() {
  marker_path <- get_marker_file_path()

  if (!file.exists(marker_path)) {
    stop(paste("Marker file not found:", marker_path))
  }

  read.csv(marker_path, row.names = 1, stringsAsFactors = FALSE)
}

#' Get list of all available cell type clusters
#'
#' @return Character vector of cluster names
get_all_clusters <- function() {
  df <- load_markers()
  df$Broad.cell.type
}

#' Get marker genes for a specific cluster
#'
#' @param cluster_name Name of the cell type cluster
#' @param n_genes Number of top genes to return (NULL = all)
#' @return Character vector of marker gene names
get_cluster_markers <- function(cluster_name, n_genes = NULL) {
  df <- load_markers()

  # Find the cluster
  cluster_row <- df[df$Broad.cell.type == cluster_name, ]

  if (nrow(cluster_row) == 0) {
    available <- get_all_clusters()
    stop(paste0("Cluster '", cluster_name, "' not found. Available: ",
                paste(available, collapse = ", ")))
  }

  # Get markers string and split
  markers_str <- cluster_row$Top.Markers[1]
  markers <- trimws(strsplit(markers_str, ",")[[1]])

  if (!is.null(n_genes)) {
    markers <- head(markers, n_genes)
  }

  return(markers)
}

#' Get a data frame formatted for runCASSIA for a single cluster
#'
#' @param cluster_name Name of the cell type cluster
#' @param n_genes Number of top genes to include
#' @return Data frame with columns gene and cell_type
get_marker_dataframe_for_cluster <- function(cluster_name, n_genes = 30) {
  markers <- get_cluster_markers(cluster_name, n_genes)
  data.frame(
    gene = markers,
    cell_type = rep(cluster_name, length(markers)),
    stringsAsFactors = FALSE
  )
}

#' Get the full marker data frame for batch processing
#'
#' @return Data frame with all clusters and their markers
get_full_marker_dataframe <- function() {
  load_markers()
}
