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

#' Add CASSIA Annotations to Seurat Object
#'
#' This function integrates CASSIA annotation results with a Seurat object by adding
#' the cell type annotations as metadata columns. The function matches cluster IDs 
#' from CASSIA results to the cluster IDs in the Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param cassia_results_path Path to the CASSIA output CSV file
#' @param cluster_col Name of the column in Seurat metadata containing cluster IDs (default: "seurat_clusters")
#' @param cassia_cluster_col Name of the column in CASSIA results containing cluster IDs (default: "True Cell Type")
#' @param prefix Prefix to add to the new metadata columns (default: "CASSIA_")
#' @param replace_existing Whether to replace existing annotations (default: FALSE)
#' @param fuzzy_match Whether to perform fuzzy matching on cluster names (default: TRUE)
#'
#' @return A Seurat object with CASSIA annotations added as metadata columns
#' @importFrom utils read.csv
#' @importFrom Seurat AddMetaData
#' @export
add_cassia_to_seurat <- function(seurat_obj, cassia_results_path, cluster_col = "seurat_clusters", 
                                cassia_cluster_col = "True Cell Type", prefix = "CASSIA_", 
                                replace_existing = FALSE, fuzzy_match = TRUE) {
  # Check if required packages are installed
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required. Please install it with 'install.packages(\"Seurat\")'")
  }
  
  # Check if file exists
  if (!file.exists(cassia_results_path)) {
    stop("CASSIA results file not found at: ", cassia_results_path)
  }
  
  # Check if the specified cluster column exists in the Seurat object
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found in Seurat metadata.")
  }
  
  # Read CASSIA results
  cassia_results <- read.csv(cassia_results_path, stringsAsFactors = FALSE)
  
  # Fix column names (replace dots with spaces)
  names(cassia_results) <- gsub("\\.", " ", names(cassia_results))
  if (grepl("\\.", cassia_cluster_col)) {
    cassia_cluster_col <- gsub("\\.", " ", cassia_cluster_col)
  }
  
  # Check if the cluster column exists in the CASSIA results
  if (!cassia_cluster_col %in% colnames(cassia_results)) {
    # Check for common alternative names
    alt_cluster_cols <- c("True Cell Type", "Cluster", "cluster_id", "ClusterID", "Cluster_ID", "cluster")
    found_cols <- alt_cluster_cols[alt_cluster_cols %in% colnames(cassia_results)]
    
    if (length(found_cols) == 0) {
      message("Available columns in CASSIA results: ", paste(colnames(cassia_results), collapse = ", "))
      stop("Could not find cluster column in CASSIA results. Please specify the correct column name using 'cassia_cluster_col'.")
    } else {
      cassia_cluster_col <- found_cols[1]
      message("Using '", cassia_cluster_col, "' as the cluster column from CASSIA results.")
    }
  }
  
  # Get unique cluster values
  seurat_clusters_unique <- unique(as.character(seurat_obj@meta.data[[cluster_col]]))
  cassia_clusters_unique <- unique(as.character(cassia_results[[cassia_cluster_col]]))
  
  # Create a mapping between CASSIA clusters and Seurat clusters if needed
  cluster_map <- NULL
  
  if (fuzzy_match && !all(seurat_clusters_unique %in% cassia_clusters_unique)) {
    # Create a normalized version of both cluster sets for comparison
    normalize_text <- function(x) {
      x <- tolower(x)
      x <- gsub("[[:punct:]]", "", x)
      x <- gsub("\\s+", " ", x)
      x <- trimws(x)
      return(x)
    }
    
    seurat_norm <- normalize_text(seurat_clusters_unique)
    cassia_norm <- normalize_text(cassia_clusters_unique)
    
    # Check if normalized versions match
    cluster_map <- data.frame(
      seurat_cluster = seurat_clusters_unique,
      cassia_cluster = NA_character_,
      stringsAsFactors = FALSE
    )
    
    # Direct match after normalization
    for (i in seq_along(seurat_norm)) {
      matches <- which(cassia_norm == seurat_norm[i])
      if (length(matches) == 1) {
        cluster_map$cassia_cluster[i] <- cassia_clusters_unique[matches]
      }
    }
    
    # For remaining unmatched clusters, try substring matching
    for (i in which(is.na(cluster_map$cassia_cluster))) {
      for (j in seq_along(cassia_norm)) {
        if (grepl(seurat_norm[i], cassia_norm[j]) || grepl(cassia_norm[j], seurat_norm[i])) {
          cluster_map$cassia_cluster[i] <- cassia_clusters_unique[j]
          break
        }
      }
    }
    
    # Print the mapping if any clusters were mapped
    if (any(!is.na(cluster_map$cassia_cluster))) {
      message("Created mapping between Seurat clusters and CASSIA clusters:")
      for (i in which(!is.na(cluster_map$cassia_cluster))) {
        if (cluster_map$seurat_cluster[i] != cluster_map$cassia_cluster[i]) {
          message("  '", cluster_map$seurat_cluster[i], "' -> '", cluster_map$cassia_cluster[i], "'")
        }
      }
    }
    
    # Check for any unmatched clusters
    unmatched <- which(is.na(cluster_map$cassia_cluster))
    if (length(unmatched) > 0) {
      warning("Could not find matches for these Seurat clusters: ", 
              paste(cluster_map$seurat_cluster[unmatched], collapse = ", "))
    }
  }
  
  # Define the column mappings for CASSIA output
  column_mapping <- list(
    general_celltype = c("Predicted Main Cell Type", "General Cell Type", "General_Cell_Type"),
    sub_celltype = c("Predicted Sub Cell Types", "Sub Cell Type", "Sub_Cell_Type"),
    mixed_celltype = c("Possible Mixed Cell Types", "Mixed Cell Types", "Possible_Mixed_Cell_Types"),
    score = c("Score", "Consensus Score", "Consensus_Score")
  )
  
  # Create a lookup table from CASSIA clusters to annotations
  lookup_tables <- list()
  
  # Also create lookup tables for split sub cell types
  lookup_tables_split <- list(
    sub_celltype_all = NULL,       # Original full string
    sub_celltype_1 = NULL,         # Most likely
    sub_celltype_2 = NULL,         # Less likely
    sub_celltype_3 = NULL          # Least likely
  )
  
  for (anno_type in names(column_mapping)) {
    col_candidates <- column_mapping[[anno_type]]
    found_col <- col_candidates[col_candidates %in% colnames(cassia_results)]
    
    if (length(found_col) > 0) {
      col_name <- found_col[1]
      
      # Special handling for sub cell types - keep the original column but also split it
      if (anno_type == "sub_celltype") {
        # Store the original values
        sub_types_full <- cassia_results[[col_name]]
        lookup_tables_split$sub_celltype_all <- setNames(
          sub_types_full,
          cassia_results[[cassia_cluster_col]]
        )
        
        # Split the comma-separated values into individual entries
        split_sub_types <- lapply(sub_types_full, function(x) {
          if (is.na(x) || x == "") {
            return(c(NA_character_, NA_character_, NA_character_))
          }
          
          # Remove any brackets, clean up the string
          clean_x <- gsub("^\\s*\\[?\\s*|\\s*\\]?\\s*$", "", x)
          
          # Split by comma and trim whitespace
          parts <- trimws(strsplit(clean_x, ",")[[1]])
          
          # Return three parts (padded with NA if fewer)
          if (length(parts) >= 3) {
            return(parts[1:3])
          } else {
            return(c(parts, rep(NA_character_, 3 - length(parts))))
          }
        })
        
        # Create lookup tables for each position
        lookup_tables_split$sub_celltype_1 <- setNames(
          sapply(split_sub_types, function(x) x[1]),
          cassia_results[[cassia_cluster_col]]
        )
        
        lookup_tables_split$sub_celltype_2 <- setNames(
          sapply(split_sub_types, function(x) x[2]),
          cassia_results[[cassia_cluster_col]]
        )
        
        lookup_tables_split$sub_celltype_3 <- setNames(
          sapply(split_sub_types, function(x) x[3]),
          cassia_results[[cassia_cluster_col]]
        )
        
        # For the regular lookup table, just use the first (most likely) subtype
        lookup_tables[[anno_type]] <- lookup_tables_split$sub_celltype_1
      } else {
        # For other annotation types, just use as is
        lookup_tables[[anno_type]] <- setNames(
          cassia_results[[col_name]],
          cassia_results[[cassia_cluster_col]]
        )
      }
    }
  }
  
  # Get cluster IDs from Seurat object
  seurat_clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  
  # Apply cluster mapping if it exists
  if (!is.null(cluster_map) && any(!is.na(cluster_map$cassia_cluster))) {
    # Create a mapping function
    map_cluster <- function(x) {
      idx <- match(x, cluster_map$seurat_cluster)
      if (is.na(idx) || is.na(cluster_map$cassia_cluster[idx])) {
        return(x)
      } else {
        return(cluster_map$cassia_cluster[idx])
      }
    }
    
    # Map the clusters
    seurat_clusters_mapped <- sapply(seurat_clusters, map_cluster)
  } else {
    seurat_clusters_mapped <- seurat_clusters
  }
  
  # Create a new data frame for metadata
  new_metadata <- data.frame(row.names = colnames(seurat_obj))
  
  # For each annotation type, map cluster IDs to annotations
  for (anno_type in names(lookup_tables)) {
    # Get the lookup table
    lookup <- lookup_tables[[anno_type]]
    
    # Map cluster IDs to annotations
    annotations <- lookup[seurat_clusters_mapped]
    
    # Add to the new metadata data frame
    meta_col_name <- paste0(prefix, anno_type)
    new_metadata[[meta_col_name]] <- annotations
    
    # Optionally replace existing annotations
    if (replace_existing && anno_type == "general_celltype") {
      if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
        message("Replacing 'cell_type' column with CASSIA annotations")
        seurat_obj$cell_type <- annotations
      }
    }
  }
  
  # For sub cell types, add all the split versions
  for (split_type in names(lookup_tables_split)) {
    if (!is.null(lookup_tables_split[[split_type]])) {
      lookup <- lookup_tables_split[[split_type]]
      
      # Map cluster IDs to annotations
      annotations <- lookup[seurat_clusters_mapped]
      
      # Add to the metadata
      meta_col_name <- paste0(prefix, split_type)
      new_metadata[[meta_col_name]] <- annotations
    }
  }
  
  # Add all metadata columns at once
  seurat_obj <- Seurat::AddMetaData(seurat_obj, new_metadata)
  
  message("Added CASSIA annotations to Seurat object")
  message("New metadata columns: ", paste(names(new_metadata), collapse = ", "))
  
  return(seurat_obj)
}
