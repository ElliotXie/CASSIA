# Add any additional utility functions here

#' Validate API Keys
#'
#' Validate API keys by making a minimal test call to the provider(s).
#' Uses intelligent caching to avoid redundant API calls. Each unique key is
#' only tested once per session unless the key changes or force_revalidate=TRUE.
#'
#' The validation makes a minimal test call using the cheapest model for each
#' provider (costing ~$0.000001 per validation):
#' \itemize{
#'   \item OpenAI: gpt-4o-mini
#'   \item Anthropic: claude-3-haiku-20240307
#'   \item OpenRouter: openai/gpt-4o-mini
#' }
#'
#' @param provider Specific provider to validate ('openai', 'anthropic', 'openrouter').
#'   If NULL, validates all providers with keys set in environment variables.
#' @param api_key Specific API key to validate. If NULL, reads from environment variables.
#' @param force_revalidate Force revalidation even if key was previously validated. Default: FALSE.
#' @param verbose Print validation status and results. Default: TRUE.
#'
#' @return If provider specified: logical (TRUE if valid, FALSE if invalid).
#'   If provider=NULL: named list mapping provider names to logical values.
#'
#' @examples
#' \dontrun{
#' # Validate all configured providers
#' validate_api_keys()
#'
#' # Validate specific provider
#' validate_api_keys("openai")
#'
#' # Validate a specific key without setting it
#' validate_api_keys("openai", api_key = "sk-test...")
#'
#' # Force revalidation (skip cache)
#' validate_api_keys("openai", force_revalidate = TRUE)
#' }
#'
#' @export
validate_api_keys <- function(provider = NULL, api_key = NULL,
                              force_revalidate = FALSE, verbose = TRUE) {
  # Import CASSIA Python module
  cassia <- tryCatch({
    reticulate::import_from_path("CASSIA", path = system.file("python", package = "CASSIA"))
  }, error = function(e) {
    stop("Failed to import CASSIA Python module: ", e$message)
  })

  # Call the Python validate_api_keys function
  result <- cassia$validate_api_keys(
    provider = provider,
    api_key = api_key,
    force_revalidate = force_revalidate,
    verbose = verbose
  )

  # Convert Python dict to R named list if needed
  if (is.null(provider) && inherits(result, "python.builtin.dict")) {
    result <- reticulate::py_to_r(result)
  }

  return(result)
}

#' Clear API Key Validation Cache
#'
#' Clear the API key validation cache. This can be useful if:
#' \itemize{
#'   \item You've rotated your API keys and want to force revalidation
#'   \item You suspect the cache is stale
#'   \item You're debugging validation issues
#' }
#'
#' @param provider Specific provider to clear from cache. If NULL, clears
#'   entire cache for all providers.
#'
#' @examples
#' \dontrun{
#' # Clear cache for specific provider
#' clear_validation_cache("openai")
#'
#' # Clear entire cache
#' clear_validation_cache()
#' }
#'
#' @export
clear_validation_cache <- function(provider = NULL) {
  # Import CASSIA Python module
  cassia <- tryCatch({
    reticulate::import_from_path("CASSIA", path = system.file("python", package = "CASSIA"))
  }, error = function(e) {
    stop("Failed to import CASSIA Python module: ", e$message)
  })

  # Call the Python clear_validation_cache function
  cassia$clear_validation_cache(provider = provider)

  invisible(NULL)
}

#' Check Python Environment
#'
#' @return TRUE if the Python environment is set up correctly, FALSE otherwise.
#' @importFrom stats setNames
#' @export
check_python_env <- function() {
  tryCatch({
    # Test import of CASSIA Python package
    test_cassia <- reticulate::import_from_path("CASSIA", path = system.file("python", package = "CASSIA"))
    return(TRUE)
  }, error = function(e) {
    warning(paste("Python environment not set up correctly:", e$message))
    return(FALSE)
  })
}

#' Set Python Environment
#'
#' @param conda_env The name of the environment to use (works with both conda and virtualenv).
#' @return TRUE if the Python environment is set successfully, FALSE otherwise.
#' @export
set_python_env <- function(conda_env) {
  # Maintain backward compatibility - support old parameter name
  env_name <- conda_env
  tryCatch({
    # Check if it's a virtualenv first
    virtualenv_exists <- tryCatch({
      existing_virtualenvs <- reticulate::virtualenv_list()
      env_name %in% existing_virtualenvs
    }, error = function(e) FALSE)
    
    # Check if it's a conda environment
    conda_exists <- tryCatch({
      existing_condas <- reticulate::conda_list()$name
      env_name %in% existing_condas
    }, error = function(e) FALSE)
    
    if (virtualenv_exists) {
      reticulate::use_virtualenv(env_name, required = TRUE)
      options(CASSIA.env_name = env_name)
      options(CASSIA.env_method = "virtualenv")
      message("Using virtualenv: ", env_name)
    } else if (conda_exists) {
      reticulate::use_condaenv(env_name, required = TRUE)
      options(CASSIA.env_name = env_name)
      options(CASSIA.env_method = "conda")
      # Maintain backward compatibility
      options(CASSIA.conda_env = env_name)
      message("Using conda environment: ", env_name)
    } else {
      stop("Environment '", env_name, "' not found in either virtualenv or conda")
    }
    
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
#' @param cassia_cluster_col Name of the column in CASSIA results containing cluster IDs (default: "Cluster ID")
#' @param prefix Prefix to add to the new metadata columns (default: "CASSIA_")
#' @param replace_existing Whether to replace existing annotations (default: FALSE)
#' @param fuzzy_match Whether to perform fuzzy matching on cluster names (default: TRUE)
#' @param columns_to_include Level of columns to include: 1 for only merged groupings, 2 for all metrics (default: 2)
#'
#' @return A Seurat object with CASSIA annotations added as metadata columns
#' @importFrom utils read.csv
#' @importFrom Seurat AddMetaData
#' @export
add_cassia_to_seurat <- function(seurat_obj, cassia_results_path, cluster_col = "seurat_clusters",
                                cassia_cluster_col = "Cluster ID", prefix = "CASSIA_", 
                                replace_existing = FALSE, fuzzy_match = TRUE,
                                columns_to_include = 1) {
  # Check if required packages are installed
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required. Please install it with 'install.packages(\"Seurat\")'")
  }
  
  # Check if file exists
  if (!file.exists(cassia_results_path)) {
    stop("CASSIA results file not found at: ", cassia_results_path)
  }
  
  # Validate columns_to_include parameter
  if (!columns_to_include %in% c(1, 2)) {
    warning("Invalid value for columns_to_include. Using default value 2 (all metrics).")
    columns_to_include <- 2
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
    alt_cluster_cols <- c("Cluster ID", "True Cell Type", "Cluster", "cluster_id", "ClusterID", "Cluster_ID", "cluster")
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
    
    # Using a functional approach to match clusters
    # First create data frames for both normalized sets
    seurat_df <- data.frame(
      id = seurat_norm,
      seurat = seurat_clusters_unique,
      stringsAsFactors = FALSE
    )
    
    cassia_df <- data.frame(
      id = cassia_norm,
      cassia = cassia_clusters_unique,
      stringsAsFactors = FALSE
    )
    
    # First, do exact matches
    exact_matches <- merge(seurat_df, cassia_df, by = "id", all = TRUE)
    
    # Define a function to find the first matching pattern
    find_first_match <- function(pattern, target_vector) {
      matches <- vapply(target_vector, function(x) 
        grepl(pattern, x) || grepl(x, pattern), logical(1))
      
      if (any(matches)) {
        return(target_vector[which(matches)[1]])
      } else {
        return(NA_character_)
      }
    }
    
    # For unmatched Seurat clusters, find matching CASSIA clusters
    missing_seurat_idx <- which(is.na(exact_matches$cassia))
    if (length(missing_seurat_idx) > 0) {
      exact_matches$cassia[missing_seurat_idx] <- sapply(
        exact_matches$id[missing_seurat_idx],
        find_first_match,
        target_vector = cassia_norm
      )
      
      # Replace the index values with actual CASSIA cluster names
      for (i in missing_seurat_idx) {
        if (!is.na(exact_matches$cassia[i])) {
          match_idx <- which(cassia_norm == exact_matches$cassia[i])[1]
          exact_matches$cassia[i] <- cassia_clusters_unique[match_idx]
        }
      }
    }
    
    # For unmatched CASSIA clusters, find matching Seurat clusters
    # (Less important but included for completeness)
    missing_cassia_idx <- which(is.na(exact_matches$seurat))
    if (length(missing_cassia_idx) > 0) {
      exact_matches$seurat[missing_cassia_idx] <- sapply(
        exact_matches$id[missing_cassia_idx],
        find_first_match,
        target_vector = seurat_norm
      )
      
      # Replace the index values with actual Seurat cluster names
      for (i in missing_cassia_idx) {
        if (!is.na(exact_matches$seurat[i])) {
          match_idx <- which(seurat_norm == exact_matches$seurat[i])[1]
          exact_matches$seurat[i] <- seurat_clusters_unique[match_idx]
        }
      }
    }
    
    # Create final cluster mapping
    cluster_map <- exact_matches[, c("seurat", "cassia")]
    
    # Print the mapping if any clusters were mapped
    mapped_clusters <- cluster_map[!is.na(cluster_map$cassia) & 
                                    !is.na(cluster_map$seurat) & 
                                    cluster_map$seurat != cluster_map$cassia, ]
    
    if (nrow(mapped_clusters) > 0) {
      message("Created mapping between Seurat clusters and CASSIA clusters:")
      mapply(function(s, c) message("  '", s, "' -> '", c, "'"), 
             mapped_clusters$seurat, mapped_clusters$cassia)
    }
    
    # Check for any unmatched clusters
    unmatched <- cluster_map$seurat[is.na(cluster_map$cassia)]
    if (length(unmatched) > 0) {
      warning("Could not find matches for these Seurat clusters: ", 
              paste(unmatched, collapse = ", "))
    }
  }
  
  # Define the column mappings for CASSIA output (new names first, then old names for backward compatibility)
  column_mapping <- list(
    general_celltype = c("Predicted General Cell Type", "Predicted Main Cell Type", "General Cell Type", "General_Cell_Type"),
    sub_celltype = c("Predicted Detailed Cell Type", "Predicted Sub Cell Types", "Sub Cell Type", "Sub_Cell_Type"),
    mixed_celltype = c("Possible Mixed Cell Types", "Mixed Cell Types", "Possible_Mixed_Cell_Types"),
    score = c("Score", "Consensus Score", "Consensus_Score")
  )
  
  # Add merged grouping columns to the mapping
  merged_grouping_mapping <- list(
    merged_grouping_1 = c("Merged Grouping 1", "Merged_Grouping_1"),
    merged_grouping_2 = c("Merged Grouping 2", "Merged_Grouping_2"),
    merged_grouping_3 = c("Merged Grouping 3", "Merged_Grouping_3")
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
  
  # Check which merged grouping columns exist in the results
  available_merged_groups <- list()
  for (group_name in names(merged_grouping_mapping)) {
    col_candidates <- merged_grouping_mapping[[group_name]]
    found_col <- col_candidates[col_candidates %in% colnames(cassia_results)]
    
    if (length(found_col) > 0) {
      col_name <- found_col[1]
      available_merged_groups[[group_name]] <- setNames(
        cassia_results[[col_name]],
        cassia_results[[cassia_cluster_col]]
      )
    }
  }
  
  # If columns_to_include is 1, only include merged groupings
  if (columns_to_include == 1) {
    if (length(available_merged_groups) == 0) {
      warning("No merged grouping columns found in CASSIA results. Will include standard metrics instead.")
      columns_to_include <- 2
    } else {
      lookup_tables <- available_merged_groups
    }
  }
  
  # If columns_to_include is 2 or no merged groupings were found, include all metrics
  if (columns_to_include == 2) {
    # Process standard metrics
    std_lookup_tables <- lapply(names(column_mapping), function(anno_type) {
      col_candidates <- column_mapping[[anno_type]]
      found_col <- col_candidates[col_candidates %in% colnames(cassia_results)]
      
      if (length(found_col) > 0) {
        col_name <- found_col[1]
        
        # For annotation types other than sub_celltype, return a simple named vector
        if (anno_type != "sub_celltype") {
          return(setNames(
            cassia_results[[col_name]],
            cassia_results[[cassia_cluster_col]]
          ))
        }
        
        # For sub_celltype, return NULL (handled separately)
        return(NULL)
      }
      return(NULL)
    })
    names(std_lookup_tables) <- names(column_mapping)
    
    # Filter out NULL entries
    std_lookup_tables <- std_lookup_tables[!sapply(std_lookup_tables, is.null)]
    
    # Combine with merged groupings if available
    lookup_tables <- c(std_lookup_tables, available_merged_groups)
  }
  
  # Special handling for sub cell types (only for columns_to_include == 2)
  if (columns_to_include == 2 && 
      any(sapply(column_mapping$sub_celltype, function(x) x %in% colnames(cassia_results)))) {
    col_name <- column_mapping$sub_celltype[column_mapping$sub_celltype %in% colnames(cassia_results)][1]
    
    # Store the original values
    sub_types_full <- cassia_results[[col_name]]
    lookup_tables_split$sub_celltype_all <- setNames(
      sub_types_full,
      cassia_results[[cassia_cluster_col]]
    )
    
    # Process all subtypes with a single lapply call
    split_results <- lapply(sub_types_full, function(x) {
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
    
    # Extract each position into a separate lookup table
    lookup_tables_split$sub_celltype_1 <- setNames(
      sapply(split_results, `[`, 1),
      cassia_results[[cassia_cluster_col]]
    )
    
    lookup_tables_split$sub_celltype_2 <- setNames(
      sapply(split_results, `[`, 2),
      cassia_results[[cassia_cluster_col]]
    )
    
    lookup_tables_split$sub_celltype_3 <- setNames(
      sapply(split_results, `[`, 3),
      cassia_results[[cassia_cluster_col]]
    )
    
    # Set the main sub_celltype to the most likely one
    lookup_tables$sub_celltype <- lookup_tables_split$sub_celltype_1
  }
  
  # Filter out NULL entries from lookup_tables
  lookup_tables <- lookup_tables[!sapply(lookup_tables, is.null)]
  
  # Get cluster IDs from Seurat object
  seurat_clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  
  # Apply cluster mapping if it exists
  seurat_clusters_mapped <- if (!is.null(cluster_map) && any(!is.na(cluster_map$cassia))) {
    # Use vectorized match operation instead of sapply with function
    matched_indices <- match(seurat_clusters, cluster_map$seurat)
    matched_cassia <- cluster_map$cassia[matched_indices]
    
    # Where there's no match, use the original value
    ifelse(is.na(matched_cassia), seurat_clusters, matched_cassia)
  } else {
    seurat_clusters
  }
  
  # Create a new data frame for metadata
  new_metadata <- data.frame(row.names = colnames(seurat_obj))
  
  # For each annotation type, map cluster IDs to annotations (vectorized)
  metadata_list <- lapply(names(lookup_tables), function(anno_type) {
    lookup <- lookup_tables[[anno_type]]
    annotations <- lookup[seurat_clusters_mapped]

    meta_col_name <- paste0(prefix, anno_type)
    result <- data.frame(col = annotations, stringsAsFactors = FALSE)
    names(result) <- meta_col_name

    # Handle replacement of existing annotations
    if (replace_existing && anno_type == "general_celltype") {
      if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
        message("Replacing 'cell_type' column with CASSIA annotations")
        seurat_obj$cell_type <- annotations
      }
    }

    return(result)
  })
  new_metadata <- do.call(cbind, metadata_list)
  rownames(new_metadata) <- colnames(seurat_obj)
  
  # For sub cell types, add all the split versions (only for columns_to_include == 2)
  if (columns_to_include == 2) {
    sub_metadata <- lapply(names(lookup_tables_split), function(split_type) {
      lookup <- lookup_tables_split[[split_type]]
      if (is.null(lookup)) return(NULL)
      
      annotations <- lookup[seurat_clusters_mapped]
      meta_col_name <- paste0(prefix, split_type)
      result <- data.frame(col = annotations, stringsAsFactors = FALSE)
      names(result) <- meta_col_name
      return(result)
    })
    sub_metadata <- Filter(Negate(is.null), sub_metadata)
    
    if (length(sub_metadata) > 0) {
      sub_metadata <- do.call(cbind, sub_metadata)
      rownames(sub_metadata) <- colnames(seurat_obj)
      new_metadata <- cbind(new_metadata, sub_metadata)
    }
    
    # Add combined general celltype + first subcelltype column
    if (paste0(prefix, "general_celltype") %in% names(new_metadata) && 
        paste0(prefix, "sub_celltype_1") %in% names(new_metadata)) {
      
      general_col <- paste0(prefix, "general_celltype")
      sub_col <- paste0(prefix, "sub_celltype_1")
      combined_col <- paste0(prefix, "combined_celltype")
      
      # Create combined column with vectorized operation - using a more visible separator
      new_metadata[[combined_col]] <- paste(
        new_metadata[[general_col]],
        new_metadata[[sub_col]],
        sep = " :: "
      )
      
      # Handle cases where subcelltype is NA
      na_idx <- which(is.na(new_metadata[[sub_col]]))
      if (length(na_idx) > 0) {
        new_metadata[[combined_col]][na_idx] <- new_metadata[[general_col]][na_idx]
      }
      
      message("Added combined cell type column: ", combined_col)
    }
  }
  
  # Add all metadata columns at once
  seurat_obj <- Seurat::AddMetaData(seurat_obj, new_metadata)
  
  message("Added CASSIA annotations to Seurat object")
  message("New metadata columns: ", paste(names(new_metadata), collapse = ", "))
  
  return(seurat_obj)
}


