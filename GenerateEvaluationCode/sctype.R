# Define a function to process single-cell data
process_data_sc5k <- function(data_sc5k, known_tissue_type = "Immune system", 
                              custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", 
                              resolution = 0.8, dims = 1:10) {
  
  # Normalize the data using log-normalization
  data_sc5k <- NormalizeData(data_sc5k, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identify highly variable features
  data_sc5k <- FindVariableFeatures(data_sc5k, selection.method = "vst", nfeatures = 2000)
  
  # Scale the data and perform principal component analysis (PCA)
  data_sc5k <- ScaleData(data_sc5k)
  data_sc5k <- RunPCA(data_sc5k)
  
  # Find neighbors and identify clusters
  data_sc5k <- FindNeighbors(data_sc5k, dims = dims)
  data_sc5k <- FindClusters(data_sc5k, resolution = resolution)
  
  # Perform UMAP dimensionality reduction
  data_sc5k <- RunUMAP(data_sc5k, dims = dims)
  
  # Generate a UMAP plot
  DimPlot(data_sc5k, reduction = "umap")
  
  # Run ScType for cell type classification
  data_sc5k <- run_sctype(data_sc5k, known_tissue_type = known_tissue_type, 
                          custom_marker_file = custom_marker_file, 
                          name = "sctype_classification", plot = TRUE)
  
  # Display the ScType classification results
  sctype_table <- data_sc5k$sctype_classification %>% table
  print(sctype_table)
  
  # Return the processed Seurat object
  return(data_sc5k)
}


# Example usage of the process_data_sc5k function on multiple datasets
data_sc5k <- process_data_sc5k(data_sc5k)
