---
title: Subclustering Analysis (Optional)
---

Subclustering analysis is a powerful technique for studying specific cell populations in greater detail. This tutorial walks you through the process of analyzing subclustered populations, such as T cells or fibroblasts, using Cassia and Seurat.

## Prerequisites
- A Seurat object containing your single-cell data
- The Cassia package installed and loaded
- Basic familiarity with R and single-cell analysis

## Workflow Summary
1. Initial Cassia analysis
2. Subcluster extraction and processing
3. Marker identification
4. Cassia subclustering analysis
5. Uncertainty assessment (optional)

## Detailed Steps

### 1. Initial Analysis
First, run the default Cassia pipeline on your complete dataset to identify major cell populations.

### 2. Subcluster Processing
Extract and process your target cluster using Seurat:

```r
# Extract target population (example using CD8+ T cells)
cd8_cells <- subset(large, cell_ontology_class == "cd8-positive, alpha-beta t cell")

# Normalize data
cd8_cells <- NormalizeData(cd8_cells)

# Identify variable features
cd8_cells <- FindVariableFeatures(cd8_cells, 
    selection.method = "vst", 
    nfeatures = 2000)

# Scale data
all.genes <- rownames(cd8_cells)
cd8_cells <- ScaleData(cd8_cells, features = all.genes)

# Run PCA
cd8_cells <- RunPCA(cd8_cells, 
    features = VariableFeatures(object = cd8_cells),
    npcs = 30)

# Perform clustering
cd8_cells <- FindNeighbors(cd8_cells, dims = 1:20)
cd8_cells <- FindClusters(cd8_cells, resolution = 0.3)

# Generate UMAP visualization
cd8_cells <- RunUMAP(cd8_cells, dims = 1:20)
```

### 3. Marker Identification
Identify markers for each subcluster:

```r
# Find markers
cd8_markers <- FindAllMarkers(cd8_cells,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25)

# Filter significant markers
cd8_markers <- cd8_markers %>% filter(p_val_adj < 0.05)

# Save results
write.csv(cd8_markers, "cd8_subcluster_markers.csv")
```

### 4. Cassia Subclustering Analysis
Run Cassia analysis on the subclusters:

```r
# Basic analysis
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-4.5-sonnet",
    provider = "openrouter",
    temperature = 0,
    n_genes = 50
)
```

**Parameter Details:**
- **`marker`**: Marker genes for the subclusters (data frame or file path).
- **`major_cluster_info`**: Description of the parent cluster or context (e.g., "CD8+ T cells").
- **`output_name`**: Base name for the output CSV file.
- **`model`**: LLM model to use.
- **`provider`**: API provider.
- **`temperature`**: Sampling temperature (0-1, default 0).
- **`n_genes`**: Number of top marker genes to use (default 50).

```r
# For mixed populations
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell mixed with other celltypes",
    output_name = "subclustering_results2",
    model = "anthropic/claude-4.5-sonnet",
    provider = "openrouter"
)
```

### 5. Uncertainty Assessment (Optional)
For more confident results, calculate CS scores:

```r
# Run multiple iterations
runCASSIA_n_subcluster(
    n = 5,
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    base_output_name = "subclustering_results_n", # Note: parameter name is base_output_name
    model = "anthropic/claude-4.5-sonnet",
    temperature = 0,
    provider = "openrouter",
    max_workers = 5,
    n_genes = 50
)

# Calculate similarity scores
similarity_scores <- runCASSIA_similarity_score_batch(
    marker = marker_sub,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "anthropic/claude-4.5-sonnet",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

## Tips and Recommendations
- Always run the default Cassia analysis first before subclustering
- Adjust clustering resolution based on your data's complexity
- When dealing with mixed populations, specify this in the `major_cluster_info` parameter
- Use the uncertainty assessment for more robust results

## Output Files
The analysis generates several output files:
- `cd8_subcluster_markers.csv`: Marker genes for each subcluster
- `subclustering_results.csv`: Basic Cassia analysis results
- `subclustering_uncertainty.csv`: Similarity scores (if uncertainty assessment is performed)
