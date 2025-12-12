---
title: Subclustering Analysis (Optional)
---

## Overview

Subclustering analysis is a powerful technique for studying specific cell populations in greater detail. This feature allows you to analyze subclustered populations, such as T cells or fibroblasts, using Cassia and Seurat.

## Quick Start

```r
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

## Input

### Prerequisites
- A Seurat object containing your single-cell data
- The Cassia package installed and loaded
- Basic familiarity with R and single-cell analysis

### Workflow Summary
1. Initial Cassia analysis on full dataset
2. Subcluster extraction and processing
3. Marker identification for subclusters
4. Cassia subclustering analysis
5. Uncertainty assessment (optional)

### Preparing Subclusters

First, run the default Cassia pipeline on your complete dataset to identify major cell populations. Then extract and process your target cluster using Seurat:

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

### Marker Identification

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

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `marker` | Marker genes for the subclusters (data frame or file path) |
| `major_cluster_info` | Description of the parent cluster or context (e.g., "CD8+ T cells" or "cd8 t cell mixed with other celltypes") |
| `output_name` | Base name for the output CSV file |
| `model` | LLM model to use |
| `provider` | API provider |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `temperature` | 0 | Sampling temperature (0-1) |
| `n_genes` | 50 | Number of top marker genes to use |

### Example: Mixed Populations

```r
# For mixed populations
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell mixed with other celltypes",
    output_name = "subclustering_results2",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

### Uncertainty Assessment Functions

For more confident results, calculate CS scores using multiple iterations:

**`runCASSIA_n_subcluster()`** - Run multiple annotation iterations:

```r
runCASSIA_n_subcluster(
    n = 5,
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    base_output_name = "subclustering_results_n",
    model = "anthropic/claude-sonnet-4.5",
    temperature = 0,
    provider = "openrouter",
    max_workers = 5,
    n_genes = 50
)
```

| Parameter | Description |
|-----------|-------------|
| `n` | Number of iterations to run |
| `base_output_name` | Base name for output files (appended with iteration number) |
| `max_workers` | Number of parallel workers |

**`runCASSIA_similarity_score_batch()`** - Calculate similarity scores across iterations:

```r
similarity_scores <- runCASSIA_similarity_score_batch(
    marker = marker_sub,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `file_pattern` | - | Glob pattern matching iteration result files |
| `main_weight` | 0.5 | Weight for main cell type similarity |
| `sub_weight` | 0.5 | Weight for subtype similarity |

## Output

| File | Description |
|------|-------------|
| `cd8_subcluster_markers.csv` | Marker genes for each subcluster |
| `{output_name}.csv` | Basic Cassia analysis results |
| `{output_name}.html` | HTML report with visualizations |
| `{output_name}_uncertainty.csv` | Similarity scores (if uncertainty assessment is performed) |

> **Automatic Report Generation**: An HTML report is automatically generated alongside the CSV output for easy visualization of subclustering results.

### Tips and Recommendations
- Always run the default Cassia analysis first before subclustering
- Adjust clustering resolution based on your data's complexity
- When dealing with mixed populations, specify this in the `major_cluster_info` parameter
- Use the uncertainty assessment for more robust results
