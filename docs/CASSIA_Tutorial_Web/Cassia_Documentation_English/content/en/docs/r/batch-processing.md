---
title: Batch Processing
---

CASSIA supports batch processing to analyze multiple clusters simultaneously. This guide explains how to prepare your data and run batch analysis efficiently.


### Model Recommendations
For detailed information on model settings and recommendations, please refer to the **[How to Select Models and Providers](setting-up-cassia.md#how-to-select-models-and-providers)** section.

### Preparing Marker Data
You have three options for providing marker data:

1. Create a data frame or CSV file containing your clusters and marker genes
2. Use Seurat's `findAllMarkers` function output directly
3. Use CASSIA's example marker data

```R
# Option 1: Load your own marker data
markers <- read.csv("path/to/your/markers.csv")

# Option 2: Use Seurat's findAllMarkers output directly
# (assuming you already have a Seurat object)
markers <- FindAllMarkers(seurat_obj)

# Option 3: Load example marker data
markers <- loadExampleMarkers()

# Preview the data
head(markers)
```

#### Marker Data Format
CASSIA accepts two formats:

1. **FindAllMarkers Output**: The standard output from Seurat's FindAllMarkers function
2. **Simplified Format**: A two-column data frame where:
   - First column: cluster identifier
   - Second column: comma-separated ranked marker genes

### Running Batch Analysis

#### Setting Up Parameters

```R
# Detect available CPU cores
available_cores <- parallel::detectCores()

# Calculate recommended workers (75% of available cores)
recommended_workers <- floor(available_cores * 0.75)

runCASSIA_batch(
    # Required parameters
    marker = markers,                    # Marker data (data frame or file path)
    output_name = "my_annotation",       # Base name for output files
    model = "anthropic/claude-4.5-sonnet", # Model to use
    tissue = "brain",                    # Tissue type
    species = "human",                   # Species
    
    # Optional parameters
    max_workers = recommended_workers,    # Number of parallel workers
    n_genes = 50,                        # Number of top marker genes to use
    additional_info = "",                # Additional context
    provider = "openrouter",              # API provider
    
    # Advanced parameters
    ranking_method = "avg_log2FC",       # Ranking method: "avg_log2FC", "p_val_adj", "pct_diff", or "Score"
    ascending = NULL,                    # Sort direction (NULL uses method default)
    temperature = 0,                     # Model temperature
    celltype_column = NULL,              # Column name for clusters (default: first column)
    gene_column_name = NULL,             # Column name for genes (default: second column)
    max_retries = 1,                     # Max retries for failed calls
    validator_involvement = "v1"         # Validator level: "v1" (moderate) or "v0" (high)
)
```

### Parameter Details

1. **Marker Gene Selection**:
   - Default: top 50 genes per cluster
   - `ranking_method`: Controls how marker genes are ranked and selected
     - `"avg_log2FC"` (default): Rank by average log2 fold change
     - `"p_val_adj"`: Rank by adjusted p-value
     - `"pct_diff"`: Rank by difference in percentage expression
     - `"Score"`: Rank by custom score
   - Filtering criteria:
     - Adjusted p-value < 0.05
     - Average log2 fold change > 0.25
     - Minimum percentage > 0.1
   - If fewer than 50 genes pass filters, all passing genes are used

2. **Parallel Processing**:
   - `max_workers`: Controls parallel processing threads
   - Recommended: 80% of available CPU cores
   - Example: For a 16-core machine, set to 13

3. **Additional Context** (optional):
   - Use `additional_info` to provide experimental context
   - Examples:
     - Treatment conditions: "Samples were antibody-treated"
     - Analysis focus: "Please carefully distinguish between cancer and non-cancer cells"
   - Tip: Compare results with and without additional context
   
4. **Model Selection**:
   - Default is `anthropic/claude-4.5-sonnet` for best performance.
   - You can use `google/gemini-2.5-flash` for a faster, preliminary look at your data.

### Output Files

The analysis generates three files:
1. `my_annotation_full.csv`: Complete conversation history
2. `my_annotation_summary.csv`: Condensed results summary
3. `my_annotation_report.html`: An interactive HTML report visualizing the results

### Add CASSIA Results Back to Seurat Object

You can easily add the annotation results back to your Seurat object using `add_cassia_to_seurat`. This function maps the CASSIA results to your Seurat object based on cluster identifiers.

```R
seurat_corrected <- add_cassia_to_seurat(
  seurat_obj = seurat_corrected,                         # The seurat object you want to add the CASSIA results to
  cassia_results_path = "my_annotation_summary.csv",     # Path to the summary results file
  cluster_col = "seurat_clusters",                       # Column in Seurat object with cluster IDs
  cassia_cluster_col = "Cluster"                         # Column in the CASSIA results with cluster IDs
)
```

This will add multiple new columns to your Seurat object, including:
- `General Cell Type`: The broad cell type category
- `True Cell Type`: The most likely specific cell type
- `Alternative Cell Type 1 & 2`: Other possible cell types
- `Mixed Cell Type`: Information on potential mixed populations
- `Quality Score`: If available, the confidence score for the annotation

### Tips for Optimal Results

1. **Resource Management**:
   - Monitor system resources when setting `max_workers`
   - Start with recommended 75% of cores and adjust if needed

2. **Marker Gene Selection**:
   - Default 50 genes works well for most cases
   - Increase for more complex cell types
   - Decrease if running into API rate limits

3. **Context Optimization**:
   - Test runs with and without additional context
   - Keep context concise and relevant
