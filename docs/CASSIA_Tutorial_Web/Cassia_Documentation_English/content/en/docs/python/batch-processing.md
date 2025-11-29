---
title: Batch Processing
---

CASSIA supports batch processing to analyze multiple clusters simultaneously. This guide explains how to run batch analysis efficiently.

### Preparing Marker Data

For Python, marker data is typically loaded from CSV or processed dataframes.

```python
# Load your marker data
# Can be Seurat's FindAllMarkers output format
markers = CASSIA.loadmarker(marker_type="unprocessed")
```

### Running Batch Analysis

```python
output_name="intestine_detailed"

# Run batch analysis
CASSIA.runCASSIA_batch(
    marker = markers,
    output_name = output_name,
    model = "anthropic/claude-sonnet-4.5",
    tissue = "large intestine",
    species = "human",
    max_workers = 6,  # Matching cluster count
    n_genes = 50,
    additional_info = None,
    provider = "openrouter"
)
```

### Parameter Details

- **`marker`**: Marker data (DataFrame or file path). The data should contain cluster IDs and gene names.
- **`output_name`**: Base name for output files.
- **`model`**: Model to use (e.g., "anthropic/claude-sonnet-4.5").
- **`tissue`**: Tissue type (e.g., "brain").
- **`species`**: Species (e.g., "human").
- **`max_workers`**: Number of parallel workers. Recommended to use ~75-80% of available CPU cores.
- **`n_genes`**: Number of top marker genes to use per cluster (default 50).
- **`additional_info`**: Additional context about the experiment (optional).
- **`provider`**: API provider ("openrouter", "openai", "anthropic").
- **`ranking_method`**: Method to rank genes ("avg_log2FC", "p_val_adj", "pct_diff", "Score").
- **`ascending`**: Sort direction (bool).
- **`temperature`**: Model temperature (0-1).
- **`celltype_column`**: Column name for cluster IDs.
- **`gene_column_name`**: Column name for gene symbols.
- **`max_retries`**: Max retries for failed API calls.
- **`validator_involvement`**: Validation level ("v1" or "v0").

### Output Files

The analysis generates the following files:
1. `{output_name}_full.csv`: Complete conversation history and reasoning.
2. `{output_name}_summary.csv`: Condensed results summary.
3. `{output_name}_report.html`: An interactive HTML report visualizing the results (generated separately).

### Tips for Optimal Results

1. **Resource Management**: Monitor system resources when setting `max_workers`. Start with a conservative number if unsure.
2. **Context Optimization**: Use `additional_info` to provide key experimental details (e.g., "tumor sample", "treated condition") to guide the model.
