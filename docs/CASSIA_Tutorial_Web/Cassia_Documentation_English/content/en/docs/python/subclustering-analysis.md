---
title: Subclustering Analysis (Optional)
---

Subclustering analysis is a powerful technique for studying specific cell populations in greater detail. This tutorial walks you through the process of analyzing subclustered populations, such as T cells or fibroblasts.

### Workflow Summary
1. Initial CASSIA analysis
2. Subcluster extraction and processing (using Seurat or Scanpy)
3. Marker identification
4. CASSIA subclustering analysis
5. Uncertainty assessment (optional)

### Running Subclustering Analysis

We recommend applying the default CASSIA first. Then, on a target cluster, apply standard pipelines (Seurat/Scanpy) to subcluster and get marker results.

```python
CASSIA.runCASSIA_subclusters(
    marker = subcluster_results,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

#### Parameter Details

- **`marker`**: Marker genes for the subclusters (data frame or file path).
- **`major_cluster_info`**: Description of the parent cluster or context (e.g., "CD8+ T cells" or "cd8 t cell mixed with other celltypes").
- **`output_name`**: Base name for the output CSV file.
- **`model`**: LLM model to use.
- **`provider`**: API provider.
- **`temperature`**: Sampling temperature (0-1).
- **`n_genes`**: Number of top marker genes to use.

### Uncertainty Assessment

For more confident results, calculate consistency scores (CS):

```python
# Run multiple iterations
CASSIA.runCASSIA_n_subcluster(
    n=5, 
    marker=subcluster_results,
    major_cluster_info="cd8 t cell",
    base_output_name="subclustering_results_n",
    model="anthropic/claude-sonnet-4.5",
    temperature=0,
    provider="openrouter",
    max_workers=5,
    n_genes=50
)

# Calculate similarity scores
CASSIA.runCASSIA_similarity_score_batch(
    marker = subcluster_results,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

### Output Files
- `{output_name}.csv`: Basic Cassia analysis results.
- `{output_name}_uncertainty.csv`: Similarity scores (if uncertainty assessment is performed).
