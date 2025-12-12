---
title: Subclustering Analysis (Optional)
---

## Overview

Subclustering analysis is a powerful technique for studying specific cell populations in greater detail. This feature allows you to analyze subclustered populations, such as T cells or fibroblasts, after initial CASSIA annotation.

## Quick Start

```python
CASSIA.runCASSIA_subclusters(
    marker = subcluster_results,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

## Input

### Workflow Summary
1. Initial CASSIA analysis on full dataset
2. Subcluster extraction and processing (using Seurat or Scanpy)
3. Marker identification for subclusters
4. CASSIA subclustering analysis
5. Uncertainty assessment (optional)

### Required Input
- **Marker genes**: DataFrame or file path containing marker genes for each subcluster (output from `FindAllMarkers` in Seurat or `sc.tl.rank_genes_groups` in Scanpy)

We recommend applying the default CASSIA first. Then, on a target cluster, apply standard pipelines (Seurat/Scanpy) to subcluster and get marker results.

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `marker` | Marker genes for the subclusters (DataFrame or file path) |
| `major_cluster_info` | Description of the parent cluster or context (e.g., "CD8+ T cells" or "cd8 t cell mixed with other celltypes") |
| `output_name` | Base name for the output CSV file |
| `model` | LLM model to use |
| `provider` | API provider |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `temperature` | 0 | Sampling temperature (0-1) |
| `n_genes` | 50 | Number of top marker genes to use |

### Uncertainty Assessment Functions

For more confident results, calculate consistency scores (CS) using multiple iterations:

**`runCASSIA_n_subcluster()`** - Run multiple annotation iterations:

```python
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
```

| Parameter | Description |
|-----------|-------------|
| `n` | Number of iterations to run |
| `base_output_name` | Base name for output files (appended with iteration number) |
| `max_workers` | Number of parallel workers |

**`runCASSIA_similarity_score_batch()`** - Calculate similarity scores across iterations:

```python
CASSIA.runCASSIA_similarity_score_batch(
    marker = subcluster_results,
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
| `{output_name}.csv` | Basic CASSIA analysis results |
| `{output_name}.html` | HTML report with visualizations |
| `{output_name}_uncertainty.csv` | Similarity scores (if uncertainty assessment is performed) |

> **Automatic Report Generation**: An HTML report is automatically generated alongside the CSV output for easy visualization of subclustering results.
