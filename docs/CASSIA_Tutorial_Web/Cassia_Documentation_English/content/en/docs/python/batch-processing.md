---
title: Batch Processing
---

## Overview

Batch processing allows you to annotate all clusters in your single-cell dataset in one run. CASSIA processes multiple clusters in parallel, generating cell type predictions with detailed reasoning for each cluster.

---

## Quick Start

```python
CASSIA.runCASSIA_batch(
    marker = markers,
    output_name = "my_annotation",
    model = "anthropic/claude-sonnet-4.5",
    tissue = "brain",
    species = "human",
    provider = "openrouter"
)
```

---

## Input

### Marker Data Formats

CASSIA accepts three input formats:

**1. Seurat FindAllMarkers Output**

Standard output from Seurat's `FindAllMarkers` function with differential expression statistics:

| p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | cluster | gene |
|-------|------------|-------|-------|-----------|---------|------|
| 0 | 3.02 | 0.973 | 0.152 | 0 | 0 | CD79A |
| 0 | 2.74 | 0.938 | 0.125 | 0 | 0 | MS4A1 |
| 0 | 2.54 | 0.935 | 0.138 | 0 | 0 | CD79B |
| 0 | 1.89 | 0.812 | 0.089 | 0 | 1 | IL7R |
| 0 | 1.76 | 0.756 | 0.112 | 0 | 1 | CCR7 |

**2. Scanpy rank_genes_groups Output (Recommended for Python)**

Output from Scanpy's `sc.tl.rank_genes_groups()` function, typically exported using `sc.get.rank_genes_groups_df()`:

| group | names | scores | pvals | pvals_adj | logfoldchanges |
|-------|-------|--------|-------|-----------|----------------|
| 0 | CD79A | 28.53 | 0 | 0 | 3.02 |
| 0 | MS4A1 | 25.41 | 0 | 0 | 2.74 |
| 0 | CD79B | 24.89 | 0 | 0 | 2.54 |
| 1 | IL7R | 22.15 | 0 | 0 | 1.89 |
| 1 | CCR7 | 20.87 | 0 | 0 | 1.76 |

**3. Simplified Format**

A two-column DataFrame with cluster ID and comma-separated marker genes:

| cluster | marker_genes |
|---------|--------------|
| 0 | CD79A,MS4A1,CD79B,HLA-DRA,TCL1A |
| 1 | IL7R,CCR7,LEF1,TCF7,FHIT,MAL |
| 2 | CD8A,CD8B,GZMK,CCL5,NKG7 |

### Loading Marker Data

```python
import CASSIA
import scanpy as sc
import pandas as pd

# Option 1: Load Seurat FindAllMarkers output (exported as CSV)
markers = pd.read_csv("seurat_markers.csv")

# Option 2: Use Scanpy rank_genes_groups output directly (Recommended for Python)
# (assuming you already have an AnnData object with rank_genes_groups computed)
markers = sc.get.rank_genes_groups_df(adata, group=None)  # Get all groups

# Option 3: Load your own simplified marker data
markers = pd.read_csv("path/to/your/markers.csv")

# Load example marker data for testing
markers = CASSIA.loadmarker(marker_type="unprocessed")

# Preview the data
print(markers.head())
```

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `marker` | Marker data (DataFrame or file path) |
| `output_name` | Base name for output files |
| `model` | LLM model ID (see Model Selection below) |
| `tissue` | Tissue type (e.g., `"brain"`, `"blood"`) |
| `species` | Species (e.g., `"human"`, `"mouse"`) |
| `provider` | API provider (`"openrouter"`, `"openai"`, `"anthropic"`) or custom base URL |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_workers` | 4 | Number of parallel workers. Recommend ~75% of CPU cores. |
| `n_genes` | 50 | Top marker genes per cluster |
| `additional_info` | `None` | Extra experimental context (see below) |
| `temperature` | 0 | Output randomness (0=deterministic, 1=creative). Keep at 0 for reproducible results. |
| `validator_involvement` | `"v1"` | Validation intensity: `"v1"` (moderate) or `"v0"` (high, slower) |
| `reasoning` | `None` | Reasoning depth for GPT-5 series via OpenRouter only: `"low"`, `"medium"`, `"high"` |
| `ranking_method` | `"avg_log2FC"` | Gene ranking: `"avg_log2FC"`, `"p_val_adj"`, `"pct_diff"`, `"Score"` |
| `ascending` | `None` | Sort direction (uses method default) |
| `celltype_column` | `None` | Column name for cluster IDs |
| `gene_column_name` | `None` | Column name for gene symbols |
| `max_retries` | 1 | Max retries for failed API calls |

### Parameter Details

**Model Selection**
- Default is `anthropic/claude-sonnet-4.5` for best performance
- Use `google/gemini-2.5-flash` for faster, preliminary analysis
- For detailed model recommendations, see [How to Select Models and Providers](setting-up-cassia.md#how-to-select-models-and-providers)

**Marker Gene Selection**
- Default: top 50 genes per cluster
- `ranking_method` controls how marker genes are ranked and selected:
  - `"avg_log2FC"` (default): Rank by average log2 fold change
  - `"p_val_adj"`: Rank by adjusted p-value
  - `"pct_diff"`: Rank by difference in percentage expression
  - `"Score"`: Rank by custom score
- Filtering criteria: adjusted p-value < 0.05, avg_log2FC > 0.25, min percentage > 0.1
- If fewer than 50 genes pass filters, all passing genes are used

**Additional Context**
- Use `additional_info` to provide experimental context
- Examples:
  - Treatment conditions: `"Samples were antibody-treated"`
  - Analysis focus: `"Please carefully distinguish between cancer and non-cancer cells"`
- Tip: Compare results with and without additional context

---

## Output

### Files Generated

| File | Description |
|------|-------------|
| `{output_name}_summary.csv` | Annotation results with cell types, markers, and metadata |
| `{output_name}_conversations.json` | Complete conversation history for debugging |
| `{output_name}_report.html` | Interactive HTML report with visualizations |
