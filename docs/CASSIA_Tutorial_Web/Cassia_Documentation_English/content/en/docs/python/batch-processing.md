---
title: Batch Processing
---

CASSIA supports batch processing to analyze multiple clusters simultaneously. This guide explains how to run batch analysis efficiently.

### Preparing Marker Data

You have three options for providing marker data:

1. Use Seurat's `FindAllMarkers` output (exported as CSV)
2. Use Scanpy's `rank_genes_groups` output directly (Recommended for Python)
3. Use a simplified format with cluster ID and comma-separated marker genes

You can also use CASSIA's example marker data for testing.

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

#### Marker Data Format
CASSIA accepts three formats:

**1. Seurat FindAllMarkers Output**

Standard output from Seurat's `FindAllMarkers` function with differential expression statistics:

```
p_val  avg_log2FC  pct.1  pct.2  p_val_adj  cluster  gene
0      3.02        0.973  0.152  0          0        CD79A
0      2.74        0.938  0.125  0          0        MS4A1
0      2.54        0.935  0.138  0          0        CD79B
0      1.89        0.812  0.089  0          1        IL7R
0      1.76        0.756  0.112  0          1        CCR7
```

**2. Scanpy rank_genes_groups Output (Recommended for Python)**

Output from Scanpy's `sc.tl.rank_genes_groups()` function, typically exported using `sc.get.rank_genes_groups_df()`:

```
group  names   scores  pvals  pvals_adj  logfoldchanges
0      CD79A   28.53   0      0          3.02
0      MS4A1   25.41   0      0          2.74
0      CD79B   24.89   0      0          2.54
1      IL7R    22.15   0      0          1.89
1      CCR7    20.87   0      0          1.76
```

**3. Simplified Format**

A two-column DataFrame with cluster ID and comma-separated marker genes:

```
cluster  marker_genes
0        CD79A,MS4A1,CD79B,HLA-DRA,TCL1A
1        IL7R,CCR7,LEF1,TCF7,FHIT,MAL
2        CD8A,CD8B,GZMK,CCL5,NKG7
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
    provider = "openrouter",
    reasoning = "medium"  # Optional: "high", "medium", "low" for compatible models
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
- **`reasoning`**: (Optional) Controls reasoning depth ("high", "medium", "low"). See [Reasoning Effort Parameter](setting-up-cassia.md#reasoning-effort-parameter) for details.

### Output Files

The analysis generates three files:

1. **`{output_name}_summary.csv`**: Clean CSV with annotation results including:
   - Cluster ID, Predicted General/Detailed Cell Type, Possible Mixed Cell Types
   - Marker Number and Marker List
   - Iterations, Model, Provider, Tissue, Species metadata

2. **`{output_name}_conversations.json`**: Structured JSON file containing complete conversation history for each cluster:
   - Annotation attempts (all iterations if validation failed)
   - Validation responses and Formatting agent responses
   - Useful for debugging and understanding the reasoning process

3. **`{output_name}_report.html`**: Interactive HTML report with:
   - Summary statistics and visualizations
   - Expandable conversation history for each cluster
   - Searchable and filterable results table

### Tips for Optimal Results

1. **Resource Management**: Monitor system resources when setting `max_workers`. Start with a conservative number if unsure.
2. **Context Optimization**: Use `additional_info` to provide key experimental details (e.g., "tumor sample", "treated condition") to guide the model.
