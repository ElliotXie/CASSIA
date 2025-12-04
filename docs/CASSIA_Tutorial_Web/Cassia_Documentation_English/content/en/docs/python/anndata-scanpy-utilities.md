---
title: "AnnData and Scanpy Integration Utilities"
---

CASSIA provides utility functions to streamline integration with Scanpy workflows and AnnData objects. These functions simplify marker extraction and annotation integration, making it easier to work with single-cell analysis pipelines.

## Overview

Two utility functions are available for Scanpy/AnnData integration:

1. **`enhance_scanpy_markers()`**: Extracts marker genes from Scanpy's differential expression results and adds `pct.1` and `pct.2` values (percentage of cells expressing each gene within/outside each cluster).

2. **`add_cassia_to_anndata()`**: Automatically adds CASSIA annotation results to your AnnData object with fuzzy cluster matching and comprehensive metadata.

## Installation

These functions require optional dependencies:

```bash
pip install scanpy anndata
```

Import in your Python script:

```python
import scanpy as sc
import CASSIA
```

## enhance_scanpy_markers()

### Description

Extracts and enhances Scanpy marker genes with percentage expression values (`pct.1` and `pct.2`). This function is designed for workflows that use Scanpy's `rank_genes_groups` for differential expression analysis.

**IMPORTANT:** This function takes an **AnnData object** as input, NOT a DataFrame. It reads directly from `adata.uns['rank_genes_groups']` (the results stored by `sc.tl.rank_genes_groups()`) and returns a new enhanced DataFrame. It is a **REPLACEMENT** for `sc.get.rank_genes_groups_df()`, not a post-processor for its output.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | *required* | Annotated data matrix with `rank_genes_groups` results in `.uns` |
| `cluster_col` | str or None | None | Column name in `adata.obs` containing cluster assignments. If None, auto-detects (tries 'leiden', then 'louvain') |
| `n_genes` | int or None | None | Number of top genes per cluster to include. If None, includes all genes |
| `min_expression` | float | 0.0 | Threshold above which a cell is considered "expressing" a gene |
| `include_stats` | bool | True | Whether to include additional statistics (logfoldchanges, pvals, scores) |
| `key` | str | "rank_genes_groups" | Key in `adata.uns` where `rank_genes_groups` results are stored |

### Returns

Returns a pandas DataFrame with the following columns:

| Column | Description |
|--------|-------------|
| `cluster` | Cluster ID |
| `gene` | Gene name |
| `pct.1` | Fraction of cells in cluster expressing gene (0.0-1.0) |
| `pct.2` | Fraction of cells outside cluster expressing gene (0.0-1.0) |
| `avg_log2FC` | Log fold change (if `include_stats=True` and available) |
| `p_val_adj` | Adjusted p-value (if `include_stats=True` and available) |
| `scores` | Scanpy score (if `include_stats=True` and available) |

### Basic Usage

```python
import scanpy as sc
import CASSIA

# Load data and perform clustering
adata = sc.read_h5ad("your_data.h5ad")
sc.tl.leiden(adata, resolution=0.5)

# Run differential expression
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# Extract enhanced markers
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)
print(markers.head())
```

### Advanced Usage

#### Custom cluster column and expression threshold

```python
# Use custom cluster column and expression threshold
markers = CASSIA.enhance_scanpy_markers(
    adata,
    cluster_col="my_custom_clusters",
    n_genes=100,
    min_expression=0.1,  # Only count cells with expression > 0.1
    include_stats=True
)
```

#### Integration with CASSIA annotation

```python
# Extract markers and run CASSIA annotation
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)

# Run CASSIA batch annotation
results = CASSIA.runCASSIA_batch(
    marker=markers,
    tissue="brain",
    species="mouse",
    output_name="brain_annotation"
)
```

## add_cassia_to_anndata()

### Description

Integrates CASSIA annotation results into an AnnData object. This function automatically:
- Matches cluster IDs between CASSIA results and AnnData (with fuzzy matching support)
- Adds multiple annotation columns to `adata.obs`
- Stores cluster-level summary in `adata.uns['CASSIA']`
- Handles missing or mismatched cluster names gracefully

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | *required* | Annotated data matrix |
| `cassia_results` | str or DataFrame | *required* | Path to CASSIA results CSV file OR DataFrame from `runCASSIA_batch()` |
| `cluster_col` | str or None | None | Column name in `adata.obs` containing cluster assignments. If None, auto-detects (tries 'leiden', then 'louvain') |
| `cassia_cluster_col` | str | "Cluster ID" | Column name in CASSIA results for cluster IDs |
| `prefix` | str | "CASSIA_" | Prefix for new column names |
| `replace_existing` | bool | False | Whether to overwrite existing CASSIA columns |
| `fuzzy_match` | bool | True | Enable fuzzy matching for cluster ID alignment (handles variations like "0" vs "cluster_0") |
| `columns_to_include` | int | 2 | 1 = merged groupings only, 2 = all metrics |
| `inplace` | bool | True | If True, modify adata in place. If False, return copy |

### Output

Adds the following columns to `adata.obs` (with specified prefix):

| Column | Description |
|--------|-------------|
| `CASSIA_general_celltype` | Main cell type prediction |
| `CASSIA_sub_celltype` | Primary sub-celltype (first from ranked list) |
| `CASSIA_sub_celltype_all` | Full comma-separated string of sub-celltypes |
| `CASSIA_sub_celltype_1/2/3` | Split ranked alternative sub-celltypes |
| `CASSIA_mixed_celltype` | Possible mixed cell types |
| `CASSIA_score` | Quality/consensus score (0-100) |
| `CASSIA_merged_grouping_1/2/3` | Hierarchical groupings (if available) |
| `CASSIA_combined_celltype` | Format "General :: Subtype" |

Also stores cluster-level summary in `adata.uns['CASSIA']`.

### Basic Usage

```python
import scanpy as sc
import CASSIA

# Load AnnData and CASSIA results
adata = sc.read_h5ad("your_data.h5ad")
cassia_results = "brain_annotation_FINAL_RESULTS.csv"

# Add CASSIA annotations
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    cluster_col="leiden"
)

# View annotations
print(adata.obs[['leiden', 'CASSIA_general_celltype', 'CASSIA_sub_celltype']].head())

# View cluster-level summary
print(adata.uns['CASSIA'])
```

### Advanced Usage

#### Using DataFrame directly from runCASSIA_batch

```python
# Run CASSIA and add results in one workflow
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)

results_df = CASSIA.runCASSIA_batch(
    marker=markers,
    tissue="lung",
    species="human",
    output_name="lung_annotation"
)

# Add annotations directly from DataFrame (no need to read CSV)
CASSIA.add_cassia_to_anndata(
    adata,
    results_df,
    cluster_col="leiden"
)
```

#### Custom prefix and column selection

```python
# Add only merged groupings with custom prefix
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    cluster_col="seurat_clusters",
    prefix="CellType_",
    columns_to_include=1  # Only merged groupings
)
```

#### Handling mismatched cluster names

```python
# Fuzzy matching handles variations automatically
# e.g., "0" in adata matches "cluster_0" in CASSIA results
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    fuzzy_match=True  # Default is True
)
```

## Complete Workflow Example

Here's a complete workflow from clustering to annotation:

```python
import scanpy as sc
import CASSIA

# Set API key
CASSIA.set_api_key("your-api-key", provider="openrouter")

# Load and preprocess data
adata = sc.read_h5ad("pbmc3k.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Dimensionality reduction and clustering
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Find marker genes (stores in adata.uns)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# Extract enhanced markers with pct.1/pct.2
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)
print(f"Extracted {len(markers)} marker genes")

# Run CASSIA annotation
results = CASSIA.runCASSIA_batch(
    marker=markers,
    tissue="blood",
    species="human",
    output_name="pbmc_annotation",
    model="anthropic/claude-sonnet-4.5",
    provider="openrouter"
)

# Add annotations to AnnData
CASSIA.add_cassia_to_anndata(adata, results)

# Visualize
sc.pl.umap(adata, color=['leiden', 'CASSIA_general_celltype'], ncols=2)

# Access annotations
print(adata.obs[['leiden', 'CASSIA_general_celltype', 'CASSIA_score']].head(20))
print("\nCluster-level summary:")
print(adata.uns['CASSIA'])
```

## Troubleshooting

### Error: "rank_genes_groups results not found in adata.uns"

**Cause:** `enhance_scanpy_markers()` was called before running `sc.tl.rank_genes_groups()`.

**Solution:** Make sure to run differential expression first:

```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
markers = CASSIA.enhance_scanpy_markers(adata)
```

### Error: "Could not auto-detect cluster column"

**Cause:** Neither 'leiden' nor 'louvain' columns exist in `adata.obs`.

**Solution:** Specify the cluster column explicitly:

```python
markers = CASSIA.enhance_scanpy_markers(adata, cluster_col="my_clusters")
```

### Warning: "Could not find matches for clusters"

**Cause:** Cluster names in AnnData don't match those in CASSIA results.

**Solution:** Check cluster naming and enable fuzzy matching (default):

```python
# View cluster names
print("AnnData clusters:", adata.obs['leiden'].unique())
print("CASSIA clusters:", cassia_results['Cluster ID'].unique())

# Fuzzy matching is enabled by default
CASSIA.add_cassia_to_anndata(adata, cassia_results, fuzzy_match=True)
```

### Error: "CASSIA columns already exist"

**Cause:** Trying to add annotations when columns already exist.

**Solution:** Set `replace_existing=True` to overwrite:

```python
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    replace_existing=True
)
```

### Missing pct.1 or pct.2 values (NaN)

**Cause:** Gene not found in the AnnData object.

**Solution:** Verify gene names match between marker results and AnnData:

```python
# Check if genes exist
missing_genes = [g for g in markers['gene'].unique() if g not in adata.var_names]
print(f"Missing genes: {len(missing_genes)}")
```

## Tips for Best Results

1. **Use sufficient marker genes**: 50-100 genes per cluster typically provides good annotation quality.
2. **Check expression thresholds**: Adjust `min_expression` if your data has different scaling.
3. **Verify cluster matching**: Always check that cluster IDs align between your analysis and CASSIA results.
4. **Save your work**: After adding annotations, save the AnnData object:
   ```python
   adata.write_h5ad("annotated_data.h5ad")
   ```
