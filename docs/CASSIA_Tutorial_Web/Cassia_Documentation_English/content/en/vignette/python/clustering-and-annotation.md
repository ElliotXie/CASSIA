---
title: "Clustering and Annotation with Scanpy"
---

This vignette demonstrates how to perform clustering on a pre-processed AnnData object using Scanpy and then annotate the resulting clusters using CASSIA. We'll assume quality control steps have already been completed.

## 1. Installation and Setup

### 1.1 Required Packages

```bash
pip install scanpy leidenalg CASSIA
```

### 1.2 Import Packages

```python
import scanpy as sc
import pandas as pd
import CASSIA
import os
```

### 1.3 Set API Keys

**You only need to choose one provider.** OpenRouter is recommended as it provides access to multiple models.

```python
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")
# CASSIA.set_api_key("your-openai-key", provider="openai")
# CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
```

## 2. Load Your Data

In this tutorial, we'll use a breast tissue dataset from the GTEX project as an example. This dataset provides a comprehensive reference for breast tissue cell types.

Download the dataset: [GTEx_breast_minimal.h5ad](https://drive.google.com/file/d/1HhGX0AD6tfUYzuYc2LToxghTGk8LZeUp/view?usp=sharing)

```python
# Load the GTEX breast dataset (assuming .h5ad format)
adata = sc.read("GTEx_breast_minimal.h5ad")
```

Explore the metadata to understand your dataset:

```python
print(adata.obs.columns)
```

This dataset includes gold standard cell type labels in columns `Broad cell type` and `Granular cell type` for reference.

## 3. Dimensional Reduction and Clustering

We'll follow a standard Scanpy workflow with each step explained.

### 3.1 Normalization

First, normalize the data to account for differences in sequencing depth between cells:

```python
# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

### 3.2 Feature Selection

Identify highly variable genes that will be used for clustering:

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
```

### 3.3 Dimensionality Reduction

```python
sc.tl.pca(adata)
```

PCA reduces the data to principal components that capture the most variance.

### 3.4 Neighbor Graph and UMAP

Build a k-nearest neighbors graph and compute UMAP for visualization:

```python
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

```python
sc.pl.umap(
    adata,
    color="Broad cell type",
    size=2,
)
```

The neighbor graph connects similar cells and forms the basis for clustering.

### 3.5 Clustering

Perform Leiden clustering to identify cell populations:

```python
sc.tl.leiden(adata, resolution=0.4, key_added='leiden_res_0.4')
```

The `resolution` parameter controls cluster granularity. Lower values = fewer, broader clusters.

### 3.6 Visualize Clusters

```python
sc.pl.umap(adata, color=['leiden_res_0.4'])
```

## 4. Finding Marker Genes

CASSIA requires a list of marker genes for each cluster. We'll use the Wilcoxon rank-sum test to identify differentially expressed genes.

### 4.1 Run Differential Expression

```python
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.4", method="wilcoxon",use_raw=False)
```

### 4.2 Extract Enhanced Markers

Use CASSIA's `enhance_scanpy_markers()` to extract marker genes with percentage expression values (`pct.1` and `pct.2`). This function reads directly from `adata.uns['rank_genes_groups']` and adds important metadata:

```python
# Extract enhanced markers with pct.1/pct.2 values
markers = CASSIA.enhance_scanpy_markers(adata, cluster_col="leiden_res_0.4", n_genes=50)
print(markers.head())
```

The returned DataFrame includes:
- `pct.1`: Percentage of cells in the cluster expressing each gene
- `pct.2`: Percentage of cells outside the cluster expressing each gene
- Standard statistics: `avg_log2FC`, `p_val_adj`, `scores`

These percentage values help CASSIA better understand marker specificity and are particularly useful for annotation boost.

## 5. Annotating Clusters with CASSIA

### 5.1 Run the CASSIA Pipeline

Now we can annotate our clusters using CASSIA:

```python
results = CASSIA.runCASSIA_pipeline(
    output_file_name = "gtex_breast_annotation",
    tissue = "Breast",
    species = "Human",
    marker = markers,
    max_workers = 6,
    annotation_model = "openai/gpt-5.1",
    annotation_provider = "openrouter",
    score_model = "anthropic/claude-sonnet-4.5",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model = "openai/gpt-5.1",
    annotationboost_provider = "openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

To see how each model performs, visit our benchmark website: [sc-llm-benchmark.com/methods/cassia](https://sc-llm-benchmark.com/methods/cassia)

### 5.2 Output Files

The pipeline creates an output folder named `CASSIA_Pipeline_{tissue}_{species}_{timestamp}/` with three subfolders:

- `01_annotation_report/` - Interactive HTML reports for the analysis
- `02_annotation_boost/` - Annotation boost results for low-scoring clusters
- `03_csv_files/` - Summary CSV files including final results

### 5.3 Load Results

```python
# Replace the folder name with the actual output folder (includes timestamp)
cassia_results = pd.read_csv("CASSIA_Pipeline_Breast_Human_XXXXXX/03_csv_files/gtex_breast_annotation_FINAL_RESULTS.csv")
```

### 5.4 Add Annotations to AnnData

Use CASSIA's `add_cassia_to_anndata()` to automatically integrate all annotation results into your AnnData object:

```python
# Add CASSIA annotations with automatic cluster matching
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    cluster_col="leiden_res_0.4",
    prefix="CASSIA_"
)
```

This function automatically:
- Matches cluster IDs between CASSIA results and AnnData (with fuzzy matching)
- Adds multiple annotation columns to `adata.obs`
- Stores cluster-level summary in `adata.uns['CASSIA']`

### 5.5 View Annotations

The function adds several useful columns to `adata.obs`:

```python
# View added columns
cassia_cols = [col for col in adata.obs.columns if col.startswith('CASSIA_')]
print("Added columns:", cassia_cols)

# View sample annotations
print(adata.obs[['leiden_res_0.4', 'CASSIA_general_celltype', 'CASSIA_sub_celltype', 'CASSIA_score']].head(10))

# View cluster-level summary
print("\nCluster-level summary:")
print(adata.uns['CASSIA'])
```

Available annotation columns include:
- `CASSIA_general_celltype`: Broad cell type annotations
- `CASSIA_sub_celltype`: Detailed sub-celltype annotations
- `CASSIA_score`: Confidence score (0-100)
- `CASSIA_combined_celltype`: Format "General :: Subtype"
- Additional columns for merged groupings and alternative predictions

## 6. Visualize Annotations

Compare CASSIA annotations with gold standard labels:

```python
# Visualize CASSIA annotations alongside original clusters
sc.pl.umap(adata, color=['leiden_res_0.4', 'CASSIA_general_celltype'], ncols=2, size=2)

# Compare with gold standard labels (if available)
sc.pl.umap(adata, color=['CASSIA_general_celltype', 'Broad cell type'], ncols=2, size=2)
```

You can also create summary statistics:

```python
# Count cells per cell type
celltype_counts = adata.obs['CASSIA_general_celltype'].value_counts()
print(celltype_counts)

# Compare cluster assignments with CASSIA annotations
comparison = adata.obs.groupby(['leiden_res_0.4', 'CASSIA_general_celltype']).size().unstack(fill_value=0)
print(comparison)
```

CASSIA automatically summarizes annotations into different levels of granularity, which you can find in the output CSVs and the various `CASSIA_*` columns in `adata.obs`.

## 7. Next Steps

After completing clustering and annotation:

- Subset major clusters and repeat the analysis for finer resolution
- Use `CASSIA.runCASSIA_annotationboost()` for clusters with low annotation quality
- Try `CASSIA.symphonyCompare()` to distinguish between similar cell types
- Explore `CASSIA.runCASSIA_subclusters()` for detailed analysis of specific populations
- Implement `CASSIA.runCASSIA_batch_n_times()` for uncertainty quantification

See the "Extended Analysis" vignette for details on these advanced techniques.
