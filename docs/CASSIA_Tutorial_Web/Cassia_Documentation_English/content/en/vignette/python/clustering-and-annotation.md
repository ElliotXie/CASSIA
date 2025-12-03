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
os.environ["OPENROUTER_API_KEY"] = "your_openrouter_key"
# os.environ["OPENAI_API_KEY"] = "your_openai_key"
# os.environ["ANTHROPIC_API_KEY"] = "your_anthropic_key"
```

## 2. Load Your Data

In this tutorial, we'll use a breast tissue dataset from the GTEX project as an example. This dataset provides a comprehensive reference for breast tissue cell types.

```python
# Load the GTEX breast dataset (assuming .h5ad format)
adata = sc.read("gtex_ref.h5ad")
```

Explore the metadata to understand your dataset:

```python
print(adata.obs.columns)
```

We assume the dataset includes gold standard cell type labels in columns like `Broad_cell_type` and `Granular_cell_type` for reference.

## 3. Dimensional Reduction and Clustering

We'll follow a standard Scanpy workflow with each step explained.

### 3.1 Normalization

First, normalize the data to account for differences in sequencing depth between cells:

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

This scales each cell to 10,000 total counts and applies log transformation to reduce the impact of highly expressed genes.

### 3.2 Feature Selection

Identify highly variable genes that will be used for clustering:

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
```

Subset to only highly variable genes (keep raw data for later visualization):

```python
adata_raw = adata
adata = adata[:, adata.var.highly_variable]
```

### 3.3 Dimensionality Reduction

Scale the data and perform PCA:

```python
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
```

Scaling ensures each gene contributes equally. PCA reduces the data to principal components that capture the most variance.

### 3.4 Neighbor Graph and UMAP

Build a k-nearest neighbors graph and compute UMAP for visualization:

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
sc.tl.umap(adata)
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

![UMAP visualization of GTEX breast tissue clusters](/images/gtex-umap-clusters.png)

## 4. Finding Marker Genes

CASSIA requires a list of marker genes for each cluster. We'll use the Wilcoxon rank-sum test to identify differentially expressed genes.

### 4.1 Run Differential Expression

```python
sc.tl.rank_genes_groups(adata, 'leiden_res_0.4', method='wilcoxon')
```

### 4.2 Extract Results

Extract the results into a DataFrame:

```python
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.head())
```

### 4.3 Format for CASSIA

Rename columns to match CASSIA's expected format:

```python
markers = markers.rename(columns={
    'names': 'gene',
    'logfoldchanges': 'avg_log2FC',
    'pvals_adj': 'p_val_adj',
    'group': 'cluster'
})
```

Filter for positive markers (upregulated genes):

```python
markers = markers[markers['avg_log2FC'] > 0]
```

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
    annotation_model = "anthropic/claude-sonnet-4.5",
    annotation_provider = "openrouter",
    score_model = "openai/gpt-5.1",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model = "anthropic/claude-sonnet-4.5",
    annotationboost_provider = "openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

To see how each model performs, visit our benchmark website: [sc-llm-benchmark.com/methods/cassia](https://sc-llm-benchmark.com/methods/cassia)

### 5.2 Output Files

The pipeline creates a folder with the following files:

| File | Description |
|------|-------------|
| `gtex_breast_annotation_summary.csv` | Summary of annotation results |
| `gtex_breast_annotation_full.csv` | Full results with conversation history |
| `gtex_breast_annotation_scored.csv` | Scored annotation results |
| `gtex_breast_annotation_report.html` | Interactive HTML report |

![CASSIA annotation report](/images/gtex-breast-annotation-report.webp)

### 5.3 Load Results

```python
cassia_results = pd.read_csv("gtex_breast_annotation_scored.csv")
```

### 5.4 Create Annotation Mapping

Create a dictionary mapping cluster IDs to cell type annotations:

```python
annotation_map = dict(zip(
    cassia_results['cluster'].astype(str),
    cassia_results['celltype_1']
))
```

You can choose different granularity levels:
- `celltype_1`: Most detailed annotations
- `CASSIA_merged_grouping_1`: Broadest categories

### 5.5 Add to AnnData

Map the annotations to your AnnData object:

```python
adata.obs['CASSIA_annotation'] = adata.obs['leiden_res_0.4'].map(annotation_map)
```

Verify the mapping:

```python
print(adata.obs[['leiden_res_0.4', 'CASSIA_annotation']].head())
```

## 6. Visualize Annotations

Compare CASSIA annotations with gold standard labels:

```python
sc.pl.umap(adata, color=['Broad_cell_type', 'CASSIA_annotation'], legend_loc='on data')
```

![Figure 1: Gold Standard Clustering](/images/Figure1_GoldStandard.webp)

![Figure 2: Annotation Results](/images/Figure2_CASSIA.webp)

CASSIA automatically summarizes annotations into different levels of granularity, which you can find in the output CSVs.

## 7. Next Steps

After completing clustering and annotation:

- Subset major clusters and repeat the analysis for finer resolution
- Use `CASSIA.runCASSIA_annotationboost()` for clusters with low annotation quality
- Try `CASSIA.symphonyCompare()` to distinguish between similar cell types
- Explore `CASSIA.runCASSIA_subclusters()` for detailed analysis of specific populations
- Implement `CASSIA.runCASSIA_batch_n_times()` for uncertainty quantification

See the "Extended Analysis" vignette for details on these advanced techniques.
