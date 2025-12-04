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

### 4.2 Extract Results

Extract the results into a DataFrame:

```python
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.head())
```

```
  group     names     scores  logfoldchanges          pvals      pvals_adj
0     0     KRT15  48.403599        5.089063   0.000000e+00   0.000000e+00
1     0   TFCP2L1  36.434315        4.719573  1.218868e-290  1.078393e-286
2     0       KIT  35.415207        4.628579  9.961057e-275  5.875364e-271
3     0      NFIB  34.773083        2.135615  6.208985e-265  2.746700e-261
4     0  ANKRD36C  34.635426        3.357123  7.404516e-263  2.620458e-259
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

The pipeline creates an output folder with three subfolders:

- `01_html_reports/` - Interactive HTML reports for each cluster
- `02_cluster_annotations/` - Individual annotation results per cluster
- `03_csv_files/` - Summary CSV files including final results

### 5.3 Load Results

```python
cassia_results = pd.read_csv("CASSIA_Pipeline_output/gtex_breast_annotation_FINAL_RESULTS.csv")
```

### 5.4 Create Annotation Mapping

Create a dictionary mapping cluster IDs to cell type annotations:

```python
annotation_map = dict(zip(
    cassia_results['Cluster ID'].astype(str),
    cassia_results['Predicted General Cell Type']
))
```

You can choose different granularity levels:
- `Predicted General Cell Type`: Broad annotations

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

CASSIA automatically summarizes annotations into different levels of granularity, which you can find in the output CSVs.

## 7. Next Steps

After completing clustering and annotation:

- Subset major clusters and repeat the analysis for finer resolution
- Use `CASSIA.runCASSIA_annotationboost()` for clusters with low annotation quality
- Try `CASSIA.symphonyCompare()` to distinguish between similar cell types
- Explore `CASSIA.runCASSIA_subclusters()` for detailed analysis of specific populations
- Implement `CASSIA.runCASSIA_batch_n_times()` for uncertainty quantification

See the "Extended Analysis" vignette for details on these advanced techniques.
