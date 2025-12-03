---
title: "Clustering and Annotation with Scanpy"
---

This vignette demonstrates how to perform clustering on a pre-processed AnnData object using Scanpy and then annotate the resulting clusters using CASSIA. We'll assume quality control steps have already been completed.

## 1. Installation and Setup

### 1.1 Required Packages

```bash
# Install required packages
pip install scanpy leidenalg CASSIA
```

### 1.2 Import Packages

```python
import scanpy as sc
import pandas as pd
import CASSIA
import os

# Set up API keys (needed for CASSIA)
# Replace with your actual keys
os.environ["OPENROUTER_API_KEY"] = "your_openrouter_key"
# os.environ["OPENAI_API_KEY"] = "your_openai_key"
# os.environ["ANTHROPIC_API_KEY"] = "your_anthropic_key"
```

## 2. Exploring Your Pre-processed AnnData Object

In this tutorial, we'll use a breast tissue dataset from the GTEX project as an example. This dataset provides a comprehensive reference for breast tissue cell types.

```python
# Load the GTEX breast dataset (assuming .h5ad format)
# You can adapt this to read_csv or read_10x depending on your data
adata = sc.read("gtex_ref.h5ad")

# Check the metadata
print(adata.obs.columns)
```

We assume the dataset includes gold standard cell type labels in columns like `Broad_cell_type` and `Granular_cell_type` for reference.

## 3. Dimensional Reduction and Clustering

### 3.1 Preprocessing and Clustering

We'll follow a standard Scanpy workflow: normalization, log-transformation, feature selection, scaling, PCA, neighbor graph construction, and clustering.

```python
# Normalize and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_raw = adata  # keep raw data for visualization if needed
adata = adata[:, adata.var.highly_variable]

# Scale and PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# Neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
sc.tl.umap(adata)

# Clustering (Leiden)
# Use resolution 0.4 to match the R example
sc.tl.leiden(adata, resolution=0.4, key_added='leiden_res_0.4')

# Visualize
sc.pl.umap(adata, color=['leiden_res_0.4'])
```

![UMAP visualization of GTEX breast tissue clusters](/images/gtex-umap-clusters.png)

### 3.2 Finding Marker Genes

CASSIA requires a list of marker genes for each cluster. We can generate this using `scanpy.tl.rank_genes_groups`.

```python
# Find markers for the clusters
sc.tl.rank_genes_groups(adata, 'leiden_res_0.4', method='wilcoxon')

# Extract markers into a DataFrame
markers = sc.get.rank_genes_groups_df(adata, group=None)

# Rename columns to match CASSIA's expected format (Seurat-like)
# Scanpy outputs: names, scores, logfoldchanges, pvals, pvals_adj
markers = markers.rename(columns={
    'names': 'gene',
    'logfoldchanges': 'avg_log2FC', 
    'pvals_adj': 'p_val_adj', 
    'group': 'cluster'
})

# Filter for positive markers (optional, CASSIA handles this but good practice)
markers = markers[markers['avg_log2FC'] > 0]

print(markers.head())
```

## 4. Annotating Clusters with CASSIA

### 4.1 Basic Annotation

Now that we have our clusters and marker genes, we can use CASSIA for annotation:

```python
# Run the CASSIA pipeline
results = CASSIA.runCASSIA_pipeline(
    output_file_name = "gtex_breast_annotation",
    tissue = "Breast",
    species = "Human",
    marker_path = markers, # Pass the DataFrame directly
    max_workers = 6,  # Matches the number of clusters in dataset
    annotation_model = "anthropic/claude-sonnet-4.5",
    annotation_provider = "openrouter",
    score_model = "openai/gpt-5.1",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model="anthropic/claude-sonnet-4.5",
    annotationboost_provider="openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

The default provider is OpenRouter, and the default models are selected to optimize for annotation quality. To see how each new model performs, visit our benchmark website: [sc-llm-benchmark.com/methods/cassia](https://sc-llm-benchmark.com/methods/cassia)

The output file is saved in a folder named after the tissue and species. Inside the folder, you will see the following files:

- `gtex_breast_annotation_summary.csv`: The summary of the annotation results
- `gtex_breast_annotation_full.csv`: The full annotation results with full conversation history and merged clusters at different levels
- `gtex_breast_annotation_scored.csv`: The scored annotation results
- `gtex_breast_annotation_report.html`: Generated report of the annotation results, contains the router to all the clusters reports.

The `gtex_breast_annotation_scored.csv` file:

![CASSIA annotation report](/images/gtex-breast-annotation-report.png)

### 4.2 Incorporating Annotations into AnnData Object

After running CASSIA, you can integrate the annotations back into your AnnData object. We'll map the annotations from the CSV results to the cluster observations.

```python
# Load the results
cassia_results = pd.read_csv("gtex_breast_annotation_scored.csv")

# Create a mapping dictionary: Cluster ID -> Annotation
# You can choose 'celltype_1' (most detailed) or 'CASSIA_merged_grouping_1' (broadest)
annotation_map = dict(zip(cassia_results['cluster'].astype(str), cassia_results['celltype_1']))

# Map the annotations to a new observation column
adata.obs['CASSIA_annotation'] = adata.obs['leiden_res_0.4'].map(annotation_map)

# Inspect the new annotations
print(adata.obs[['leiden_res_0.4', 'CASSIA_annotation']].head())
```

### 4.3 Visualizing the Annotations

Now we can visualize the CASSIA annotations on the UMAP.

```python
sc.pl.umap(adata, color=['Broad_cell_type', 'CASSIA_annotation'], legend_loc='on data')
```

![Figure 1: Gold Standard Clustering](/images/Figure1_GoldStandard.png)

![Figure 2: Annotation Results](/images/Figure2_CASSIA.png)

CASSIA automatically summarizes annotations into different levels of granularity (from general to detailed), which you can find in the output CSVs and visualize similarly.

## 5. Next Steps

After completing clustering and annotation:

- Subset the major clusters and repeat the steps
- Use `CASSIA.runCASSIA_annotationboost()` for clusters with low annotation quality
- Try `CASSIA.symphonyCompare()` to distinguish between similar cell types
- Explore `CASSIA.runCASSIA_subclusters()` for more detailed analysis of specific populations
- Implement `CASSIA.runCASSIA_batch_n_times()` for uncertainty quantification

See the "Extended Analysis" vignette for details on these advanced techniques.

