---
title: "Clustering and Annotation with Seurat"
---

This vignette demonstrates how to perform clustering on a pre-processed Seurat object and then annotate the resulting clusters using CASSIA. We'll assume quality control steps have already been completed.

## 1. Installation and Setup

### 1.1 Required Packages

```r
# Install required packages
install.packages(c("Seurat", "dplyr", "reticulate", "devtools"))

# Install CASSIA from GitHub
library(devtools)
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")

# Load packages
library(Seurat)
library(dplyr)
library(CASSIA)

# Set up API key (needed for CASSIA)
setLLMApiKey("your_api_key", provider = "openrouter", persist = TRUE)  # or use anthropic/openai
```

## 2. Exploring Your Pre-processed Seurat Object

In this tutorial, we'll use a breast tissue dataset from the GTEX project as an example. This dataset provides a comprehensive reference for breast tissue cell types and can be downloaded from:

[***GTEX Breast Tissue Dataset***](https://figshare.com/articles/dataset/processed_breast_RDS/28941470?file=54252515)

```r
# Load the GTEX breast dataset
gtex_data <- readRDS("gtex_ref.rds")

# Take a look at the dataset metadata
gtex_data@meta.data %>% colnames()
```

>[1] "orig.ident"         "nCount_RNA"         "nFeature_RNA"       "Broad.cell.type"    "Granular.cell.type"

We can see that the dataset includes gold standard cell type labels in the `Broad.cell.type` and `Granular.cell.type` columns.

## 3. Dimensional Reduction and Clustering

### 3.1 PCA and UMAP

```r
gtex_data <- NormalizeData(object = gtex_data)
gtex_data <- FindVariableFeatures(object = gtex_data)
gtex_data <- ScaleData(object = gtex_data)
gtex_data <- RunPCA(object = gtex_data)
gtex_data <- FindNeighbors(object = gtex_data, dims = 1:25)     # Use 25 PCs for clustering
gtex_data <- FindClusters(object = gtex_data,resolution = 0.4)  # Use resolution 0.4 for clustering
gtex_data <- RunUMAP(object = gtex_data, dims = 1:25)           # Use 25 PCs for UMAP
DimPlot(object = gtex_data, reduction = "umap")                # Visualize the UMAP
```


![UMAP visualization of GTEX breast tissue clusters](/images/gtex-umap-clusters.png)

### 3.3 Finding Marker Genes

```r

Idents(gtex_data)="RNA_snn_res.0.4" # This is the default cluster name if you use the default resolution of 0.4


all_markers <- FindAllMarkers(gtex_data, 
                              only.pos = TRUE) 

```
The marker filter is implemented internally in CASSIA, and the raw marker file is needed for the annotation boost agent too; thus, no filter is needed at this step.

## 4. Annotating Clusters with CASSIA

### 4.1 Basic Annotation

Now that we have our clusters and marker genes, we can use CASSIA for annotation:


```r
results <- runCASSIA_pipeline(
    output_file_name = "gtex_breast_annotation",
    tissue = "Breast",
    species = "Human",
    marker = all_markers,
    max_workers = 6  # Depends on your CPU cores
)

```
The default provider is OpenRouter, and the default models are selected to optimize for annotation quality. To see how each new model performs, visit our benchmark website:

[***sc-llm-benchmark.com/methods/cassia***](https://sc-llm-benchmark.com/methods/cassia)


The output file is saved in a folder named after the tissue and species. Inside the folder, you will see the following files:

- `gtex_breast_annotation_summary.csv`: The summary of the annotation results
- `gtex_breast_annotation_full.csv`: The full annotation results with full conversation history and merged clusters at different levels
- `gtex_breast_annotation_scored.csv`: The scored annotation results
- `gtex_breast_annotation_report.html`: Genearted report of the annotation results, contains the router to all the clusters reports.

The gtex_breast_annotation_scored.csv file:

![CASSIA annotation report](/images/gtex-breast-annotation-report.png)



### 4.2 Incorporating Annotations into Seurat Object

After running CASSIA, you can integrate the annotations back into your Seurat object:

```r
seurat_breast=add_cassia_to_seurat(
    seurat_obj = gtex_data, 
    cassia_results_path = "gtex_breast_annotation_scored.csv", # the path to the CASSIA results csv file
    cluster_col = "RNA_snn_res.0.4", # the column name of the cluster column used for Findallmarkers in the seurat object
    columns_to_include = "1") # default is to export only the merged column back to seurat, set to "2" if you want to export all the columns
```


### 4.3 Visualizing the Annotations

```r
DimPlot(seurat_breast,group.by = "Broad.cell.type", label = TRUE, repel = TRUE)
DimPlot(seurat_breast,group.by = "Granular.cell.type", label = TRUE, repel = TRUE)
DimPlot(seurat_breast,group.by = "RNA_snn_res.0.4", label = TRUE, repel = TRUE)
DimPlot(seurat_breast,group.by = "CASSIA_merged_grouping_1", label = TRUE, repel = TRUE)
DimPlot(seurat_breast,group.by = "CASSIA_merged_grouping_2", label = TRUE, repel = TRUE)
DimPlot(seurat_breast,group.by = "CASSIA_merged_grouping_3", label = TRUE, repel = TRUE)

```



![Figure 1: Gold Standard Clustering|width=50%,height=50%](/images/Figure1_GoldStandard.png)

![Figure 2: Annotation Results](/images/Figure2_CASSIA.png)



CASSIA will automatically summarize the annotation into three different detailed levels, from the most general to the most detailed, as shown in the UMAPs.





## 5. Next Steps

After completing clustering and annotation:

- Subset the major clusters and repeat the steps
- Use `runCASSIA_annotationboost()` for clusters with low annotation quality
- Try `symphonyCompare()` to distinguish between similar cell types
- Explore `runCASSIA_subclusters()` for more detailed analysis of specific populations
- Implement `runCASSIA_batch_n_times()` for uncertainty quantification