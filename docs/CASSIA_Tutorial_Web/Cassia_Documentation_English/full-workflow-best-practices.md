---
title: "Full Workflow Best Practices for Single-Cell RNA-seq Analysis"
---

This vignette provides a comprehensive guide to best practices for single-cell RNA-seq analysis using CASSIA, including rigorous quality control steps, doublet removal, background correction, and advanced annotation techniques.

## 1. Comprehensive Setup

### 1.1 Required Packages and Environment

```r
# Install required packages
install.packages(c("Seurat", "dplyr", "patchwork", "ggplot2", "reticulate", "devtools", "DoubletFinder"))
BiocManager::install(c("scDblFinder", "scater", "scran"))

# Install CASSIA from GitHub
library(devtools)
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")

# Load essential packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(scDblFinder)
library(CASSIA)

# Set up CASSIA environment and API key
setup_cassia_env(conda_env = "cassia_env")
setLLMApiKey("your_api_key", provider = "openai", persist = TRUE)  # Or use anthropic/openrouter
```

## 2. Data Preparation and Quality Control

### 2.1 Loading Data and Initial QC

```r
# Load data (10X format example)
counts <- Read10X(data.dir = "path/to/filtered_gene_bc_matrices/")
seurat_obj <- CreateSeuratObject(counts = counts, project = "scRNA_project", min.cells = 3, min.features = 200)

# Calculate mitochondrial and ribosomal percentages
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Comprehensive QC visualization
qc_violin <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
qc_scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
qc_scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc_scatter3 <- FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
qc_violin / (qc_scatter1 | qc_scatter2 | qc_scatter3)

# Determine QC thresholds
# Instead of using fixed thresholds, use MAD (median absolute deviation) approach
nCount_mad <- median(seurat_obj$nCount_RNA) + 3*mad(seurat_obj$nCount_RNA)
nFeature_mad_min <- median(seurat_obj$nFeature_RNA) - 3*mad(seurat_obj$nFeature_RNA)
nFeature_mad_max <- median(seurat_obj$nFeature_RNA) + 3*mad(seurat_obj$nFeature_RNA)
mt_mad <- median(seurat_obj$percent.mt) + 2*mad(seurat_obj$percent.mt)

# Create a more detailed QC plot with thresholds
ggplot(seurat_obj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c() +
  geom_hline(yintercept = nFeature_mad_min, linetype = "dashed", color = "red") +
  geom_hline(yintercept = nFeature_mad_max, linetype = "dashed", color = "red") +
  geom_vline(xintercept = nCount_mad, linetype = "dashed", color = "red") +
  labs(title = "QC Metrics with Adaptive Thresholds", 
       x = "Number of UMIs", 
       y = "Number of Genes",
       color = "% MT Genes")

# Filter cells based on adaptive thresholds
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > nFeature_mad_min & 
                              nFeature_RNA < nFeature_mad_max & 
                              nCount_RNA < nCount_mad & 
                              percent.mt < mt_mad)
```

### 2.2 Detecting and Removing Doublets

Multiple doublet detection methods for increased confidence:

```r
# Initial preprocessing for doublet detection
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# Method 1: DoubletFinder
# pK identification (no ground-truth)
sweep.res <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]

# Run DoubletFinder with optimal pK
nExp <- round(0.08 * ncol(seurat_obj))  # Assuming ~8% doublet rate for most platforms
seurat_obj <- doubletFinder_v3(seurat_obj, 
                              PCs = 1:20, 
                              pN = 0.25, 
                              pK = optimal_pk, 
                              nExp = nExp, 
                              reuse.pANN = FALSE, 
                              sct = FALSE)

# Get the DoubletFinder results column name (it has a variable name)
df_column <- colnames(seurat_obj@meta.data)[grepl("DF.classification", colnames(seurat_obj@meta.data))]

# Method 2: scDblFinder
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce)
seurat_obj$scDblFinder_class <- sce$scDblFinder.class

# Combine results from both methods for higher confidence
seurat_obj$doublet_consensus <- "Singlet"
seurat_obj$doublet_consensus[seurat_obj@meta.data[,df_column] == "Doublet" | 
                           seurat_obj$scDblFinder_class == "doublet"] <- "Doublet"

# Visualize doublets on UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, group.by = "doublet_consensus", cols = c("Singlet" = "lightgrey", "Doublet" = "red"))

# Remove consensus doublets
seurat_clean <- subset(seurat_obj, subset = doublet_consensus == "Singlet")
```

### 2.3 Ambient RNA Correction

```r
# For 10X data, you can use SoupX for ambient RNA correction
# This requires the raw (unfiltered) matrix
library(SoupX)

# Convert to SoupChannel
sc <- SoupChannel(tod = counts, toc = raw_counts)  # Assuming you have raw counts
sc <- setClusters(sc, Idents(seurat_obj))
sc <- autoEstCont(sc)
corrected_counts <- adjustCounts(sc)

# Create a new Seurat object with corrected counts
seurat_clean <- CreateSeuratObject(counts = corrected_counts, 
                                  meta.data = seurat_clean@meta.data)
```

## 3. Advanced Preprocessing and Clustering

### 3.1 Normalization and Feature Selection

```r
# SCTransform normalization (alternative to standard normalization)
seurat_clean <- SCTransform(seurat_clean, 
                           vars.to.regress = c("percent.mt", "percent.rb"),
                           verbose = FALSE)

# Dimensional reduction
seurat_clean <- RunPCA(seurat_clean, verbose = FALSE)

# Determine optimal number of PCs using Elbow plot
ElbowPlot(seurat_clean, ndims = 50)

# Alternative: More rigorous PC selection using jackstraw method
seurat_clean <- JackStraw(seurat_clean, num.replicate = 100, dims = 40)
seurat_clean <- ScoreJackStraw(seurat_clean, dims = 1:40)
JackStrawPlot(seurat_clean, dims = 1:40)
```

### 3.2 Batch Correction (If Applicable)

If your dataset contains multiple batches:

```r
# Option 1: Integration with Seurat (for multiple samples)
# Assuming seurat_obj_list is a list of Seurat objects from different batches
seurat_list <- SplitObject(seurat_clean, split.by = "batch")
seurat_list <- lapply(seurat_list, SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Option 2: Harmony batch correction
library(harmony)
seurat_clean <- RunHarmony(seurat_clean, group.by.vars = "batch")
```

### 3.3 Clustering Strategy

```r
# Generate a range of clustering resolutions
resolutions <- seq(0.2, 1.2, by = 0.2)
for(res in resolutions) {
  seurat_clean <- FindNeighbors(seurat_clean, dims = 1:30)
  seurat_clean <- FindClusters(seurat_clean, resolution = res, 
                              algorithm = 1, # Louvain algorithm
                              random.seed = 42)
  # Rename the resulting clusters to indicate resolution
  seurat_clean$clusters <- paste0("res", res, "_", seurat_clean@meta.data[,paste0("SCT_snn_res.", res)])
}

# Visualize different clustering resolutions
plot_list <- list()
for(res in resolutions) {
  column_name <- paste0("SCT_snn_res.", res)
  plot_list[[as.character(res)]] <- DimPlot(seurat_clean, 
                                          group.by = column_name, 
                                          label = TRUE) + 
                                    ggtitle(paste0("Resolution: ", res))
}
wrap_plots(plot_list, ncol = 3)

# Choose an optimal resolution based on biological knowledge and cluster stability
# For this example, we'll use 0.6
seurat_clean <- FindNeighbors(seurat_clean, dims = 1:30)
seurat_clean <- FindClusters(seurat_clean, resolution = 0.6)
seurat_clean <- RunUMAP(seurat_clean, dims = 1:30)
DimPlot(seurat_clean, reduction = "umap", label = TRUE)
```

## 4. Finding Marker Genes with Strict Filtering

```r
# Find markers with stricter thresholds than default
all_markers <- FindAllMarkers(seurat_clean, 
                             only.pos = TRUE, 
                             min.pct = 0.3,            # More stringent min.pct 
                             logfc.threshold = 0.5,    # More stringent logfc
                             test.use = "wilcox",      # Wilcoxon rank sum test
                             min.diff.pct = 0.15)      # Require minimum difference in expression percentage

# Filter markers further to get high-confidence markers
high_confidence_markers <- all_markers %>%
  filter(p_val_adj < 0.01) %>%           # Very stringent p-value cutoff
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

# Save marker files
write.csv(all_markers, "all_markers.csv", row.names = FALSE)
write.csv(high_confidence_markers, "high_confidence_markers.csv", row.names = FALSE)
```

## 5. CASSIA Comprehensive Annotation Workflow

### 5.1 Initial Annotation

```r
# Run CASSIA batch annotation with detailed tissue specification
results <- runCASSIA_batch(
    marker = high_confidence_markers,
    output_name = "high_quality_annotation",
    model = "gpt-4o",
    tissue = "blood",  # Be as specific as possible about your tissue
    species = "human", # Specify correct species
    max_workers = min(8, length(unique(Idents(seurat_clean)))),
    n_genes = 50,      # Use up to 50 top genes per cluster
    additional_info = "These cells are from peripheral blood of healthy adults", # Add relevant context
    provider = "openai"
)
```

### 5.2 Quality Assurance

```r
# Score annotations
quality_scores <- runCASSIA_score_batch(
    input_file = "high_quality_annotation_full.csv",
    output_file = "high_quality_annotation_scored.csv",
    max_workers = 6,
    model = "claude-3-5-sonnet-20241022",
    provider = "anthropic"
)

# Generate detailed HTML report
runCASSIA_generate_score_report(
    csv_path = "high_quality_annotation_scored.csv",
    output_name = "high_quality_annotation_report.html"
)
```

### 5.3 Uncertainty Analysis

```r
# Run multiple annotation iterations to assess uncertainty
iteration_results <- runCASSIA_batch_n_times(
    n = 5,  # 5 independent runs
    marker = high_confidence_markers,
    output_name = "uncertainty_assessment",
    model = "gpt-4o",
    provider = "openai",
    tissue = "blood",
    species = "human",
    max_workers = 5,
    batch_max_workers = 3
)

# Calculate similarity scores across iterations
similarity_scores <- runCASSIA_similarity_score_batch(
    marker = high_confidence_markers,
    file_pattern = "uncertainty_assessment_*_full.csv",
    output_name = "annotation_consistency",
    max_workers = 5,
    model = "claude-3-5-sonnet-20241022",
    provider = "anthropic",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

### 5.4 Detailed Analysis of Challenging Clusters

Use annotation boost for clusters with low confidence scores:

```r
# Identify low-confidence clusters (score < 75)
low_conf_clusters <- read.csv("high_quality_annotation_scored.csv") %>%
  filter(score < 75) %>%
  pull(cluster)

# Run annotation boost on each low-confidence cluster
for(cluster in low_conf_clusters) {
  validation_results <- runCASSIA_annotationboost(
    full_result_path = "high_quality_annotation_full.csv",
    marker = high_confidence_markers,
    output_name = paste0("cluster_", cluster, "_boost"),
    cluster_name = cluster,
    major_cluster_info = "Human peripheral blood mononuclear cells",
    num_iterations = 5,
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter"
  )
}
```

## 6. Integrating CASSIA Results

### 6.1 Adding Annotations to Seurat Object

```r
# Read final annotations
final_annot <- read.csv("high_quality_annotation_scored.csv")

# Create a mapping from cluster to cell type
cluster_celltype_map <- final_annot %>%
  select(cluster, celltype_1, score) %>%
  mutate(confidence = case_when(
    score >= 85 ~ "High",
    score >= 70 ~ "Medium",
    TRUE ~ "Low"
  )) %>%
  mutate(celltype = paste0(celltype_1, " (", confidence, ")"))

# Add annotations to the Seurat object
seurat_clean$celltype <- cluster_celltype_map$celltype[match(Idents(seurat_clean), cluster_celltype_map$cluster)]

# Create a clean version without confidence labels for visualization
seurat_clean$celltype_clean <- cluster_celltype_map$celltype_1[match(Idents(seurat_clean), cluster_celltype_map$cluster)]

# Visualize annotated clusters
DimPlot(seurat_clean, group.by = "celltype_clean", label = TRUE, repel = TRUE) +
  labs(title = "CASSIA Cell Type Annotations") +
  theme(legend.position = "right")
```

### 6.2 Advanced Visualization of Results

```r
# Generate a comprehensive marker heatmap
top_markers_per_celltype <- high_confidence_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Create a mapping between cluster IDs and cell types
cluster_id_map <- final_annot %>%
  select(cluster, celltype_1)

# Add cell type information to markers
top_markers_per_celltype <- top_markers_per_celltype %>%
  left_join(cluster_id_map, by = "cluster") %>%
  arrange(celltype_1, desc(avg_log2FC))

# Create a heatmap with cell type labels
DoHeatmap(seurat_clean, 
          features = top_markers_per_celltype$gene, 
          group.by = "celltype_clean",
          angle = 90) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))

# Generate dotplot of canonical markers
canonical_markers <- c(
  "CD3D", "CD3E", "CD4", "CD8A",  # T cells
  "MS4A1", "CD79A", "CD79B",      # B cells
  "FCGR3A", "CD14", "LYZ",        # Monocytes
  "FCER1A", "CST3",               # Dendritic cells
  "PPBP", "PF4"                   # Platelets
)

DotPlot(seurat_clean, 
        features = canonical_markers, 
        group.by = "celltype_clean") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### 6.3 Subclustering Analysis

For more detailed analysis of specific populations:

```r
# Extract a specific cell population for subclustering
# Example: T cells
t_cells <- subset(seurat_clean, subset = celltype_clean %in% c("CD4+ T cell", "CD8+ T cell", "T cell"))

# Re-process the subset
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells)
t_cells <- ScaleData(t_cells)
t_cells <- RunPCA(t_cells)
t_cells <- FindNeighbors(t_cells, dims = 1:20)
t_cells <- FindClusters(t_cells, resolution = 0.4)
t_cells <- RunUMAP(t_cells, dims = 1:20)

# Find markers for subclusters
t_cell_markers <- FindAllMarkers(t_cells, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)

# Use CASSIA for subclustering analysis
runCASSIA_subclusters(
    marker = t_cell_markers,
    major_cluster_info = "T cells from human PBMC",
    output_name = "t_cell_subclusters",
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter"
)
```

## 7. Final Quality Control and Validation

### 7.1 Cross-validation Against Reference Data

```r
# Use the Retrieve Augmented Agent to access public reference data
retrieval_results <- runCASSIA_retrieve(
    marker = high_confidence_markers,
    output_name = "reference_validation",
    model = "openai/gpt-4o-2024-11-20",
    provider = "openrouter",
    tissue = "blood",
    species = "human"
)
```

### 7.2 Save Final Analysis

```r
# Save the fully annotated Seurat object
saveRDS(seurat_clean, "fully_annotated_seurat_object.rds")

# Create a data package with all important results
dir.create("final_results", showWarnings = FALSE)
file.copy("high_quality_annotation_scored.csv", "final_results/")
file.copy("high_quality_annotation_report.html", "final_results/")
file.copy("annotation_consistency_similarity.csv", "final_results/")

# Create a summary report of the analysis
summary_df <- data.frame(
  cluster = unique(Idents(seurat_clean)),
  cell_count = table(Idents(seurat_clean)),
  cell_type = unique(seurat_clean$celltype_clean),
  confidence_score = final_annot$score
)

write.csv(summary_df, "final_results/cluster_summary.csv", row.names = FALSE)
```

## 8. Troubleshooting and Optimization

### 8.1 Common Issues and Solutions

When working with CASSIA in a comprehensive workflow, consider these troubleshooting steps:

1. **API Rate Limits**: If you encounter rate limit errors:
   ```r
   # Reduce parallel workers and add delays
   max_workers <- 2  # Lower worker count
   Sys.sleep(2)      # Add pauses between API calls
   ```

2. **Memory Issues with Large Datasets**:
   ```r
   # For very large datasets, process in chunks
   # Randomly subset the data for initial analysis
   small_subset <- subset(seurat_obj, cells = sample(colnames(seurat_obj), 5000))
   ```

3. **Cluster Ambiguity Resolution**:
   ```r
   # For ambiguous clusters, use the Compare Cell Types function
   compareCelltypes(
       tissue = "blood",
       celltypes = c("NK cell", "CD8+ T cell", "NKT cell"),
       marker = "your_marker_gene_list_as_string",
       species = "human",
       output_file = "ambiguous_cluster_resolution"
   )
   ```

### 8.2 Best Practices Summary

1. Start with rigorous quality control including adaptive thresholds
2. Apply multiple doublet detection methods
3. Use ambient RNA correction when possible
4. Try multiple clustering resolutions to find optimal granularity
5. Use stricter thresholds for marker detection
6. Provide detailed tissue and context information to CASSIA
7. Validate results with multiple iterations
8. Use annotation boost for challenging clusters
9. Cross-reference with canonical markers and public data
10. Document the entire analysis workflow for reproducibility 