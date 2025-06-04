---
title: "R Notebook"
output: html_notebook
---

#install the gem and load the data

```{r}
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
devtools::install_github("ElliotXie/CASSIA/CASSIA_R", ref = "dev_refactor")

library(CASSIA)
library(dplyr)
library(reticulate)


markers_unprocessed <- loadExampleMarkers(processed = FALSE)

# Get unique cluster names in their original order
cluster_levels <- unique(markers_unprocessed$cluster)
# Re-factor the cluster column using these levels
markers_unprocessed$cluster <- factor(markers_unprocessed$cluster, levels = cluster_levels, labels = seq_along(cluster_levels))


```



#update the code

```{r}

###只需要copy就可以了！！！！


devtools::install_github("ElliotXie/CASSIA/CASSIA_R")


packagePath <- find.package("CASSIA")

# Copy all Python modules to the installed package
pythonScriptPath <- file.path(packagePath, "python", "tools_function.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/tools_function.py", pythonScriptPath, overwrite = TRUE)

pythonScriptPath <- file.path(packagePath, "python", "main_function_code.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/main_function_code.py", pythonScriptPath, overwrite = TRUE)

pythonScriptPath <- file.path(packagePath, "python", "annotation_boost.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/annotation_boost.py", pythonScriptPath, overwrite = TRUE)

# Copy additional Python modules
pythonScriptPath <- file.path(packagePath, "python", "subclustering.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/subclustering.py", pythonScriptPath, overwrite = TRUE)

pythonScriptPath <- file.path(packagePath, "python", "Uncertainty_quantification.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/Uncertainty_quantification.py", pythonScriptPath, overwrite = TRUE)

pythonScriptPath <- file.path(packagePath, "python", "merging_annotation_code.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/merging_annotation_code.py", pythonScriptPath, overwrite = TRUE)

pythonScriptPath <- file.path(packagePath, "python", "cell_type_comparison.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/cell_type_comparison.py", pythonScriptPath, overwrite = TRUE)

pythonScriptPath <- file.path(packagePath, "python", "llm_utils.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/llm_utils.py", pythonScriptPath, overwrite = TRUE)

pythonScriptPath <- file.path(packagePath, "python", "generate_comparison.py")
file.copy("D:/newgit/CASSIA/CASSIA_R/inst/python/generate_comparison.py", pythonScriptPath, overwrite = TRUE)

detach("package:CASSIA", unload = TRUE)

library(CASSIA)

```




#test the runpipeline
```{r}
# Run the CASSIA pipeline in fast mode
fast_results <- runCASSIA_pipeline(
    output_file_name = "TEST2",
    tissue = "large intestine",
    species = "human",
    marker = markers_unprocessed,
    max_workers = 6,
    score_threshold = 92
)






fast_results <- py_tools$run_cell_analysis_pipeline(
    output_file_name = "TEST1",
    tissue = "large intestine",
    species = "human",
    marker_path = markers_unprocessed,
    score_model = "deepseek/deepseek-chat-v3-0324",
    score_provider = "openrouter",
    max_workers = 6,
    score_threshold = 95
)


```


#test the runcassia_batch

```{r}

runCASSIA_batch(
    # Required parameters
    marker = markers_unprocessed,                    # Marker data (data frame or file path)
    output_name = "my_annotation2",       # Base name for output files
    tissue = "large intestine",                    # Tissue type
    species = "human",                   # Species
    # Optional parameters
    max_workers = 6,    # Number of parallel workers
    n_genes = 50,                        # Number of top marker genes to use
    provider = "openrouter",
    model="google/gemini-2.5-flash-preview",
)

```


#test the merge
```{r}

```



#test the scoring
```{r}


quality_scores <- runCASSIA_score_batch(
  input_file = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Uncertainty_gemini_4_full.csv",
  output_file = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Uncertainty_gemini_4_full_scored.csv"
)




# Generate quality report
runCASSIA_generate_score_report(
  csv_path = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Uncertainty_gemini_4_full_scored.csv",
  output_name = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Uncertainty_gemini_4_full_scored_report.csv")
)

```

#test subcluster
```{r}

marker_sub=loadExampleMarkers_subcluster()


runCASSIA_subclusters(marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "google/gemini-2.5-flash-preview",
    provider = "openrouter")


```


#test compare celltype

```{r}

marker="IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"


compareCelltypes(
      tissue = "large intestine",
      celltypes = c("Plasma Cells","IgA-secreting Plasma Cells","IgG-secreting Plasma Cells","IgM-secreting Plasma Cells"),
      marker = marker,
      species = "human",
      output_file = "plasama_cell_subtype"
    )
  
  
  
```

#test the cs score

```{r}

iteration_results <- runCASSIA_batch_n_times(
    n = 3,
    marker = markers_unprocessed,
    output_name = "my_annotation",
    model = "meta-llama/llama-4-maverick",
    provider = "openrouter",
    tissue = "large intestine",
    species = "human",
    max_workers = 6,
    batch_max_workers = 3
)


similarity_scores <- runCASSIA_similarity_score_batch(
    marker = markers_unprocessed,
    file_pattern = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Uncertainty_gemini_*_full.csv", # The file pattern of the uncertainty results
    output_name = "intestine_similarity_test6",
    max_workers = 6,
    provider = "openrouter",
    model = "google/gemini-2.5-flash-preview",
)




similarity_scores <- runCASSIA_similarity_score_batch(
    marker = markers_unprocessed,
    file_pattern = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Uncertainty_gemini_*_full.csv", # The file pattern of the uncertainty results
    output_name = "intestine_similarity_test3",
    max_workers = 6,
    provider = "openai",
    model = "gpt-4o-mini",
)




```


#test the annotation boost
```{r}
library(CASSIA)


# Define cluster information
cluster_info <- "Human Breast"

#Specify the cluster you want to validate
target_cluster = "11"

# Run validation
runCASSIA_annotationboost(
    full_result_path = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Breast_Human_20250504_221205/cassia_test3_full.csv",
    marker = marker_cluster_filter2,
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    output_name = "Cluster1_report",
    num_iterations = 5, # Number of validation rounds
    model ="google/gemini-2.5-flash-preview",
    provider = "openrouter"
)



# Run validation plus for the high mitochondrial content cluster
runCASSIA_annotationboost(
    full_result_path = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/CASSIA_Uncertainty_gemini_4_full.csv",
    marker = markers_unprocessed,
    output_name="monocyte_annotationboost2",
    cluster_name = "monocyte",
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "google/gemini-2.5-flash-preview",
    provider = "openrouter"
)




```




#Actual seurat object test

```{r}
gtex_data=readRDS("processed_breast.RDS")

gtex_data@meta.data %>% colnames()

# If starting with raw data, perform QC
# Calculate mitochondrial percentage
gtex_data$mitoRatio <- PercentageFeatureSet(object = gtex_data, pattern = "^MT-")


genes=rownames(gtex_data) %>% as.data.frame()

gtex_data$mitoRatio %>% summary



# Filter cells based on QC metrics
gtex_data <- subset(gtex_data, subset = nFeature_RNA > 200 & nCount_RNA >= 500 & percent.mt < 10)





```






#ts full test quality control included.

```{r}
seurat_obj=readRDS("C:/Users/ellio/Downloads/seurat_object (10).rds")
```



```{r}

# Calculate mitochondrial and ribosomal percentages
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")



seurat_obj$percent.mt %>% summary

qc_violin


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














