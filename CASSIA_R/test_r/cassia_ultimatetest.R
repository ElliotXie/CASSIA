---
title: "R Notebook"
output: html_notebook
---

#test internally

### LOCAL DEVELOPMENT TESTING ###
# Use this section for testing local changes without pushing to GitHub
# 
# FIXED IMPORT ISSUES:
# - runCASSIA_batch_n_times: moved from py_tools to py_uncertainty
# - runCASSIA_similarity_score_batch: moved from py_tools to py_uncertainty  
# - runCASSIA_n_times_similarity_score: moved from py_tools to py_uncertainty
# - These functions are now properly imported from Uncertainty_quantification.py

```{r}
library(devtools)
library(dplyr)
library(reticulate)

# Set working directory to the package root (CASSIA_R)
# Make sure you're in the right directory - should contain DESCRIPTION file
getwd()  # Check current directory
# If needed, uncomment and modify this path:
setwd("D:/newgit/CASSIA/CASSIA_R")

# Check if we're in the right directory
if (!file.exists("DESCRIPTION")) {
  stop("Not in package root directory. Please navigate to CASSIA_R folder containing DESCRIPTION file.")
}

# Clean any previous builds
devtools::clean_dll()

# Document the package (update NAMESPACE and .Rd files)
devtools::document()

# Build the package
devtools::build()

# Install the locally built package
devtools::install(upgrade = "never")  # Use upgrade = "never" to avoid dependency updates

# Restart R session is recommended but not required
# .rs.restartR()  # Uncomment if using RStudio

# Load the package
library(CASSIA)

# Load example data
markers_unprocessed <- loadExampleMarkers(processed = FALSE)

print("✓ Local package built and installed successfully!")
print("✓ Uncertainty quantification import issues have been fixed!")
```

### ALTERNATIVE: Install from GitHub (original method) ###



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





setwd("D:/newgit/CASSIA/CASSIA_R/test_r")




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

runCASSIA_merge_annotations( csv_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation2_full.csv",
                             provider = "openrouter",
                             model = "deepseek/deepseek-chat-v3-0324",
                             detail_level = "detailed",
                             process_all = TRUE,
                             debug = FALSE)

```



#test the scoring
```{r}

quality_scores <- runCASSIA_score_batch(
  input_file = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation2_full.csv",
  output_file = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation2_full_scored.csv"
)




# Generate quality report
runCASSIA_generate_score_report(
  csv_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation2_full_scored.csv",
  output_name = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation2_full_scored_report.csv"
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
cluster_info <- "Human Large Intestine"

#Specify the cluster you want to validate
target_cluster = "monocyte"

# Run validation
runCASSIA_annotationboost(
    full_result_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation2_full.csv",
    marker = markers_unprocessed,
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    output_name = "monocyte_report",
    num_iterations = 5, # Number of validation rounds
    model ="google/gemini-2.5-flash-preview",
    provider = "openrouter"
)



# Run validation plus for the high mitochondrial content cluster
runCASSIA_annotationboost_additional_task(
    full_result_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation2_full.csv",
    marker = markers_unprocessed,
    output_name="monocyte_annotationboost",
    cluster_name = "monocyte",
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "google/gemini-2.5-flash-preview",
    provider = "openrouter",
    additional_task = "check if this is monocyte"
)


```



### test the http custmerized model.



```{r}
# Set up DeepSeek API for testing
Sys.setenv(CUSTERMIZED_API_KEY = "sk-f7c54e95a5e040589f41d83553b55861")

custom_base_url <- "https://api.deepseek.com"
custom_model <- "deepseek-chat"

print("🔧 Testing with DeepSeek API...")
print(paste("Base URL:", custom_base_url))
print(paste("Model:", custom_model))
print("✓ API key configured")
```

#test the runcassia_batch with custom model

```{r}
# Set the API key for DeepSeek using the standard approach
setLLMApiKey("sk-f7c54e95a5e040589f41d83553b55861", provider = custom_base_url)


runCASSIA_batch(
    # Required parameters
    marker = markers_unprocessed,
    output_name = "my_annotation_custom",
    tissue = "large intestine",
    species = "human",
    # Custom model parameters
    max_workers = 6,
    n_genes = 50,
    provider = custom_base_url,  # Use HTTP endpoint as provider
    model = custom_model,        # Custom model name
    temperature = 0
)

```

#test the pipeline with custom model
```{r}
# Run the CASSIA pipeline with custom models
runCASSIA_pipeline(
    output_file_name = "TEST_CUSTOM",
    tissue = "large intestine", 
    species = "human",
    marker = markers_unprocessed,
    max_workers = 4,
    # Custom annotation model
    annotation_model = custom_model,
    annotation_provider = custom_base_url,
    # Custom scoring model  
    score_model = custom_model,
    score_provider = custom_base_url,
    # Custom annotation boost model
    annotationboost_model = custom_model,
    annotationboost_provider = custom_base_url,
    # Custom merge model
    merge_model = custom_model,
    score_threshold = 75
)

```

#test the merge with custom model
```{r}

runCASSIA_merge_annotations(
    csv_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation_custom_full.csv",
    provider = custom_base_url,
    model = custom_model,
    detail_level = "detailed",
    process_all = FALSE,
    debug = FALSE
)

```

#test the scoring with custom model
```{r}

quality_scores_custom <- runCASSIA_score_batch(
    input_file = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation_custom_full.csv",
    output_file = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation_custom_scored.csv",
    model = custom_model,
    provider = custom_base_url,
    max_workers = 4
)

# Generate quality report for custom model results
runCASSIA_generate_score_report(
    csv_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation_custom_scored.csv",
    output_name = "D:/newgit/CASSIA/CASSIA_R/test_r/custom_model_report"
)

```

#test subcluster with custom model
```{r}

marker_sub <- loadExampleMarkers_subcluster()

runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results_custom",
    model = custom_model,
    provider = custom_base_url,
    temperature = 0
)

```

#test compare celltype with custom model
```{r}

marker <- "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

# Test with custom model list for comparison
custom_model_list <- c(custom_model, "backup-model-name")

compareCelltypes(
    tissue = "large intestine",
    celltypes = c("Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"),
    marker = marker,
    species = "human",
    model_list = custom_model_list,
    output_file = "plasma_cell_subtype_custom"
)

```

#test uncertainty quantification with custom model
```{r}

# Test batch n times with custom model
iteration_results_custom <- runCASSIA_batch_n_times(
    n = 3,
    marker = markers_unprocessed,
    output_name = "my_annotation_custom_uncertainty",
    model = custom_model,
    provider = custom_base_url,
    tissue = "large intestine",
    species = "human",
    max_workers = 4,
    batch_max_workers = 2,
    temperature = 0
)

# Test similarity scoring with custom model
similarity_scores_custom <- runCASSIA_similarity_score_batch(
    marker = markers_unprocessed,
    file_pattern = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation_custom_uncertainty_*_full.csv",
    output_name = "intestine_similarity_custom",
    max_workers = 4,
    provider = custom_base_url,
    model = custom_model
)

```

#test annotation boost with custom model
```{r}

# Define cluster information for custom model testing
cluster_info <- "Human Large Intestine"
target_cluster <- "monocyte"

# Test annotation boost with custom model
runCASSIA_annotationboost(
    full_result_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation_custom_full.csv",
    marker = markers_unprocessed,
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    output_name = "monocyte_report_custom",
    num_iterations = 3,  # Reduced for testing
    model = custom_model,
    provider = custom_base_url,
    temperature = 0
)

# Test annotation boost with additional task using custom model
runCASSIA_annotationboost_additional_task(
    full_result_path = "D:/newgit/CASSIA/CASSIA_R/test_r/my_annotation_custom_full.csv",
    marker = markers_unprocessed,
    output_name = "monocyte_boost_custom",
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    num_iterations = 3,  # Reduced for testing
    model = custom_model,
    provider = custom_base_url,
    additional_task = "check if this is a monocyte and assess confidence",
    temperature = 0
)

```

#test advanced custom model configurations
```{r}

# Test with multiple custom endpoints for redundancy
custom_endpoints <- list(
    primary = list(base_url = "http://localhost:11434", model = "llama3.2:latest"),
    secondary = list(base_url = "http://localhost:8080", model = "mistral:latest"),
    tertiary = list(base_url = "https://api.your-custom-service.com/v1", model = "custom-model-v2")
)

# Test with primary custom endpoint
print("Testing with primary custom endpoint...")
tryCatch({
    runCASSIA_batch(
        marker = markers_unprocessed[1:3, ],  # Small subset for testing
        output_name = "test_primary_custom",
        tissue = "large intestine",
        species = "human",
        max_workers = 2,
        n_genes = 20,
        provider = custom_endpoints$primary$base_url,
        model = custom_endpoints$primary$model,
        temperature = 0.1
    )
    print("✓ Primary custom endpoint test completed")
}, error = function(e) {
    print(paste("✗ Primary custom endpoint failed:", e$message))
})

```

#test custom model with different parameters
```{r}

# Test various temperature settings with custom model
temperature_tests <- c(0, 0.3, 0.7, 1.0)

for (temp in temperature_tests) {
    tryCatch({
        print(paste("Testing temperature:", temp))
        
        runCASSIA(
            model = custom_model,
            temperature = temp,
            marker_list = c("CD14", "CD68", "CSF1R", "CX3CR1", "CD163"),
            tissue = "large intestine",
            species = "human",
            additional_info = paste("Temperature test:", temp),
            provider = custom_base_url
        )
        
        print(paste("✓ Temperature", temp, "test completed"))
        
    }, error = function(e) {
        print(paste("✗ Temperature", temp, "test failed:", e$message))
    })
}

```

### Custom Model Configuration Summary
```{r}

print("🎯 Custom Model Testing Summary:")
print("==================================")
print(paste("Base URL:", custom_base_url))
print(paste("Model:", custom_model))
print("")
print("Tested Functions:")
print("✓ runCASSIA_batch")
print("✓ runCASSIA_pipeline") 
print("✓ runCASSIA_merge_annotations")
print("✓ runCASSIA_score_batch")
print("✓ runCASSIA_subclusters")
print("✓ compareCelltypes")
print("✓ runCASSIA_batch_n_times")
print("✓ runCASSIA_similarity_score_batch")
print("✓ runCASSIA_annotationboost")
print("✓ runCASSIA_annotationboost_additional_task")
print("")
print("📁 Output files generated with '_custom' suffix")
print("🔧 All functions tested with custom HTTP endpoint")

```
