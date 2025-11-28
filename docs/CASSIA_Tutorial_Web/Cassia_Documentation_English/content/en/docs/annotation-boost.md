---
title: Annotation Boost (Optional)
---


Annotation Boost is an advanced validation tool that enhances annotation confidence through multiple iterations of analysis. It's particularly useful for:
- Validating low-confidence annotations
- Getting detailed insights into specific cell clusters
- Resolving ambiguous cell type assignments
- Generating comprehensive validation reports

### Required Components

1. **Input Data**:
   - Full results CSV from CASSIA batch analysis
   - Original marker gene file (***Note: The marker file should not be filtered!***)
   - Cluster context information
   - Specific cluster identifier

2. **Model Configuration**:
   - Recommended: `google/gemini-2.5-flash-preview` via OpenRouter

### Running Annotation Boost

```R
# Setup parameters
validation_config <- list(
    model = "google/gemini-2.5-flash-preview",
    provider = "openrouter"
)

# Define cluster information
cluster_info <- "Human PBMC"

#Specify the cluster you want to validate
target_cluster = "CD4+ T cell"

# Run validation
runCASSIA_annotationboost(

    # Required parameters
    full_result_path = "cell_type_analysis_results.csv",
    marker = marker_data,
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    output_name = "Cluster1_report",
    num_iterations = 5, # Number of validation rounds

    # Model configuration
    model = validation_config$model,
    provider = validation_config$provider,
)
```

### Parameter Details

   - `full_result_path`: Path to original CASSIA results
   - `marker`: Marker gene data (same as used in initial analysis)
   - `cluster_name`: Target cluster name
   - `major_cluster_info`: Dataset context
   - `num_iterations`: Number of validation rounds (default: 5)


### Troubleshooting

1. **Low Confidence Results**:
   - Review the quality of the clusters, in terms of doublet, mixed, or low quality clusters.
   - Review marker gene quality

3. **Inconsistent Results**:
   - Check marker gene consistency
   - Verify input data quality
