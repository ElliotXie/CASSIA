---
title: Batch Processing
---

CASSIA supports batch processing to analyze multiple clusters simultaneously. This guide explains how to prepare your data and run batch analysis efficiently.


### Model Recommendations
If you're using OpenRouter as your provider, you can specify models like `"openai/gpt-4o-2024-11-20"` or `"anthropic/claude-3.5-sonnet"`. Here are some model recommendations:

- **Gemini 2.5 Flash** (Best performance, cost-effective)
    - Model ID: `"google/gemini-2.5-flash-preview"`
- **Deepseek v3** (Open source, excellent performance, cost-effective)
    - Model ID: `"deepseek/deepseek-chat-v3-0324"`
- **Deepseek v3 Free** (Free tier, slower)
    - Model ID: `"deepseek/deepseek-chat-v3-0324:free"`

### Preparing Marker Data
You have three options for providing marker data:

1. Create a data frame or CSV file containing your clusters and marker genes
2. Use Seurat's `findAllMarkers` function output directly
3. Use CASSIA's example marker data

```R
# Option 1: Load your own marker data
markers <- read.csv("path/to/your/markers.csv")

# Option 2: Use Seurat's findAllMarkers output directly
# (assuming you already have a Seurat object)
markers <- FindAllMarkers(seurat_obj)

# Option 3: Load example marker data
markers <- loadExampleMarkers()

# Preview the data
head(markers)
```

#### Marker Data Format
CASSIA accepts two formats:

1. **FindAllMarkers Output**: The standard output from Seurat's FindAllMarkers function
2. **Simplified Format**: A two-column data frame where:
   - First column: cluster identifier
   - Second column: comma-separated ranked marker genes

### Running Batch Analysis

#### Setting Up Parameters

```R
# Detect available CPU cores
available_cores <- parallel::detectCores()

# Calculate recommended workers (75% of available cores)
recommended_workers <- floor(available_cores * 0.75)

runCASSIA_batch(
    # Required parameters
    marker = markers,                    # Marker data (data frame or file path)
    output_name = "my_annotation",       # Base name for output files
    model = "google/gemini-2.5-flash-preview", # Model to use
    tissue = "brain",                    # Tissue type
    species = "human",                   # Species
    
    # Optional parameters
    max_workers = recommended_workers,    # Number of parallel workers
    n_genes = 50,                        # Number of top marker genes to use
    additional_info = "",                # Additional context
    provider = "openrouter"              # API provider
)
```

### Parameter Details

1. **Marker Gene Selection**:
   - Default: top 50 genes per cluster
   - Filtering criteria:
     - Adjusted p-value < 0.05
     - Average log2 fold change > 0.25
     - Minimum percentage > 0.1
   - If fewer than 50 genes pass filters, all passing genes are used

2. **Parallel Processing**:
   - `max_workers`: Controls parallel processing threads
   - Recommended: 80% of available CPU cores
   - Example: For a 16-core machine, set to 13

3. **Additional Context** (optional):
   - Use `additional_info` to provide experimental context
   - Examples:
     - Treatment conditions: "Samples were antibody-treated"
     - Analysis focus: "Please carefully distinguish between cancer and non-cancer cells"
   - Tip: Compare results with and without additional context

### Output Files

The analysis generates two files:
1. `my_annotation_full.csv`: Complete conversation history
2. `my_annotation_summary.csv`: Condensed results summary

### Tips for Optimal Results

1. **Resource Management**:
   - Monitor system resources when setting `max_workers`
   - Start with recommended 75% of cores and adjust if needed

2. **Marker Gene Selection**:
   - Default 50 genes works well for most cases
   - Increase for more complex cell types
   - Decrease if running into API rate limits

3. **Context Optimization**:
   - Test runs with and without additional context
   - Keep context concise and relevant
