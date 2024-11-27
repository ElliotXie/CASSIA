# CASSIA Tutorial
Welcome to the CASSIA R tutorial! This guide will walk you through using the CASSIA package in R for cell type annotation using Large Language Models (LLMs). We'll cover the following key features:

- **Single Cluster Analysis (One Cluster CASSIA)**
- **Batch Processing (Batch CASSIA)**
- **Quality Scoring**
- **Report Generation**
- **Annotation Boost (Validation Plus)**
- **Uncertainty Quantification (Run n times / Similarity Score)**

Let's get started!

---

## Table of Contents

1. [Introduction](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#introduction)
2. [Setup and Installation](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#setup-and-installation)
    - [Installing Required Packages](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#installing-required-packages)
    - [Setting Up the Python Environment](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#setting-up-the-python-environment)
    - [Setting API Keys](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#setting-api-keys)
3. [Single Cluster Analysis](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#single-cluster-analysis)
4. [Batch Processing](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#batch-processing)
5. [Quality Scoring](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#quality-scoring)
6. [Report Generation](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#report-generation)
7. [Annotation Boost (Validation Plus)](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#annotation-boost-validation-plus)
8. [Uncertainty Quantification](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#uncertainty-quantification)
9. [Conclusion](https://chatgpt.com/c/67475af8-2f78-800b-bac0-9593da8d3f6a#conclusion)

---

## Introduction

CASSIA (Collaborative Agent System for Single cell Interpretable Annotation) is a tool designed to enhance cell type annotation by leveraging the power of multi-agent Large Language Models (LLMs).

This tutorial will guide you through the essential functionalities of CASSIA.

---

## Setup and Installation

### Installing Required Packages

First, ensure you have the `reticulate` package installed, which allows R to interface with Python.

```R
install.packages("reticulate")
```

Next, you'll need to install the CASSIA package. You can install it directly:

```R
install.packages("CASSIA")
```

### Setting Up the Python Environment

CASSIA relies on Python for some of its backend processing. When you load the CASSIA package, it attempts to set up the required Python environment automatically. However, if you encounter issues, you can use the `setup_cassia_env()` function to create and configure the necessary Python environment manually.

```R
library(CASSIA)

# Manually set up the Python environment if needed
setup_cassia_env(conda_env = "cassia_env")
```

This function will:

- Create a new Conda environment named `cassia_env` if it doesn't already exist.
- Install the required Python packages: `openai`, `pandas`, `numpy`, `scikit-learn`, `requests`, and `anthropic`.

### Setting API Keys

To use LLMs like OpenAI's GPT-4, Anthropic's Claude, or models via OpenRouter, you need to set your API keys using the `setLLMApiKey()` function.

```R
# For OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# For Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)

# For OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```

- Replace `"your_api_key"` with your actual API key.
- Set `provider` to `"openai"`, `"anthropic"`, or `"openrouter"` depending on your provider.
- Setting `persist = TRUE` saves the key in your `.Renviron` file for future sessions.

---

## Single Cluster Analysis

**One Cluster CASSIA**

Analyze a single cluster of marker genes to determine the cell type.

### Example

If you're using OpenRouter as your provider, you can specify models like `"openai/gpt-4o-2024-11-20"` or `"anthropic/claude-3.5-sonnet"`. Here are some model recommendations:

- **Claude 3.5 Sonnet** (Best performance, slightly more expensive)
    - Model ID: `"anthropic/claude-3.5-sonnet"`
- **GPT-4o** (Balanced option)
    - Model ID: `"openai/gpt-4o-2024-11-20"`
- **Llama 3.2** (Open source, cost-effective)
    - Model ID: `"meta-llama/llama-3.2-90b-vision-instruct"`

#### Example Code

```R
# Parameters
model <- "openai/gpt-4o-2024-11-20"  # Model ID when using OpenRouter
temperature <- 0
marker_list <- c("CD3D", "CD3E", "CD2", "TRAC")
tissue <- "blood"
species <- "human"
additional_info <- NULL
provider <- "openrouter"  # or "openai", "anthropic"

# Run the analysis
result <- runCASSIA(
  model = model,
  temperature = temperature,
  marker_list = marker_list,
  tissue = tissue,
  species = species,
  additional_info = additional_info,
  provider = provider
)

# View structured output
print(result$structured_output)

# View conversation history
print(result$conversation_history)
```

_Note:_ When using OpenRouter, you need to specify the full model ID.

---

## Batch Processing in CASSIA

CASSIA supports batch processing to analyze multiple clusters simultaneously. This guide explains how to prepare your data and run batch analysis efficiently.

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
    model = "gpt-4o",                     # Model to use
    tissue = "brain",                    # Tissue type
    species = "human",                   # Species
    
    # Optional parameters
    max_workers = recommended_workers,    # Number of parallel workers
    n_genes = 50,                        # Number of top marker genes to use
    additional_info = "",                # Additional context
    provider = "openai"                  # API provider
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
   - Document any context-dependent variations in results


---

## Quality Scoring in CASSIA

Quality scoring helps evaluate the reliability of cell type annotations. CASSIA provides automated scoring functionality through the `runCASSIA_score_batch` function, which analyzes the reasoning and evidence behind each annotation.

### Running Quality Scoring

#### Basic Usage
```R
runCASSIA_score_batch(
    input_file = "my_annotation_full.csv",
    output_file = "my_annotation_scored.csv",
    max_workers = 4,
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter"
)
```

#### Parameter Details

1. **Input/Output Files**:
   - `input_file`: Path to the full annotation results (from `runCASSIA_batch`)
   - `output_file`: Where to save the scored results
   
2. **Processing Parameters**:
   - `max_workers`: Number of parallel scoring threads
   - Recommended: Use fewer workers than annotation step to avoid API limits if provider set to anthropic

3. **Model Configuration**:
   - Recommended model: `anthropic/claude-3.5-sonnet`
   - Recommended provider: `openrouter`

### API Provider Considerations

#### OpenRouter
- **Advantages**:
  - Higher rate limits
  - Easy to switch models
- **Setup**:
  ```R
  provider <- "openrouter"
  model <- "anthropic/claude-3.5-sonnet"
  ```

#### Anthropic Direct
- **Considerations**:
  - New users have usage limits
  - May need to reduce `max_workers`
  - Better for smaller datasets
- **Setup**:
  ```R
  provider <- "anthropic"
  model <- "claude-3-5-sonnet-20241022"
  ```

### Output Format
The scored output file contains:
- Original annotation data
- Quality scores (0-100)
- Confidence metrics
- Detailed reasoning for scores

### Interpreting Scores

- **90-100**: High confidence, strong evidence
- **76-89**: Good confidence, adequate evidence
- **<75**: Low confidence, need to run through Annotation Boost Agent and Compare Agent

## Report Generation

Generate detailed reports from your analysis. This step typically follows after quality scoring.

The score report includes all outputs from CASSIA, including structured outputs, conversation histories, and quality scores.

### Batch Reports from Scored Results

```R
runCASSIA_generate_score_report(
  csv_path = "my_annotation_scored.csv",
  output_name = "CASSIA_reports_summary"
)
```

_Generates individual reports and an index page from `scored_results.csv`._

---


## Annotation Boost (ValidationPlus)

Annotation Boost is an advanced validation tool that enhances annotation confidence through multiple iterations of analysis. It's particularly useful for:
- Validating low-confidence annotations
- Getting detailed insights into specific cell clusters
- Resolving ambiguous cell type assignments
- Generating comprehensive validation reports

### Required Components

1. **Input Data**:
   - Full results CSV from CASSIA batch analysis
   - Original marker gene file (Seurat output or custom marker file)
   - Cluster context information
   - Specific cluster identifier

2. **Model Configuration**:
   - Recommended: `anthropic/claude-3.5-sonnet` via OpenRouter
   - Alternative models may be used but may provide different levels of detail

### Running Annotation Boost

```R
# Setup parameters
validation_config <- list(
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter"
)

# Define cluster information
cluster_info <- "Human PBMC"

#Specify the cluster you want to validate
target_cluster = "CD4+ T cell"

# Run validation
runCASSIA_validatorplus(

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
   - Increase `num_iterations`
   - Review marker gene quality

3. **Inconsistent Results**:
   - Check marker gene consistency
   - Verify input data quality
   - Consider biological variability
---


## Uncertainty Quantification

Uncertainty quantification in CASSIA helps assess annotation reliability through multiple analysis iterations and similarity scoring. This process is crucial for:
- Identifying robust cell type assignments
- Detecting mixed or ambiguous clusters
- Quantifying annotation confidence
- Understanding prediction variability

### Multiple Iteration Analysis

#### Basic Usage
```R

# Run multiple analyses
runCASSIA_batch_n_times(
    # Core parameters
    n = 5, #number of iteratioins
    marker = marker_data,
    output_name = "my_annottaion_repeat",
    
    # Model settings
    model = "gpt-4o,
    provider = "openai",
    
    # Context information
    tissue = "brain",
    species = "human",
    additional_info = NULL,


    # Processing control
    max_workers = 4,        # Total parallel workers
    batch_max_workers = 2   # Workers per batch
)
```

#### Parameter Details

1. **Iteration Control**:
   - `n`: Number of analysis iterations
   - Recommended: 5 iterations for standard analysis
   - Consider more iterations for critical applications

2. **Resource Management**:
   - `max_workers`: Overall parallel processing limit
   - `batch_max_workers`: Workers per iteration
   - max_workers * batch_max_workers to match your number of cores.


### Similarity Score Calculation

#### Running Similarity Analysis
```R

# Calculate similarity scores
runCASSIA_similarity_score_batch(
    # Input parameters
    marker = marker_data,
    file_pattern = "my_annottaion_repeat_*_full.csv",
    output_name = "similarity_results",
    
    
    # Processing parameters
    max_workers = 4,
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter",
    
    # Scoring weights
    main_weight = 0.5, # Weight for main cell type
    sub_weight = 0.5  # Weight for subtypes
)
```

#### Scoring Parameters

1. **Weight Configuration**:
   - `main_weight`: Importance of main cell type match (0-1)
   - `sub_weight`: Importance of subtype match (0-1)
   - Weights should sum to 1.0

2. **File Management**:
   - `file_pattern`: Pattern to match iteration results
   - Uses * to match iteration numbers
   - Example:  if you have "my_annottaion_repeat_1_full.csv", "my_annottaion_repeat_2_full.csv", and "my_annottaion_repeat_3_full.csv", use "my_annottaion_repeat__full.csv" to match the pattern.

#### Output Interpretation

1. **Similarity Scores**:
   - Range: 0 (completely different) to 1 (identical)
   - Interpretation guidelines:
     - >0.9: High consistency
     - 0.75-0.9: Moderate consistency
     - <0.75: Low consistency

#### Troubleshooting

1. **Performance Issues**:
   - Reduce worker counts
   - Process in smaller batches

2. **Low Similarity Scores**:
   - Review marker gene quality
   - Use Annotation Boost function
   - Review cluster heterogeneity
   - Consider biological variability
   - Increase iteration count
   - Try subclustering

---

## Conclusion

This tutorial has covered the core functionalities of CASSIA, from basic cell type analysis to batch processing, quality evaluation, and uncertainty quantification. You've learned to leverage both basic and advanced features to generate reliable cell type annotations for your single-cell data.

CASSIA is actively evolving, and we're excited to expand its capabilities. Future developments will include integration with clustering, multi-modal data analysis, and enhanced visualization tools.

Stay tuned for updates as we continue enhancing CASSIA's capabilities! Follow our GitHub repository for the latest developments.

Happy annotating!

---

## Additional Resources
- [Refer to the paper for more detail and example]()
- **API Provider Guides:**
	- [How to get an OpenAI api key](https://platform.openai.com/api-keys)
	- [How to get an Anthropic api key](https://console.anthropic.com/settings/keys)
	- [How to get an OpenRouter api key](https://openrouter.ai/settings/keys)
    - [OpenAI API Documentation](https://beta.openai.com/docs/)
    - [Anthropic API Documentation](https://docs.anthropic.com/)
    - [OpenRouter API documentatioon](https://openrouter.ai/docs/quick-start)

---
