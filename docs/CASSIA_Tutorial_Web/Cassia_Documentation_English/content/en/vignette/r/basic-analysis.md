---
title: "Getting Started with CASSIA: Basic Analysis"
---

This vignette covers the fundamental steps to get started with CASSIA for single-cell RNA sequencing analysis. We'll focus on the essential workflow for basic cell type annotation.

## 1. Installation and Setup

Before using CASSIA, you need to install it and set up the necessary environment.

### 1.1 Package Installation

```r
# Install prerequisite packages
install.packages("reticulate")
install.packages("devtools")

# Install CASSIA from GitHub
library(devtools)
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")

# Load the CASSIA package
library(CASSIA)
```

### 1.2 Python Environment Setup

CASSIA relies on Python for some backend operations. The environment will be set up automatically when you load the package, but if you encounter issues:

```r
# Manually set up the Python environment
setup_cassia_env(conda_env = "cassia_env")
```

### 1.3 API Key Configuration

CASSIA requires at least one API key to function. You can use OpenAI, Anthropic, or OpenRouter:

```r
# Set up at least one API key (choose the provider you have an account with)
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)
# Or
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)
# Or
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```

Setting `persist = TRUE` saves the key to your `.Renviron` file for future sessions.

## 2. Preparing Your Data

CASSIA works with marker gene data, typically produced by the `FindAllMarkers` function in Seurat.

### 2.1 Sample Data

For this vignette, we'll use the example data included with CASSIA:

```r
# Load example marker data (output from Seurat's FindAllMarkers)
markers <- loadExampleMarkers(processed = FALSE)

# Preview the data
head(markers)
```

### 2.2 Using Your Own Data

If you have your own Seurat object, here's how to prepare the data:

```r
# Example of preparing your own data from a Seurat object
# seurat_obj <- your_seurat_object
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)
# seurat_obj <- FindNeighbors(seurat_obj)
# seurat_obj <- FindClusters(seurat_obj)
# markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

## 3. Basic Analysis Workflow

### 3.1 Fast Analysis (All-In-One)

For a quick analysis that performs all steps (annotation, scoring, and annotation boost):

```r
# Run the complete CASSIA pipeline in fast mode
results <- runCASSIA_pipeline(
    output_file_name = "MyAnalysis",
    tissue = "blood", # Specify the tissue type of your sample
    species = "human", # Specify "human" or "mouse"
    marker = markers,
    max_workers = 4, # Adjust based on your computer's capabilities
    annotation_model = "anthropic/claude-sonnet-4.5", # Model for initial annotation
    annotation_provider = "openrouter",
    score_model = "anthropic/claude-sonnet-4.5", # Model for quality scoring
    score_provider = "openrouter",
    score_threshold = 75, # Minimum score for high-quality annotations
    annotationboost_model = "anthropic/claude-sonnet-4.5", # Model for annotation boost
    annotationboost_provider = "openrouter",
    do_merge_annotations = TRUE, # Enable merging of annotations
    validator_involvement = "v1" # Validator involvement level (v0=high, v1=moderate)
)
```

### 3.2 Step-by-Step Analysis

If you prefer more control, you can run each step separately:

```r
# Step 1: Initial annotation
batch_results <- runCASSIA_batch(
    marker = markers,
    output_name = "StepByStep",
    model = "anthropic/claude-sonnet-4.5",
    tissue = "blood",
    species = "human",
    max_workers = 4,
    n_genes = 50, # Number of top marker genes to use
    provider = "openrouter",
    validator_involvement = "v1"
)

# Step 2: Score the annotations
quality_scores <- runCASSIA_score_batch(
    input_file = "StepByStep_full.csv",
    output_file = "StepByStep_scored.csv",
    max_workers = 4,
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)

# Step 3: Generate a report
runCASSIA_generate_score_report(
    csv_path = "StepByStep_scored.csv",
    output_name = "StepByStep_report.html"
)
```

## 4. Interpreting Results

CASSIA produces several output files:

### 4.1 Full Results CSV

The `_full.csv` file contains detailed annotation information for each cluster:

```r
# Read the full results file
results <- read.csv("MyAnalysis_full.csv")
head(results)
```

Key columns include:
- `cluster`: The cluster identifier
- `celltype_1`, `celltype_2`, `celltype_3`: Top three predicted cell types
- `reasoning`: The model's explanation for the annotation
- `confidence`: The model's confidence in the annotation

### 4.2 Scored Results CSV

The `_scored.csv` file includes quality scores for each annotation:

```r
# Read the scored results file
scored <- read.csv("MyAnalysis_scored.csv")
head(scored)
```

Key additional columns:
- `score`: Quality score from 0-100
- `score_reasoning`: Explanation of the score
- `score_category`: Category based on score threshold

### 4.3 HTML Report

The HTML report provides a visual summary of the annotations, including:
- Distribution of scores across clusters
- Detailed annotation information for each cluster
- Quality assessment and recommendations

## 5. Common Issues and Solutions

### 5.1 API Rate Limits

If you encounter rate limit errors:

```r
# Adjust the maximum number of workers to avoid rate limits
results <- runCASSIA_pipeline(
    # ... other parameters ...
    max_workers = 2, # Reduce to avoid rate limits
    # ... other parameters ...
)
```

### 5.2 Memory Issues

For large datasets:

```r
# Process fewer genes per cluster
results <- runCASSIA_pipeline(
    # ... other parameters ...
    n_genes = 30, # Reduce from default 50
    # ... other parameters ...
)
```

## 6. Next Steps

After completing basic analysis, consider:

- Use `runCASSIA_annotationboost()` for clusters with low confidence scores
- Try `compareCelltypes()` to resolve ambiguous annotations
- Implement `runCASSIA_batch_n_times()` for uncertainty quantification
- Explore `runCASSIA_subclusters()` for detailed analysis of specific cell populations

Refer to the "Complete CASSIA Workflow" vignette for details on these advanced techniques. 