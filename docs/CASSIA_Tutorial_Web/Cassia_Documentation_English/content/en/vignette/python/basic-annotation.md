---
title: Basic Annotation
---

This Python tutorial demonstrates a complete workflow using CASSIA for cell type annotation of single-cell RNA sequencing data. We'll analyze an intestinal cell dataset containing six distinct populations:

1. monocyte
2. plasma cells
3. cd8-positive, alpha-beta t cell
4. transit amplifying cell of large intestine
5. intestinal enteroendocrine cell
6. intestinal crypt stem cell

## Setup and Environment Preparation

First, let's install and import the required packages:

```bash
pip install CASSIA
```

```python
import CASSIA
```

### Set API keys

**You only need to choose one provider.** OpenRouter is recommended as it provides access to multiple models.

```python
# Set API key (choose one provider)
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")  # Recommended
# CASSIA.set_api_key("your-openai-key", provider="openai")
# CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
```

### Load Data

```python
processed_markers = CASSIA.loadmarker(marker_type="processed")
unprocessed_markers = CASSIA.loadmarker(marker_type="unprocessed")
subcluster_results = CASSIA.loadmarker(marker_type="subcluster_results")

# List available marker sets
available_markers = CASSIA.list_available_markers()
print(available_markers) 
```

## Fast Mode

Run the CASSIA pipeline in fast mode. This is a one-step process to get results quickly.

```python
# Run the CASSIA pipeline in fast mode
CASSIA.runCASSIA_pipeline(
    output_file_name = "FastAnalysisResults",
    tissue = "large intestine",
    species = "human",
    marker = unprocessed_markers,
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

## Detailed Batch Analysis

For more granular control, you can run the steps individually.

```python
output_name="intestine_detailed"

# Run batch analysis
CASSIA.runCASSIA_batch(
    marker = unprocessed_markers,
    output_name = output_name,
    model = "anthropic/claude-sonnet-4.5",
    tissue = "large intestine",
    species = "human",
    max_workers = 6,  # Matching cluster count
    n_genes = 50,
    additional_info = None,
    provider = "openrouter")
```

## Quality Scoring

After annotation, run quality scoring to assess confidence.

```python
# Run quality scoring
CASSIA.runCASSIA_score_batch(
    input_file = output_name + "_full.csv",
    output_file = output_name + "_scored.csv",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter"
)

# Generate quality report
CASSIA.runCASSIA_generate_score_report(
    csv_path = output_name + "_scored.csv",
    index_name = output_name + "_report.html"
)
```

