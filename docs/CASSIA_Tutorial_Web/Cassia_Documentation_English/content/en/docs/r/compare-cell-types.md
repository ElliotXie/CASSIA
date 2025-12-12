---
title: Symphony Compare (Optional)
---

## Overview

> **Note:** Symphony Compare requires OpenRouter as the provider. This is the only supported provider for this module.

`symphonyCompare` is an advanced module that orchestrates multiple AI models to compare cell types with automatic consensus building. It conducts a comprehensive cell type comparison using multiple AI models in parallel, automatically triggering discussion rounds when models disagree on the best matching cell type. Think of it as a virtual panel of expert biologists debating and reaching consensus.

## Quick Start

```r
results <- symphonyCompare(
  tissue = "peripheral blood",
  celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
  marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
  species = "human"
)

cat("Consensus:", results$consensus, "\n")
cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\n")
```

## Input

- **Marker genes**: A character vector or comma-separated string of marker genes from CASSIA's previous results or your own analysis
- **Candidate cell types**: A vector of 2-4 cell types to compare
- **Tissue context**: The tissue type being analyzed
- **Species**: The species of the sample (default: human)

## Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `tissue` | character | The tissue type being analyzed (e.g., "blood", "brain", "liver") |
| `celltypes` | character vector | A vector of 2-4 cell types to compare |
| `marker_set` | character vector/string | A list or string of marker genes |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `species` | character | `"human"` | The species of the sample |
| `model_preset` | character | `"premium"` | Pre-configured model ensemble. `"premium"`: High-performance ensemble (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4). `"budget"`: Cost-effective models (DeepSeek V3.2, Grok 4 Fast, Kimi K2, Gemini 2.5 Flash) |
| `enable_discussion` | logical | `TRUE` | If TRUE, models will "discuss" and reconsider their votes if initial consensus is not reached |
| `generate_report` | logical | `TRUE` | Whether to generate an HTML report of the analysis |
| `verbose` | logical | `TRUE` | Whether to print progress messages to the console |

### Advanced Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_discussion_rounds` | integer | `2` | Maximum number of discussion rounds allowed |
| `consensus_threshold` | numeric | `0.8` | The confidence threshold required to declare consensus (0-1) |

## Output

The function returns a list containing consensus results and generates output files.

**Return Value:**
- `results`: List of all model responses and scores
- `consensus`: The consensus cell type (if reached)
- `confidence`: Confidence level of the consensus (0-1)
- `csv_file`: Path to the generated CSV file with detailed results
- `html_file`: Path to the generated interactive HTML report (if enabled)
- `summary`: Summary statistics of the comparison
