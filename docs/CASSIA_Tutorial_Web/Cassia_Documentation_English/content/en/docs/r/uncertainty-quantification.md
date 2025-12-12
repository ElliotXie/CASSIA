---
title: Uncertainty Quantification (Optional)
---

## Overview

Uncertainty quantification in CASSIA helps assess annotation reliability through multiple analysis iterations and similarity scoring. This process is crucial for:
- Identifying robust cell type assignments
- Detecting mixed or ambiguous clusters
- Quantifying annotation confidence
- Understanding prediction variability

> **Cost Warning**: Running multiple iterations with LLM models can incur significant costs. Each iteration makes separate API calls, so the total cost will be approximately n times the cost of a single run.

## Quick Start

```R
library(CASSIA)

# Step 1: Run multiple iterations
runCASSIA_batch_n_times(
    n = 5,
    marker = marker_data,
    output_name = "my_annotation",
    model = "openai/gpt-5.1",
    provider = "openrouter",
    tissue = "brain",
    species = "human",
    reasoning = "medium"
)

# Step 2: Calculate similarity scores
runCASSIA_similarity_score_batch(
    marker = marker_data,
    file_pattern = "my_annotation_*_full.csv",
    output_name = "similarity_results",
    model = "openai/gpt-5.1",
    provider = "openrouter",
    reasoning = "medium"
)
```

## Input

| Input | Description | Format |
|-------|-------------|--------|
| `marker` | Marker gene data | Data frame or file path |
| `tissue` | Tissue type context | String (e.g., "brain", "large intestine") |
| `species` | Species context | String (e.g., "human", "mouse") |
| `file_pattern` | Pattern to match iteration results | Glob pattern with `*` wildcard |

## Parameters

### Batch Iteration (`runCASSIA_batch_n_times`)

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `n` | Yes | - | Number of analysis iterations (recommended: 5) |
| `marker` | Yes | - | Marker gene data (data frame or path) |
| `output_name` | Yes | - | Base name for output files |
| `model` | Yes | - | LLM model to use |
| `provider` | Yes | - | API provider |
| `tissue` | Yes | - | Tissue type |
| `species` | Yes | - | Species |
| `max_workers` | No | 4 | Overall parallel processing limit |
| `batch_max_workers` | No | 2 | Workers per iteration (max_workers * batch_max_workers should match your cores) |
| `reasoning` | No | NULL | Reasoning effort level ("low", "medium", "high") - only for GPT-5 models |

### Similarity Scoring (`runCASSIA_similarity_score_batch`)

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `marker` | Yes | - | Marker gene data |
| `file_pattern` | Yes | - | Pattern to match iteration results (e.g., `"output_*_full.csv"`) |
| `output_name` | Yes | - | Base name for results |
| `model` | Yes | - | LLM model for scoring |
| `provider` | Yes | - | API provider |
| `max_workers` | No | 4 | Number of parallel workers |
| `main_weight` | No | 0.5 | Importance of main cell type match (0-1) |
| `sub_weight` | No | 0.5 | Importance of subtype match (0-1) |
| `generate_report` | No | TRUE | Generate HTML report |
| `report_output_path` | No | "uq_batch_report.html" | Path for HTML report |
| `reasoning` | No | NULL | Reasoning effort level ("low", "medium", "high") - only for GPT-5 models |

## Output

### Files Generated

| File | Description |
|------|-------------|
| `{output_name}_{n}_full.csv` | Results from each iteration |
| `{output_name}_similarity.csv` | Similarity scores across iterations |
| `uq_batch_report.html` | HTML visualization report |

### Interpreting Similarity Scores

| Score Range | Interpretation | Action |
|-------------|----------------|--------|
| > 0.9 | High consistency | Robust annotation |
| 0.75 - 0.9 | Moderate consistency | Review recommended |
| < 0.75 | Low consistency | Use [Annotation Boost](annotation-boost.md) or [Subclustering](subclustering-analysis.md) |

### Troubleshooting Low Scores

1. **Review Data**: Check marker gene quality and cluster heterogeneity
2. **Try Advanced Agents**: Use [Annotation Boost Agent](annotation-boost.md) or [Subclustering](subclustering-analysis.md)
3. **Adjust Parameters**: Increase iteration count for more reliable consensus
