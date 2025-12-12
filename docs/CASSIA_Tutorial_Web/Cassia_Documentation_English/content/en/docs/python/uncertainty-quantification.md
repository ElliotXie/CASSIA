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

### Single Cluster Analysis

```python
from CASSIA import runCASSIA_n_times_similarity_score

result = runCASSIA_n_times_similarity_score(
    tissue="large intestine",
    species="human",
    marker_list=["CD38", "CD138", "JCHAIN", "MZB1", "SDC1"],
    model="openai/gpt-5.1",
    provider="openrouter",
    n=5,
    reasoning="medium"
)

print(f"Main cell type: {result['general_celltype_llm']}")
print(f"Similarity score: {result['similarity_score']}")
```

### Batch Analysis

```python
import CASSIA

# Step 1: Run multiple iterations
CASSIA.runCASSIA_batch_n_times(
    n=5,
    marker=marker_data,
    output_name="my_annotation",
    model="openai/gpt-5.1",
    provider="openrouter",
    tissue="large intestine",
    species="human",
    reasoning="medium"
)

# Step 2: Calculate similarity scores
CASSIA.runCASSIA_similarity_score_batch(
    marker=marker_data,
    file_pattern="my_annotation_*_full.csv",
    output_name="similarity_results",
    model="openai/gpt-5.1",
    provider="openrouter",
    reasoning="medium"
)
```

## Input

| Input | Description | Format |
|-------|-------------|--------|
| `marker_list` | Marker genes for single cluster | List of gene names |
| `marker` | Marker gene data for batch | DataFrame or file path |
| `tissue` | Tissue type context | String (e.g., "brain", "large intestine") |
| `species` | Species context | String (e.g., "human", "mouse") |
| `file_pattern` | Pattern to match iteration results | Glob pattern with `*` wildcard |

## Parameters

### Single Cluster (`runCASSIA_n_times_similarity_score`)

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `tissue` | Yes | - | Tissue type for context |
| `species` | Yes | - | Species for context |
| `marker_list` | Yes | - | List of marker genes |
| `model` | Yes | - | LLM model to use |
| `provider` | Yes | - | API provider ("openrouter", "openai", "anthropic") |
| `n` | No | 5 | Number of analysis iterations |
| `temperature` | No | 0.3 | LLM temperature (lower = more consistent) |
| `max_workers` | No | 3 | Parallel processing workers |
| `main_weight` | No | 0.5 | Weight for main cell type in similarity (0-1) |
| `sub_weight` | No | 0.5 | Weight for subtype in similarity (0-1) |
| `validator_involvement` | No | "v1" | Validator mode ("v0" strict, "v1" moderate) |
| `additional_info` | No | None | Additional context string |
| `generate_report` | No | True | Generate HTML report |
| `report_output_path` | No | "uq_report.html" | Path for HTML report |
| `reasoning` | No | None | Reasoning effort level ("low", "medium", "high") - only for GPT-5 models |

### Batch Iteration (`runCASSIA_batch_n_times`)

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `n` | Yes | - | Number of analysis iterations (recommended: 5) |
| `marker` | Yes | - | Marker gene data (DataFrame or path) |
| `output_name` | Yes | - | Base name for output files |
| `model` | Yes | - | LLM model to use |
| `provider` | Yes | - | API provider |
| `tissue` | Yes | - | Tissue type |
| `species` | Yes | - | Species |
| `max_workers` | No | 4 | Overall parallel processing limit |
| `batch_max_workers` | No | 2 | Workers per iteration |
| `reasoning` | No | None | Reasoning effort level ("low", "medium", "high") - only for GPT-5 models |

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
| `generate_report` | No | True | Generate HTML report |
| `report_output_path` | No | "uq_batch_report.html" | Path for HTML report |
| `reasoning` | No | None | Reasoning effort level ("low", "medium", "high") - only for GPT-5 models |

## Output

### Files Generated

| File | Description |
|------|-------------|
| `{output_name}_{n}_full.csv` | Results from each iteration |
| `{output_name}_similarity.csv` | Similarity scores across iterations |
| `uq_report.html` / `uq_batch_report.html` | HTML visualization report |

### Return Values (Single Cluster)

| Key | Description |
|-----|-------------|
| `general_celltype_llm` | Consensus main cell type |
| `sub_celltype_llm` | Consensus sub cell type |
| `similarity_score` | Overall similarity across iterations (0-1) |
| `consensus_types` | Cell types that appeared most frequently |
| `Possible_mixed_celltypes_llm` | Detected mixed cell type populations |
| `original_results` | Raw results from each iteration |

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
