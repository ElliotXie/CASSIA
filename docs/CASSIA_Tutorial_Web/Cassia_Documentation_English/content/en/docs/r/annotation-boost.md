---
title: Annotation Boost (Optional)
---

## Overview

Annotation Boost is an advanced validation tool that enhances annotation confidence through multiple iterations of analysis. It's particularly useful for:
- Validating low-confidence annotations
- Getting detailed insights into specific cell clusters
- Resolving ambiguous cell type assignments
- Generating comprehensive validation reports

## Quick Start

```R
runCASSIA_annotationboost(
    full_result_path = "batch_results_summary.csv",
    marker = marker_data,
    cluster_name = "CD4+ T cell",
    major_cluster_info = "Human PBMC",
    output_name = "Cluster1_report",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter",
)
```

## Input

- Full results CSV from CASSIA batch analysis (`_summary.csv`)
- Original marker gene file (***Note: The marker file should not be filtered!***)
- Cluster context information
- Specific cluster identifier
- (Optional) Conversations JSON file from batch annotation (`_conversations.json`)

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `full_result_path` | Path to the CASSIA results CSV file (`_summary.csv`) |
| `marker` | Marker gene data (data frame or path). Use the same marker data as the initial analysis (do not filter) |
| `cluster_name` | Exact name of the target cluster to validate |
| `major_cluster_info` | Context about the dataset (e.g., "Human PBMC", "Mouse Brain") |
| `output_name` | Base name for the output validation report |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `num_iterations` | 5 | Number of validation rounds |
| `model` | - | LLM model to use. Recommended: `anthropic/claude-sonnet-4.5` or better |
| `provider` | - | API provider for the model |
| `conversations_json_path` | - | Path to the conversations JSON file from batch annotation |
| `conversation_history_mode` | `"full"` | How to use prior conversation history: `"full"`, `"final"`, or `"none"` |

### Advanced Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `search_strategy` | `"breadth"` | Strategy for exploring hypotheses: `"breadth"` or `"depth"` |
| `report_style` | `"per_iteration"` | Format of the final report: `"per_iteration"` or `"total_summary"` |
| `validator_involvement` | `"v1"` | Level of validation strictness: `"v1"` (moderate) or `"v0"` (high) |
| `reasoning` | - | Reasoning effort level: `"low"`, `"medium"`, `"high"`. Only supported by OpenAI GPT-5 series models |

## Output

The analysis generates the following output files:
- `{output_name}_summary.html`: HTML report with detailed analysis results and visualizations.
- `{output_name}_raw_conversation.txt`: Raw conversation text containing the full analysis dialogue.
