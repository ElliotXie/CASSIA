---
title: Single Cluster Analysis
---

## Overview

The `runCASSIA` function analyzes a single cluster of marker genes to identify the cell type. This function is specifically designed for users who only have one cluster to analyze.

Note: CASSIA is designed to handle multiple clusters at once via [Batch Processing](batch-processing.md). Use this function when you need to annotate just a single cluster.

---

## Quick Start

```R
result <- runCASSIA(
    marker_list = c("CD3D", "CD3E", "CD2", "TRAC"),
    model = "anthropic/claude-sonnet-4.5",
    tissue = "blood",
    species = "human",
    provider = "openrouter"
)

# View the annotation result
print(result$structured_output)
```

For model recommendations, see [How to Select Models and Providers](setting-up-cassia.md#how-to-select-models-and-providers).

---

## Input

### Marker List Format

Provide a character vector of marker gene names for your cluster:

```R
marker_list <- c("CD3D", "CD3E", "CD2", "TRAC", "IL7R")
```

These should be the top differentially expressed genes that characterize your cluster of interest.

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `marker_list` | Character vector of marker gene names for the cluster |
| `model` | LLM model ID (e.g., `"anthropic/claude-sonnet-4.5"`) |
| `tissue` | Tissue type (e.g., `"blood"`, `"brain"`) |
| `species` | Species (e.g., `"human"`, `"mouse"`) |
| `provider` | API provider (`"openrouter"`, `"openai"`, `"anthropic"`) |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `temperature` | 0 | Output randomness (0=deterministic, 1=creative). Keep at 0 for reproducible results. |
| `additional_info` | `NULL` | Extra experimental context about the sample |
| `validator_involvement` | `"v1"` | Validation intensity: `"v1"` (moderate) or `"v0"` (high, slower) |
| `reasoning` | `NULL` | Reasoning depth for compatible models (`"low"`, `"medium"`, `"high"`). See below. |

### Parameter Details

**Model Selection**
- Default: `anthropic/claude-sonnet-4.5` for best performance
- Alternative: `google/gemini-2.5-flash` for faster analysis
- When using OpenRouter, specify the complete model ID
- See [How to Select Models and Providers](setting-up-cassia.md#how-to-select-models-and-providers) for detailed recommendations

**Reasoning Parameter**
- Controls reasoning depth for compatible models (GPT-5 series via OpenRouter)
- Options: `"low"`, `"medium"`, `"high"`
- Omit this parameter for standard mode
- See [Reasoning Effort Parameter](setting-up-cassia.md#reasoning-effort-parameter) for details

**Additional Context**
- Use `additional_info` to provide experimental context
- Example: `"Sample from tumor microenvironment, focus on immune infiltration"`

---

## Output

The function returns a list with two components:

| Component | Description |
|-----------|-------------|
| `structured_output` | The annotation result containing predicted cell type and reasoning |
| `conversation_history` | Complete conversation log for debugging and transparency |

