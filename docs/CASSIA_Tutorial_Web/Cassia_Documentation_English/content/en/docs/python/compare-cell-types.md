---
title: Symphony Compare (Optional)
---

## Overview

> **Note:** Symphony Compare requires OpenRouter as the provider. This is the only supported provider for this module.

The `symphonyCompare` agent acts as a virtual panel of experts to resolve ambiguous cell type annotations. It orchestrates multiple AI models (a "Symphony" of agents) to compare potential cell types, debate their findings in discussion rounds, and reach a consensus based on marker gene evidence.

This agent is particularly useful after running the default CASSIA pipeline if you are unsure about a specific cluster's identity. For example, distinguishing between different subtypes of Plasma Cells.

## Quick Start

```python
results = CASSIA.symphonyCompare(
    tissue = "large intestine",
    celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells"],
    marker_set = "IGLL5, IGLV6-57, JCHAIN, IGKC, TNFRSF17, IGHG1, MZB1",
    species = "human"
)

print(f"Consensus: {results['consensus']} (confidence: {results['confidence']:.1%})")
```

## Input

- **Marker genes**: A string of comma-separated marker genes from CASSIA's previous results or your own analysis
- **Candidate cell types**: A list of 2-4 cell types to compare
- **Tissue context**: The tissue type being analyzed
- **Species**: The species of the sample (default: human)

## Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `tissue` | string | The tissue type being analyzed (e.g., "large intestine") |
| `celltypes` | list | A list of 2-4 cell types to compare |
| `marker_set` | string | A string of comma-separated marker genes |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `species` | string | `"human"` | The species of the sample |
| `model_preset` | string | `"premium"` | Configuration of models to use. `"premium"`: High-performance ensemble (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4). `"budget"`: Cost-effective models (DeepSeek V3.2, Grok 4 Fast, Kimi K2, Gemini 2.5 Flash) |
| `output_basename` | string | - | Base name for output files |
| `enable_discussion` | bool | `True` | Whether to enable multi-round debate between models |

### Advanced Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_discussion_rounds` | int | `2` | Maximum number of discussion rounds |
| `consensus_threshold` | float | `0.8` | Fraction of models required for consensus (0-1) |

## Output

The function returns a dictionary containing consensus results and generates output files.

**Return Value:**
- `consensus`: The consensus cell type reached by the model panel
- `confidence`: Confidence level of the consensus (0-1)

**Generated Files:**
- `{output_basename}.csv`: Detailed comparison results, reasoning, and scores from all models and rounds
- `{output_basename}_report.html`: An interactive HTML report visualizing the debate and consensus process
