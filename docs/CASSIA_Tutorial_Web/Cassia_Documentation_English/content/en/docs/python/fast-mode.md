---
title: Fast Mode
---

## Overview

CASSIA's Fast Mode offers a streamlined, one-line solution for running the complete analysis pipeline. This mode combines annotation, scoring, and annotation boost for correcting low quality annotations in a single function call, using optimized default parameters.

## Quick Start

```python
CASSIA.runCASSIA_pipeline(
    output_file_name = "my_analysis",
    tissue = "brain",
    species = "human",
    marker = marker_data,
    max_workers = 4
)
```

## Input

| Input | Type | Description |
|-------|------|-------------|
| `marker` | DataFrame or path | Marker gene data with cluster and gene columns |
| `tissue` | String | Tissue type of the sample (e.g., "brain", "lung") |
| `species` | String | Species of the sample (e.g., "human", "mouse") |

## Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `output_file_name` | String | Base name for the output folder and files |
| `tissue` | String | The tissue type of the sample |
| `species` | String | The species of the sample |
| `marker` | DataFrame/path | Marker gene data |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_workers` | 6 | Number of parallel processes to use |
| `overall_provider` | "openrouter" | Main provider for all pipeline stages ("openai", "anthropic", "openrouter") |
| `score_threshold` | 75 | Annotations below this score (0-100) trigger Annotation Boost |
| `additional_info` | None | Optional experimental context (e.g., "treated with drug X") |
| `validator_involvement` | "v1" | Validation strictness ("v1" = moderate, "v0" = high) |
| `do_merge_annotations` | True | If True, merges detailed cell types into broader categories |
| `annotation_model` | "anthropic/claude-sonnet-4.5" | Model for initial cell type annotation |
| `annotation_provider` | "openrouter" | Provider for annotation model |
| `score_model` | "openai/gpt-5.1" | Model for quality scoring |
| `score_provider` | "openrouter" | Provider for score model |
| `annotationboost_model` | "anthropic/claude-sonnet-4.5" | Model for refining low-confidence annotations |
| `annotationboost_provider` | "openrouter" | Provider for annotation boost model |
| `merge_model` | "google/gemini-2.5-flash" | Model for the merging step |
| `merge_provider` | "openrouter" | Provider for merge model |
| `overall_reasoning` | None | Reasoning effort for all stages ("low", "medium", "high") |
| `annotation_reasoning` | None | Override reasoning level for annotation stage only |
| `score_reasoning` | None | Override reasoning level for scoring stage only |
| `annotationboost_reasoning` | None | Override reasoning level for annotation boost stage only |
| `merge_reasoning` | None | Override reasoning level for merging stage only |

## Output

The pipeline generates a timestamped main folder (`CASSIA_Pipeline_{tissue}_{species}_{timestamp}`) containing three organized subfolders:

### Folder Structure

```
CASSIA_Pipeline_brain_human_20240115_143022/
├── 01_annotation_report/
│   └── {name}_report.html          # Interactive HTML report
├── 02_annotation_boost/
│   └── {cluster_name}/             # One folder per low-scoring cluster
│       └── {name}_{cluster}_boosted_report.html
└── 03_csv_files/
    ├── {name}_summary.csv          # Initial annotation results
    ├── {name}_conversations.json   # Full conversation history
    ├── {name}_scored.csv           # Results with quality scores
    ├── {name}_merged.csv           # Merged annotations (if enabled)
    └── {name}_FINAL_RESULTS.csv    # Combined final results
```

### Output Files

| Folder | File | Description |
|--------|------|-------------|
| `01_annotation_report` | `{name}_report.html` | Interactive HTML report with all annotations |
| `02_annotation_boost` | Per-cluster folders | Boost analysis for clusters scoring below threshold |
| `03_csv_files` | `{name}_FINAL_RESULTS.csv` | **Main output** - combined results with scores and merged annotations |
| `03_csv_files` | `{name}_summary.csv` | Initial cell type annotations |
| `03_csv_files` | `{name}_scored.csv` | Annotations with quality scores |
| `03_csv_files` | `{name}_merged.csv` | Broader category groupings (if `do_merge_annotations = True`) |
| `03_csv_files` | `{name}_conversations.json` | Full LLM conversation history for reproducibility |

### Performance Tips

- For optimal performance, adjust `max_workers` based on your system's CPU cores
- Use `additional_info` to provide relevant experimental context
- Monitor `score_threshold` to balance stringency with throughput

---

Next we introduce each function in detail...
