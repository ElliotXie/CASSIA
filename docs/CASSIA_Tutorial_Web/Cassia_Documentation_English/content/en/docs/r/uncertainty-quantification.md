---
title: Uncertainty Quantification (Optional)
---


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
    n = 5,
    marker = marker_data,
    output_name = "my_annotation_repeat",
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter",
    tissue = "brain",
    species = "human",
    max_workers = 4,
    batch_max_workers = 2
)
```
> **⚠️ Cost Warning**: Running multiple iterations with LLM models can incur significant costs. Each iteration makes separate API calls, so the total cost will be approximately n times the cost of a single run. Consider starting with a smaller number of iterations for testing purposes.

#### Parameter Details

- **`n`**: Number of analysis iterations (Recommended: 5).
- **`marker`**: Marker gene data (data frame or path).
- **`output_name`**: Base name for output files.
- **`model`**: LLM model to use.
- **`provider`**: API provider.
- **`tissue`**: Tissue type.
- **`species`**: Species.
- **`max_workers`**: Overall parallel processing limit.
- **`batch_max_workers`**: Workers per iteration (max_workers * batch_max_workers should match your number of cores).

### Similarity Score Calculation

#### Running Similarity Analysis
```R
# Calculate similarity scores
runCASSIA_similarity_score_batch(
    marker = marker_data,
    file_pattern = "my_annotation_repeat_*_full.csv",
    output_name = "similarity_results",
    max_workers = 4,
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

#### Parameter Details

- **`marker`**: Marker gene data.
- **`file_pattern`**: Pattern to match iteration results (use `*` wildcard).
    - Example: `"my_annotation_repeat_*_full.csv"` matches `my_annotation_repeat_1_full.csv`, `my_annotation_repeat_2_full.csv`, etc.
- **`output_name`**: Base name for results.
- **`max_workers`**: Number of parallel workers.
- **`model`**: LLM model for scoring.
- **`provider`**: API provider.
- **`main_weight`**: Importance of main cell type match (0-1).
- **`sub_weight`**: Importance of subtype match (0-1). (Weights should sum to 1.0).
- **`generate_report`**: Whether to generate an HTML report (default: TRUE).
- **`report_output_path`**: Path for the HTML report (default: 'uq_batch_report.html').

> **📊 Automatic Report Generation**: By default, an HTML report is automatically generated with visualizations of the uncertainty analysis results.

### Output Interpretation & Troubleshooting

#### Similarity Scores
- **Range**: 0 (completely different) to 1 (identical).
- **> 0.9**: High consistency (Robust annotation).
- **0.75 - 0.9**: Moderate consistency.
- **< 0.75**: Low consistency (Ambiguous).

#### Troubleshooting Low Scores
1. **Review Data**: Check marker gene quality and cluster heterogeneity.
2. **Try Advanced Agents**: Use the [Annotation Boost Agent](annotation-boost.md) or [Subclustering](subclustering-analysis.md).
3. **Parameters**: Increase iteration count.
