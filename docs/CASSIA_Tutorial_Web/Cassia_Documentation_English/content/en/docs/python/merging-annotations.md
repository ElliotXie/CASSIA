---
title: Merging Annotations (Optional)
---

Merging annotations groups detailed cell type annotations into broader categories at multiple levels of granularity. This is useful for visualization and downstream analysis when you need hierarchical groupings of cell types.

### Overview

CASSIA provides two merging functions:
- `merge_annotations()`: Merge at a single detail level
- `merge_annotations_all()`: Merge at all three levels simultaneously

### Detail Levels

| Level | Column Name | Description | Example |
|-------|-------------|-------------|---------|
| Broad | Merged_Grouping_1 | General categories | "T cells", "Myeloid cells" |
| Detailed | Merged_Grouping_2 | Intermediate specificity | "CD4 T cells", "Macrophages" |
| Very Detailed | Merged_Grouping_3 | Normalized specific names | "CD4+ naive T cells" |

### Single Level Merging

Use `merge_annotations()` to merge at a specific detail level:

```python
from CASSIA import merge_annotations

# Merge at broad level
result_df = merge_annotations(
    csv_path="annotation_results_summary.csv",
    output_path="merged_broad.csv",
    provider="openrouter",
    model="google/gemini-2.5-flash",
    detail_level="broad",  # Options: "broad", "detailed", "very_detailed"
    batch_size=10
)

# Check results
print(result_df[['Cluster ID', 'main_cell_type', 'Merged_Grouping_1']])
```

### Multi-Level Merging (All at Once)

Use `merge_annotations_all()` to merge at all three levels simultaneously with parallel processing:

```python
from CASSIA import merge_annotations_all

# Merge at all levels simultaneously
result_df = merge_annotations_all(
    csv_path="annotation_results_summary.csv",
    output_path="merged_all.csv",
    provider="openrouter",
    model="google/gemini-2.5-flash",
    batch_size=10,
    reasoning="low"  # Optional: "low", "medium", "high"
)

# Result contains all three grouping columns
print(result_df[['main_cell_type', 'Merged_Grouping_1', 'Merged_Grouping_2', 'Merged_Grouping_3']])
```

### Parameter Details

- **`csv_path`**: Path to CASSIA batch results CSV file (output from `runCASSIA_batch`).
- **`output_path`**: Path for output CSV file with merged groupings.
- **`provider`**: API provider ("openrouter", "openai", "anthropic").
- **`model`**: LLM model to use for grouping decisions.
- **`detail_level`**: Level of grouping ("broad", "detailed", "very_detailed") - only for `merge_annotations()`.
- **`batch_size`**: Number of cell types to process per LLM call (default: 10).
- **`reasoning`**: (Optional) Reasoning effort level ("low", "medium", "high"). Controls how much the model "thinks" before responding. Only supported by OpenAI GPT-5 series models (e.g., `gpt-5.1`). Via OpenRouter, no additional verification needed. Via direct OpenAI API, identity verification (KYC) is required.

### Example Output

For a plasma cell annotation:

| Detail Level | Grouping |
|--------------|----------|
| Broad | B cells / Plasma cells |
| Detailed | Plasma cells |
| Very Detailed | Plasma cells |

For a CD8-positive, alpha-beta T cell:

| Detail Level | Grouping |
|--------------|----------|
| Broad | T cells |
| Detailed | CD8 T cells |
| Very Detailed | CD8+ alpha-beta T cells |

### Use Cases

1. **Visualization**: Create summary plots at different granularity levels.
2. **Cross-dataset comparison**: Compare results using normalized cell type names.
3. **Downstream analysis**: Group cells by broad categories for statistical analysis.

### Integration with Pipeline

Merging is automatically available in `runCASSIA_pipeline()` via the `do_merge_annotations` parameter:

```python
CASSIA.runCASSIA_pipeline(
    ...,
    do_merge_annotations=True,  # Enable automatic merging
    merge_model="google/gemini-2.5-flash",
    merge_provider="openrouter"
)
```

### Notes

- The merging process uses LLM-based grouping, so results may vary slightly between runs.
- Lower temperature (e.g., 0.3) is recommended for more consistent groupings.
- Parallel processing in `merge_annotations_all()` uses ThreadPoolExecutor with 3 workers by default.
