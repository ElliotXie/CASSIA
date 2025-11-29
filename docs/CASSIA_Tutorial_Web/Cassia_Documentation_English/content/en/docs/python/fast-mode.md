---
title: Fast Mode
---

CASSIA's Fast Mode offers a streamlined, one-line solution for running the complete analysis pipeline. This mode combines annotation, scoring, and annotation boost for correcting low quality annotations in a single function call, using optimized default parameters.

### Basic Usage

```python
# Run the CASSIA pipeline in fast mode
CASSIA.runCASSIA_pipeline(
    output_file_name = "FastAnalysisResults",
    tissue = "large intestine",
    species = "human",
    marker_path = unprocessed_markers,
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

### Parameter Details

- **`output_file_name`**: Base name for the output folder and files.
- **`tissue`**: The tissue type of the sample (e.g., "brain").
- **`species`**: The species of the sample (e.g., "human").
- **`marker_path`**: Marker gene data (data frame or path to CSV).
- **`max_workers`**: Number of parallel processes to use.
- **`annotation_model`**: Model used for the initial cell type annotation step.
- **`annotation_provider`**: Provider for the annotation model.
- **`score_model`**: Model used for quality scoring. **Recommendation**: Use a high-capability model like `claude-4.5-sonnet` for accurate scoring.
- **`score_provider`**: Provider for the score model.
- **`score_threshold`**: Annotations with a quality score below this threshold (0-100) will trigger the Annotation Boost process. Default is 75.
- **`annotationboost_model`**: Model used for refining low-confidence annotations.
- **`annotationboost_provider`**: Provider for the annotation boost model.
- **`do_merge_annotations`**: Logical. If `True`, merges detailed cell types into broader categories.
- **`merge_model`**: Model used for the merging step.
- **`merge_provider`**: Provider for the merge model.
- **`additional_info`**: Optional experimental context (e.g., "treated with drug X").
- **`validator_involvement`**: Controls validation strictness ("v1" = moderate, "v0" = high).

### Key Features Explanation

#### Merging Annotations
The pipeline includes an automatic merging step (`do_merge_annotations = True`) that groups annotated clusters into broader categories (e.g., grouping "CD4+ T cells" and "CD8+ T cells" into "T cells"). This provides a hierarchical view of your cell types, making it easier to understand the major populations in your dataset.

#### Validator Involvement
The `validator_involvement` parameter controls the intensity of the validation process:
- `"v0"`: High involvement. Stronger validation checks are applied, which may be slower but more rigorous.
- `"v1"`: Moderate involvement (Default). Balanced validation suitable for most standard analyses.

### Output Files
The pipeline generates a folder which contains the following files:
1. Annotation results csv files
2. Scored results csv files
3. Merged annotation results (if enabled)
4. Basic CASSIA report
5. Annotation boost report

### Performance Tips
- For optimal performance, adjust `max_workers` based on your system's CPU cores
- Use `additional_info` to provide relevant experimental context
- Monitor `score_threshold` to balance stringency with throughput
