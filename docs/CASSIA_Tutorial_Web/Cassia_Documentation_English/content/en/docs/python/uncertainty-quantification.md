---
title: Uncertainty Quantification (Optional)
---

Uncertainty quantification in CASSIA helps assess annotation reliability through multiple analysis iterations and similarity scoring. This process is crucial for:
- Identifying robust cell type assignments
- Detecting mixed or ambiguous clusters
- Quantifying annotation confidence
- Understanding prediction variability

### Single Cluster Uncertainty Analysis

For analyzing uncertainty in a single cluster, use `runCASSIA_n_times_similarity_score()`:

```python
from CASSIA import runCASSIA_n_times_similarity_score

# Run multiple iterations on a single cluster with similarity scoring
result = runCASSIA_n_times_similarity_score(
    tissue="large intestine",
    species="human",
    marker_list=["CD38", "CD138", "JCHAIN", "MZB1", "SDC1"],
    model="google/gemini-2.5-flash",
    provider="openrouter",
    n=5,  # Number of iterations
    temperature=0.3,
    max_workers=3,
    main_weight=0.5,
    sub_weight=0.5,
    validator_involvement="v1"
)

# Access results
print(f"Main cell type: {result['general_celltype_llm']}")
print(f"Sub cell type: {result['sub_celltype_llm']}")
print(f"Similarity score: {result['similarity_score']}")
print(f"Consensus types: {result['consensus_types']}")

# Check for mixed cell types
if result.get('Possible_mixed_celltypes_llm'):
    print(f"Possible mixed types: {result['Possible_mixed_celltypes_llm']}")
```

#### Parameter Details (Single Cluster)

- **`tissue`**: Tissue type for context.
- **`species`**: Species for context.
- **`marker_list`**: List of marker genes for the cluster.
- **`model`**: LLM model to use.
- **`provider`**: API provider ("openrouter", "openai", "anthropic").
- **`n`**: Number of analysis iterations (default: 5).
- **`temperature`**: LLM temperature (lower = more consistent).
- **`max_workers`**: Parallel processing workers.
- **`main_weight`**: Weight for main cell type in similarity (0-1).
- **`sub_weight`**: Weight for subtype in similarity (0-1).
- **`validator_involvement`**: Validator mode ("v0" strict, "v1" moderate).
- **`additional_info`**: Optional additional context string.

#### Return Values (Single Cluster)

The function returns a dictionary containing:
- **`general_celltype_llm`**: Consensus main cell type.
- **`sub_celltype_llm`**: Consensus sub cell type.
- **`similarity_score`**: Overall similarity across iterations (0-1).
- **`consensus_types`**: Cell types that appeared most frequently.
- **`Possible_mixed_celltypes_llm`**: Detected mixed cell type populations.
- **`original_results`**: Raw results from each iteration.

### Batch Iteration Analysis

For batch processing across multiple clusters, use `runCASSIA_batch_n_times`:

```python
# Run multiple iterations
iteration_results = CASSIA.runCASSIA_batch_n_times(
    n=2,
    marker=unprocessed_markers,
    output_name=output_name + "_Uncertainty",
    model="anthropic/claude-sonnet-4.5",
    provider="openrouter",
    tissue="large intestine",
    species="human",
    max_workers=6,
    batch_max_workers=3  # Conservative setting for API rate limits
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
- **`batch_max_workers`**: Workers per iteration.

### Similarity Score Calculation

```python
# Calculate similarity scores
similarity_scores = CASSIA.runCASSIA_similarity_score_batch(
    marker=unprocessed_markers,
    file_pattern=output_name + "_Uncertainty_*_full.csv",
    output_name="intestine_uncertainty",
    max_workers=6,
    model="openai/gpt-5.1",
    provider="openrouter",
    main_weight=0.5,
    sub_weight=0.5
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
