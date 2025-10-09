# Test: Uncertainty Quantification

## Overview

**Method**: Multiple stochastic runs of `runCASSIA_batch()`
**Purpose**: Test annotation stability and confidence through repeated analysis
**Category**: Quality Assessment
**Priority**: Medium
**Expected Runtime**: 15-20 minutes

## Description

This test validates **annotation stability and uncertainty quantification** by running batch annotation multiple times with stochastic settings (temperature > 0) and analyzing consistency across runs. This is critical for:

- Assessing annotation confidence
- Identifying unstable/ambiguous cell types
- Understanding LLM variability
- Building consensus annotations
- Validating annotation reliability

### What This Test Validates

- ✅ Multiple independent annotation runs
- ✅ Pairwise similarity calculation
- ✅ Per-cluster consistency metrics
- ✅ Overall stability scoring
- ✅ Consensus annotation generation
- ✅ Confidence score calculation
- ✅ Annotation variance analysis

## Prerequisites

### Required

- Python 3.8+
- CASSIA package installed
- OpenRouter API key
- ~15-20 minutes execution time (5 runs × 3-4 min each)
- Sample data with at least 3 clusters

### Environment Setup

```bash
# Set API key
export OPENROUTER_API_KEY='your-api-key-here'

# Verify CASSIA is installed
python -c "import CASSIA; print(CASSIA.__version__)"
```

## Quick Start

```bash
# Navigate to test directory
cd test/05_uncertainty_quantification/

# Run test (will run 5 independent stochastic runs)
python test_uq_batch.py
```

## Configuration

The test is configured via `config.json`:

```json
{
  "test_name": "test_uncertainty_quantification",
  "description": "Test annotation stability through multiple stochastic runs",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "processed.csv",
  "method_params": {
    "n_iterations": 5,
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4,
    "n_genes": 50,
    "n_clusters_subset": 3,
    "output_name": "uq_test"
  },
  "validation": {
    "check_all_runs_complete": true,
    "check_similarity_scores": true,
    "min_similarity": 0.5,
    "max_similarity": 1.0,
    "check_stability_metrics": true
  }
}
```

### Key Parameters

- **n_iterations** (5): Number of independent runs
- **temperature** (0.7): Stochasticity level (0.5-1.0 recommended for UQ)
- **n_clusters_subset** (3): Use subset for faster testing
- **min_similarity/max_similarity**: Expected similarity range for validation

### Customization

Edit `config.json` to:

- **Increase runs** for better statistics (`n_iterations: 10`)
- **Adjust temperature** to control variability (0.5 = less, 1.0 = more)
- **Test more clusters** by changing `n_clusters_subset`
- **Change similarity thresholds** based on expected stability

## Expected Outputs

### Console Output

```
================================================================================
CASSIA Test: Uncertainty Quantification (Multiple Stochastic Runs)
================================================================================
✓ API key validated for provider: openrouter
✓ Loaded subset: 3 rows (3 clusters)
  - Clusters: ['monocyte', 'plasma cell', 'cd8_positive_alpha_beta_t_cell']

================================================================================
Test configuration:
  Test: test_uncertainty_quantification
  Model: google/gemini-2.5-flash-preview
  Temperature: 0.7
================================================================================

================================================================================
RUNNING MULTIPLE STOCHASTIC ITERATIONS
================================================================================

UQ configuration:
  - Number of runs: 5
  - Temperature: 0.7 (stochastic)
  - Clusters: 3
  - Model: google/gemini-2.5-flash-preview

  Running iteration 1/5...
    ✓ Iteration 1 completed: 3 clusters annotated
  Running iteration 2/5...
    ✓ Iteration 2 completed: 3 clusters annotated
  Running iteration 3/5...
    ✓ Iteration 3 completed: 3 clusters annotated
  Running iteration 4/5...
    ✓ Iteration 4 completed: 3 clusters annotated
  Running iteration 5/5...
    ✓ Iteration 5 completed: 3 clusters annotated

Starting: All 5 UQ iterations
Completed: All 5 UQ iterations (645.23s = 10.75min)

✓ All 5 iterations completed

================================================================================
ANALYZING ANNOTATION STABILITY
================================================================================

Calculating pairwise similarity scores:
  Run 1 vs Run 2: 1.000 (3/3 matches)
  Run 1 vs Run 3: 0.667 (2/3 matches)
  Run 1 vs Run 4: 1.000 (3/3 matches)
  Run 1 vs Run 5: 0.667 (2/3 matches)
  Run 2 vs Run 3: 0.667 (2/3 matches)
  Run 2 vs Run 4: 1.000 (3/3 matches)
  Run 2 vs Run 5: 0.667 (2/3 matches)
  Run 3 vs Run 4: 0.667 (2/3 matches)
  Run 3 vs Run 5: 1.000 (3/3 matches)
  Run 4 vs Run 5: 0.667 (2/3 matches)

Calculating stability metrics:
  - Average pairwise similarity: 0.767
  - monocyte: 1.000 consistency (5/5 agreement on 'Classical Monocyte')
  - plasma cell: 1.000 consistency (5/5 agreement on 'Plasma Cell')
  - cd8_positive_alpha_beta_t_cell: 0.600 consistency (3/5 agreement on 'Cytotoxic T Cell')

  Overall stability score: 0.867

Creating consensus annotations:
  - monocyte: 'Classical Monocyte' (confidence: 1.000)
  - plasma cell: 'Plasma Cell' (confidence: 1.000)
  - cd8_positive_alpha_beta_t_cell: 'Cytotoxic T Cell' (confidence: 0.600)

================================================================================
VALIDATING UQ RESULTS
================================================================================
  ✓ All 5 runs completed
  ✓ Average similarity within expected range: 0.767 ∈ [0.5, 1.0]
  ✓ Stability score calculated: 0.867

================================================================================
ARCHIVING RESULTS
================================================================================
  ✓ Archived run 1: 20251007_150234_run_1_annotations.csv
  ✓ Archived run 2: 20251007_150234_run_2_annotations.csv
  ✓ Archived run 3: 20251007_150234_run_3_annotations.csv
  ✓ Archived run 4: 20251007_150234_run_4_annotations.csv
  ✓ Archived run 5: 20251007_150234_run_5_annotations.csv
  ✓ Archived similarity matrix: 20251007_150234_similarity_matrix.csv
  ✓ Archived consensus annotations: 20251007_150234_consensus_annotations.csv
  ✓ Archived stability metrics: 20251007_150234_stability_metrics.json

================================================================================
UNCERTAINTY QUANTIFICATION SUMMARY
================================================================================
Clusters analyzed: 3
Stochastic runs: 5
Average similarity: 0.767
Overall stability: 0.867
Similarity range: 0.667 - 1.000

Per-cluster confidence:
  - monocyte: 1.000 (5/5 consensus)
  - plasma cell: 1.000 (5/5 consensus)
  - cd8_positive_alpha_beta_t_cell: 0.600 (3/5 consensus)

Test completed successfully (645.23s = 10.75min)

================================================================================
✅ UNCERTAINTY QUANTIFICATION TEST COMPLETED SUCCESSFULLY
📁 Results directory: results/
📄 Log file: results/20251007_150234_test_log.txt
================================================================================
```

### Generated Files

```
test/05_uncertainty_quantification/
│
├── config.json                                    # Test configuration
├── test_uq_batch.py                               # Test script
├── README.md                                      # This file
│
└── results/                                       # Timestamped results
    ├── 20251007_150234_run_1_annotations.csv      # Run 1 results
    ├── 20251007_150234_run_2_annotations.csv      # Run 2 results
    ├── 20251007_150234_run_3_annotations.csv      # Run 3 results
    ├── 20251007_150234_run_4_annotations.csv      # Run 4 results
    ├── 20251007_150234_run_5_annotations.csv      # Run 5 results
    ├── 20251007_150234_similarity_matrix.csv      # Pairwise similarities
    ├── 20251007_150234_consensus_annotations.csv  # Consensus results
    ├── 20251007_150234_stability_metrics.json     # Detailed metrics
    └── 20251007_150234_test_log.txt              # Detailed log
```

### Similarity Matrix CSV

| | Run_1 | Run_2 | Run_3 | Run_4 | Run_5 |
|---------|-------|-------|-------|-------|-------|
| Run_1 | 1.000 | 1.000 | 0.667 | 1.000 | 0.667 |
| Run_2 | 1.000 | 1.000 | 0.667 | 1.000 | 0.667 |
| Run_3 | 0.667 | 0.667 | 1.000 | 0.667 | 1.000 |
| Run_4 | 1.000 | 1.000 | 0.667 | 1.000 | 0.667 |
| Run_5 | 0.667 | 0.667 | 1.000 | 0.667 | 1.000 |

### Consensus Annotations CSV

| True Cell Type | Consensus_Annotation | Confidence_Score |
|----------------------------------------|----------------------|------------------|
| monocyte | Classical Monocyte | 1.000 |
| plasma cell | Plasma Cell | 1.000 |
| cd8_positive_alpha_beta_t_cell | Cytotoxic T Cell | 0.600 |

### Stability Metrics JSON

```json
{
  "n_runs": 5,
  "avg_similarity": 0.767,
  "stability_score": 0.867,
  "similarity_range": {
    "min": 0.667,
    "max": 1.000,
    "std": 0.164
  },
  "cluster_consistency": {
    "monocyte": {
      "consistency": 1.0,
      "most_common": "Classical Monocyte",
      "count": 5,
      "total_runs": 5,
      "unique_annotations": 1
    },
    "plasma cell": {
      "consistency": 1.0,
      "most_common": "Plasma Cell",
      "count": 5,
      "total_runs": 5,
      "unique_annotations": 1
    },
    "cd8_positive_alpha_beta_t_cell": {
      "consistency": 0.6,
      "most_common": "Cytotoxic T Cell",
      "count": 3,
      "total_runs": 5,
      "unique_annotations": 2
    }
  }
}
```

## Uncertainty Quantification Workflow

### How It Works

1. **Subset Selection**: Use small cluster subset for faster testing
2. **Multiple Runs**: Execute `runCASSIA_batch()` N times independently
3. **Stochastic Sampling**: Use temperature > 0 to enable variability
4. **Similarity Calculation**: Compare all pairs of runs
5. **Consensus Building**: Identify most common annotation per cluster
6. **Confidence Scoring**: Calculate agreement proportion
7. **Metrics Reporting**: Generate comprehensive stability metrics

### Metrics Explained

**Pairwise Similarity**:
- Proportion of clusters with identical annotations between two runs
- Range: [0, 1] where 1 = perfect agreement
- Formula: `matches / total_clusters`

**Per-Cluster Consistency**:
- Proportion of runs agreeing on most common annotation
- Range: [0, 1] where 1 = unanimous agreement
- Formula: `max_count / n_runs`

**Overall Stability Score**:
- Average of all per-cluster consistency scores
- Range: [0, 1] where 1 = perfect stability
- Formula: `mean(cluster_consistencies)`

**Average Similarity**:
- Mean of all pairwise similarities (excluding diagonal)
- Range: [0, 1] where 1 = all runs identical
- Formula: `mean(upper_triangle(similarity_matrix))`

## Validation Checks

### 1. Run Completion

- ✅ All N runs completed successfully
- ✅ Each run produced valid results
- ✅ No crashes or errors during execution

### 2. Similarity Range

- ✅ Average similarity within expected bounds
- ✅ Not too low (< 0.3 suggests random annotations)
- ✅ Not too high (> 0.95 suggests deterministic behavior)

### 3. Stability Metrics

- ✅ Stability score calculated correctly
- ✅ Per-cluster metrics reasonable
- ✅ Consensus annotations identified

## Sample Results

### High Stability Example

```
monocyte:
  - Run 1: "Classical Monocyte"
  - Run 2: "Classical Monocyte"
  - Run 3: "Classical Monocyte"
  - Run 4: "Classical Monocyte"
  - Run 5: "Classical Monocyte"

  Consistency: 1.0 (5/5 agreement)
  Confidence: High ✅
  Interpretation: Very stable annotation
```

### Low Stability Example

```
cd8_positive_alpha_beta_t_cell:
  - Run 1: "Cytotoxic T Cell"
  - Run 2: "CD8+ T Cell"
  - Run 3: "Cytotoxic T Cell"
  - Run 4: "Cytotoxic CD8+ T Cell"
  - Run 5: "Cytotoxic T Cell"

  Consistency: 0.6 (3/5 agreement on "Cytotoxic T Cell")
  Confidence: Medium ⚠️
  Interpretation: Unstable annotation, needs boost or more markers
```

## Troubleshooting

### 1. Too High Similarity (> 0.95)

**Issue**: Runs are too similar, not capturing variability

**Cause**: Temperature too low or deterministic behavior

**Solutions**:
```json
{
  "temperature": 1.0  // Increase from 0.7
}
```

Or verify model truly uses temperature parameter.

### 2. Too Low Similarity (< 0.3)

**Issue**: Runs are too different, annotations seem random

**Causes**:
- Poor marker quality
- Temperature too high
- Ambiguous cell types

**Solutions**:
```json
{
  "temperature": 0.5,  // Reduce from 0.7
  "n_genes": 100       // Use more markers
}
```

Or use better quality marker data (unprocessed.csv).

### 3. Test Takes Too Long

**Issue**: 5 runs × 4 min each = 20 min too slow

**Solutions**:
```json
{
  "n_iterations": 3,        // Fewer runs
  "n_clusters_subset": 2,   // Fewer clusters
  "max_workers": 8          // More parallelism
}
```

### 4. All Clusters Have Low Consistency

**Issue**: Every cluster shows < 0.7 consistency

**Possible Causes**:
- Markers not distinctive enough
- Model not suited for task
- Temperature too high

**Solutions**:
- Use unprocessed.csv with more genes
- Try different model
- Lower temperature to 0.5
- Run annotation boost on unstable clusters

### 5. Memory Issues with Many Runs

**Issue**: Out of memory with n_iterations > 10

**Solution**: Process runs sequentially or reduce data size:
```python
# In test script, clear memory after each run:
import gc
results_list.append(result_df.copy())
del result_df
gc.collect()
```

## Advanced Usage

### Increase Number of Runs

For better statistical power:

```json
{
  "n_iterations": 10  // More robust estimates
}
```

Runtime: ~30-40 minutes for 10 runs

### Test Full Dataset

Remove subset restriction:

```python
# In test_uq_batch.py, modify:
marker_data = loader.load_processed()  # Full dataset
# Remove subset call
```

### Different Temperature Values

Test multiple temperature settings:

```python
temperatures = [0.3, 0.5, 0.7, 0.9]

for temp in temperatures:
    config['temperature'] = temp
    # Run test and compare stability across temperatures
```

### Bootstrap Confidence Intervals

Calculate confidence intervals for metrics:

```python
from scipy import stats

# Bootstrap sampling from runs
n_bootstrap = 1000
bootstrap_similarities = []

for _ in range(n_bootstrap):
    # Resample runs with replacement
    sampled_runs = np.random.choice(len(results_list), size=len(results_list), replace=True)
    # Calculate similarity for this sample
    # ...

# Calculate 95% CI
ci = stats.t.interval(0.95, len(bootstrap_similarities)-1,
                      loc=np.mean(bootstrap_similarities),
                      scale=stats.sem(bootstrap_similarities))
```

## Performance Benchmarks

### Expected Performance (gemini-2.5-flash-preview)

| Configuration | Time | Operations |
|---------------|------|------------|
| 3 clusters, 5 runs | 15-20 min | 15 annotations |
| 6 clusters, 5 runs | 20-25 min | 30 annotations |
| 3 clusters, 10 runs | 30-40 min | 30 annotations |
| 6 clusters, 10 runs | 40-50 min | 60 annotations |

### Resource Usage

- **Memory**: 300-600 MB (stores all run results)
- **Disk Space**: ~50 MB per 5-run test
- **Network**: ~50-100 KB per annotation call
- **CPU**: Single-threaded (sequential runs)

### Optimization Tips

1. **Use subset** for faster testing (3-5 clusters recommended)
2. **Reduce iterations** to 3 for quick validation
3. **Increase max_workers** for faster individual runs
4. **Cache results** to avoid re-running if interrupted

## Interpreting Results

### Confidence Scores

| Score | Interpretation | Action |
|-------|----------------|--------|
| > 0.9 | Very High | Trust annotation, use with confidence |
| 0.7-0.9 | High | Reliable annotation, consider single verification |
| 0.5-0.7 | Medium | Uncertain, recommend annotation boost |
| < 0.5 | Low | Unstable, requires manual review or better markers |

### Stability Categories

**Stable Cluster** (consistency ≥ 0.8):
- All runs agree on same or very similar annotations
- High confidence in cell type identity
- Can use consensus annotation directly

**Moderately Stable** (0.6 ≤ consistency < 0.8):
- Some variation across runs
- Consensus exists but not unanimous
- Consider annotation boost for confirmation

**Unstable Cluster** (consistency < 0.6):
- High variability across runs
- No clear consensus annotation
- Requires intervention:
  - Run annotation boost
  - Add more markers
  - Manual expert review

## Use Cases

### When to Use UQ Testing

1. **Pre-publication validation** of annotation results
2. **Method comparison** (compare stability across models)
3. **Marker set optimization** (test which markers give stable results)
4. **Quality control** for large-scale annotation projects
5. **Building annotation databases** with confidence scores

### Downstream Applications

**Filter Low-Confidence Annotations**:
```python
consensus_df = pd.read_csv('consensus_annotations.csv')
high_confidence = consensus_df[consensus_df['Confidence_Score'] >= 0.8]
```

**Identify Clusters Needing Boost**:
```python
unstable_clusters = consensus_df[consensus_df['Confidence_Score'] < 0.7]
for cluster in unstable_clusters['True Cell Type']:
    CASSIA.runCASSIA_annotationboost(cluster_name=cluster, ...)
```

**Report Confidence in Publications**:
```
"Cell type annotations were determined with an average confidence of 0.87
(range: 0.60-1.00) based on 5 independent stochastic runs."
```

## Related Tests

- **01_runCASSIA_batch** - Single batch run (deterministic)
- **03_annotation_boost** - Deep-dive for unstable clusters
- **02_runCASSIA_pipeline** - Full pipeline with built-in scoring

## References

- [CASSIA Documentation](../../data/CASSIA_Package_Documentation.md)
- [Testing Framework Plan](../../../../dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md)
- [Uncertainty Quantification in ML](https://arxiv.org/abs/1906.08158)
- [Annotation Confidence in scRNA-seq](https://doi.org/10.1038/s41592-019-0535-3)

---

**Last Updated**: 2025-10-07
**Version**: 1.0
**Status**: ✅ Active
