# Test: runCASSIA_subclusters

## Overview

**Method**: `CASSIA.runCASSIA_subclusters()`
**Purpose**: Test hierarchical subcluster annotation within major cell types
**Category**: Hierarchical Analysis
**Priority**: Medium
**Expected Runtime**: 3-5 minutes

## Description

This test validates the **subclustering functionality** which annotates sub-populations within a major cluster type. This enables hierarchical cell type annotation, identifying refined cell states or subtypes within broad categories.

### Use Cases

- Identifying CD8+ T cell subtypes (naive, effector, memory, exhausted)
- Distinguishing macrophage polarization states (M1, M2, etc.)
- Characterizing epithelial cell maturation states
- Resolving heterogeneity within major cell lineages

### What This Test Validates

- ✅ Subcluster marker processing
- ✅ Top 2 cell type prediction per subcluster
- ✅ Key marker identification
- ✅ Biological reasoning generation
- ✅ HTML report creation
- ✅ CSV output format

## Prerequisites

### Required

- Python 3.8+
- CASSIA package installed
- OpenRouter API key
- Subcluster marker data (2-column format: Subcluster ID, Top Markers)

### Environment Setup

```bash
# Set API key
export OPENROUTER_API_KEY='your-api-key-here'

# Verify CASSIA is installed
python -c "import CASSIA; from CASSIA.subclustering import runCASSIA_subclusters"
```

## Quick Start

```bash
# Navigate to test directory
cd test/06_subclustering/

# Run test
python test_subclustering.py
```

## Configuration

```json
{
  "test_name": "test_subclustering",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0,
  "method_params": {
    "major_cluster_info": "CD8 T cells",
    "output_name": "subcluster_test_results",
    "n_genes": 50
  },
  "validation": {
    "expected_columns": [
      "Result ID",
      "main_cell_type",
      "sub_cell_type",
      "key_markers",
      "reason"
    ]
  }
}
```

### Key Parameters

- **major_cluster_info**: Description of parent cluster (e.g., "CD8 T cells", "macrophages")
- **n_genes**: Number of top genes to include per subcluster
- **output_name**: Base name for CSV and HTML outputs

## Expected Outputs

### Console Output

```
================================================================================
CASSIA Test: runCASSIA_subclusters (Hierarchical Annotation)
================================================================================
✓ API key validated for provider: openrouter

Preparing subcluster data:
Creating synthetic subcluster marker data...
  ✓ Created 4 synthetic subclusters
    - Subcluster_0: IL7R, CD8A, CD8B, CCL4, KLRB1, ITK, LEF1, SELL...
    - Subcluster_1: LAYN, HAVCR2, TIGIT, IKZF2, KLRC2, KLRC3, PDCD1...
    - Subcluster_2: GZMK, GZMH, PRF1, NKG7, CCR7, CD27, GZMA, GNLY...
    - Subcluster_3: WFDC2, CEACAM7, CLDN8, PPARG, MKI67, PCNA, TOP2A...

================================================================================
RUNNING SUBCLUSTER ANNOTATION
================================================================================

Subcluster configuration:
  - Major cluster: CD8 T cells
  - Number of subclusters: 4
  - Model: google/gemini-2.5-flash-preview
  - Number of genes: 50

Starting: Subcluster annotation
Completed: Subcluster annotation (125.43s = 2.09min)

✓ Subcluster annotation completed
  - Results shape: (4, 5)
  - Output file: subcluster_test_results.csv

================================================================================
VALIDATING SUBCLUSTER RESULTS
================================================================================
Validating subcluster results:
  ✓ All expected columns present: ['Result ID', 'main_cell_type', ...]
  ✓ All main_cell_type values non-null
  ✓ All sub_cell_type values non-null
  ✓ Reasons provided: 4/4 subclusters

Subcluster annotation results:
--------------------------------------------------------------------------------

1. Subcluster annotation:
   Main type: Memory CD8+ T cells
   Sub type: Central memory CD8+ T cells
   Key markers: IL7R, CD8A, CD8B, CCL4, KLRB1, ITK
   Reason: IL7R and CCL4 expression with CD8 markers indicates...

2. Subcluster annotation:
   Main type: Exhausted CD8+ T cells
   Sub type: Terminally exhausted CD8+ T cells
   Key markers: LAYN, HAVCR2, TIGIT, IKZF2, KLRC2, KLRC3
   Reason: HAVCR2 (TIM-3) and TIGIT co-expression with LAYN...

...

✓ HTML report generated: subcluster_test_results.html

================================================================================
SUBCLUSTER ANNOTATION SUMMARY
================================================================================
Major cluster: CD8 T cells
Subclusters annotated: 4
Unique main types: 3
Unique sub types: 4

✅ SUBCLUSTER ANNOTATION TEST COMPLETED SUCCESSFULLY
```

### Generated Files

```
test/06_subclustering/
│
├── config.json
├── test_subclustering.py
├── README.md
│
└── results/
    ├── 20251007_153045_subcluster_annotations.csv
    ├── 20251007_153045_subcluster_report.html
    ├── 20251007_153045_subcluster_summary.json
    └── 20251007_153045_test_log.txt
```

### Output CSV Structure

| Result ID | main_cell_type | sub_cell_type | key_markers | reason |
|-----------|----------------|---------------|-------------|---------|
| 1 | Memory CD8+ T cells | Central memory CD8+ T cells | IL7R, CD8A, CD8B, CCL4 | IL7R and CCL4 expression... |
| 2 | Exhausted CD8+ T cells | Terminally exhausted CD8+ T cells | LAYN, HAVCR2, TIGIT | HAVCR2 (TIM-3) and TIGIT... |
| 3 | Effector CD8+ T cells | Cytotoxic CD8+ T cells | GZMK, GZMH, PRF1, NKG7 | High granzyme expression... |
| 4 | Proliferating CD8+ T cells | Activated CD8+ T cells | MKI67, PCNA, TOP2A | MKI67 and PCNA indicate... |

## Workflow

1. **Input Preparation**: Marker data for subclusters (2 columns: ID, Markers)
2. **Prompt Construction**: Create annotation prompt with major cluster context
3. **LLM Analysis**: Identify top 2 cell types per subcluster with reasoning
4. **Result Extraction**: Parse LLM output into structured format
5. **CSV Export**: Save results with all metadata
6. **HTML Generation**: Create formatted report

## Validation Checks

### 1. Output Format

- ✅ CSV file created
- ✅ All expected columns present
- ✅ Correct number of rows (one per subcluster)

### 2. Cell Type Annotations

- ✅ main_cell_type populated for all subclusters
- ✅ sub_cell_type populated for all subclusters
- ✅ Biologically plausible given major cluster context

### 3. Supporting Information

- ✅ key_markers identified
- ✅ Biological reasoning provided
- ✅ HTML report generated

## Troubleshooting

### 1. Missing Data File

**Error**: `Subcluster data file not found`

**Solution**: The test will create synthetic data automatically. To use your own:
```python
# Place subcluster_results.csv in CASSIA/data/
# Format: Two columns (Subcluster, Top_Markers)
```

### 2. Incorrect Output Format

**Issue**: Missing columns in output CSV

**Solutions**:
- Verify LLM response format is parseable
- Check model supports structured output
- Try different model

### 3. Biologically Implausible Results

**Issue**: Annotations don't match major cluster

**Causes**:
- Marker quality poor
- Major cluster info unclear
- Model hallucination

**Solutions**:
```json
{
  "method_params": {
    "major_cluster_info": "CD8 T cells from human intestine",  // More specific
    "n_genes": 100  // More markers
  },
  "temperature": 0  // More deterministic
}
```

### 4. HTML Report Not Generated

**Issue**: CSV created but no HTML

**Cause**: Report generation function may have failed

**Check**: Look for errors in log file

## Advanced Usage

### Multiple Major Clusters

Test different major cluster types:

```python
major_clusters = [
    "CD8 T cells",
    "macrophages",
    "epithelial cells"
]

for cluster in major_clusters:
    config['method_params']['major_cluster_info'] = cluster
    # Run test
```

### Real Subcluster Data

Use actual scRNA-seq subcluster results:

```python
# After running clustering analysis
subclusters = pd.DataFrame({
    'Subcluster': ['0', '1', '2', '3'],
    'Top_Markers': [
        'GENE1, GENE2, GENE3, ...',
        'GENEA, GENEB, GENEC, ...',
        ...
    ]
})
subclusters.to_csv('../../data/subcluster_results.csv', index=False)
```

### Batch Subclustering

Run multiple subclustering analyses:

```python
from CASSIA.subclustering import runCASSIA_n_subcluster

# Run 5 times for UQ
runCASSIA_n_subcluster(
    n=5,
    marker=subcluster_data,
    major_cluster_info="CD8 T cells",
    base_output_name="cd8_subclusters",
    model="google/gemini-2.5-flash-preview",
    provider="openrouter",
    max_workers=5
)
```

## Performance Benchmarks

| Configuration | Time | Notes |
|---------------|------|-------|
| 4 subclusters | 2-3 min | Single LLM call |
| 10 subclusters | 3-5 min | Single LLM call |
| Batch (5 runs × 4 subclusters) | 10-15 min | Parallel execution |

## Related Tests

- **01_runCASSIA_batch** - Identify major clusters first
- **03_annotation_boost** - Deep-dive for ambiguous subclusters
- **05_uncertainty_quantification** - Multiple runs for confidence

## References

- [CASSIA Documentation](../../data/CASSIA_Package_Documentation.md)
- [Subclustering Source](../../subclustering.py:334)
- [Testing Framework Plan](../../../../dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md)

---

**Last Updated**: 2025-10-07
**Version**: 1.0
**Status**: ✅ Active
