# Test: merge_annotations_all

## Overview

**Method**: `merge_annotations_all()`
**Purpose**: Test annotation merging at multiple granularity levels
**Category**: Annotation Post-Processing
**Priority**: Medium
**Expected Runtime**: 3-5 minutes

## Description

This test validates the **annotation merging functionality** which consolidates cell type annotations at three different granularity levels using LLM-guided grouping. This is particularly useful for:

- Creating hierarchical cell type categories
- Standardizing diverse annotation terminology
- Supporting analysis at multiple levels of biological detail
- Comparing annotations across datasets at different resolutions

### Three Merging Levels

1. **Broad (Merged_Grouping_1)**: General lineage categories
   - Examples: "Myeloid cells", "T cells", "B cells"
   - Combines related cell types into major categories
   - Useful for high-level tissue composition analysis

2. **Detailed (Merged_Grouping_2)**: Intermediate specificity
   - Examples: "Macrophages", "CD4 T cells", "CD8 T cells"
   - Balances detail and generality
   - Useful for comparative analysis across studies

3. **Very Detailed (Merged_Grouping_3)**: Normalized specific annotations
   - Examples: "Inflammatory macrophages", "Naive CD4+ T cells", "Memory B cells"
   - Preserves specific states and markers
   - Uses consistent nomenclature
   - Useful for fine-grained analysis

### What This Test Validates

- ✅ Multi-level merging execution (3 parallel LLM calls)
- ✅ All three merged grouping columns created
- ✅ Proper hierarchy (Broad → Detailed → Very Detailed)
- ✅ Biologically sensible groupings at each level
- ✅ Non-null values in merged columns
- ✅ Integration with batch annotation results

## Prerequisites

### Required

- Python 3.8+
- CASSIA package installed
- OpenRouter API key
- Prerequisite batch annotation results
- ~3-5 minutes execution time

### Environment Setup

```bash
# Set API key
export OPENROUTER_API_KEY='your-api-key-here'

# Verify CASSIA is installed
python -c "import CASSIA; print(CASSIA.__version__)"
```

### Data Dependencies

**IMPORTANT**: This test requires batch annotation results with the following columns:
- `True Cell Type` (cluster identifier)
- `Predicted Main Cell Type` (general annotation)
- `Predicted Sub Cell Types` (subtype annotation)

The test will automatically create these prerequisite results if configured to do so.

## Quick Start

```bash
# Navigate to test directory
cd test/04_merge_annotations/

# Run test (will create prerequisites automatically)
python test_merge_annotations.py
```

## Configuration

The test is configured via `config.json`:

```json
{
  "test_name": "test_merge_annotations",
  "description": "Test annotation merging at multiple granularity levels",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.3,
  "data_file": "processed.csv",
  "method_params": {
    "batch_size": 20,
    "additional_context": "human large intestine tissue"
  },
  "prerequisites": {
    "requires_batch_results": true,
    "batch_results_file": "merge_test_batch_full.csv"
  },
  "validation": {
    "check_merged_columns": true,
    "expected_merged_columns": [
      "Merged_Grouping_1",
      "Merged_Grouping_2",
      "Merged_Grouping_3"
    ],
    "check_grouping_hierarchy": true
  }
}
```

### Key Parameters

- **batch_size** (20): Number of clusters processed per LLM call
- **additional_context**: Tissue/species context for better grouping
- **temperature** (0.3): Lower temperature for consistent groupings
- **requires_batch_results**: Auto-create batch results if needed

### Customization

Edit `config.json` to:

- **Use existing batch results** by setting `requires_batch_results: false`
- **Adjust batch size** for API efficiency (10-30 recommended)
- **Add domain context** via `additional_context` for better groupings
- **Change model** to test different LLM behaviors

## Expected Outputs

### Console Output

```
================================================================================
CASSIA Test: merge_annotations_all (Multi-Level Merging)
================================================================================
✓ API key validated for provider: openrouter
✓ Loaded 6 rows from processed.csv
  - Number of clusters: 6

================================================================================
Test configuration:
  Test: test_merge_annotations
  Model: google/gemini-2.5-flash-preview
  Provider: openrouter
================================================================================

================================================================================
CREATING PREREQUISITE BATCH RESULTS
================================================================================
Running batch annotation to create prerequisite data...
Output name: merge_test_batch

=== Starting cell type analysis ===
Analyzing monocyte...
Analyzing plasma cell...
...
✓ Cell type analysis completed

Starting: Prerequisite batch annotation
Completed: Prerequisite batch annotation (135.21s = 2.25min)
✓ Created prerequisite batch results: merge_test_batch_full.csv

================================================================================
RUNNING ANNOTATION MERGING
================================================================================

Merge configuration:
  - Input file: merge_test_batch_full.csv
  - Batch size: 20
  - Context: human large intestine tissue
  - Levels: Broad, Detailed, Very Detailed (3 parallel calls)

Starting: Annotation merging (all levels)

Processing all three detail levels in parallel...
Processed clusters 1-6 out of 6
Completed processing for broad detail level
Processed clusters 1-6 out of 6
Completed processing for detailed detail level
Processed clusters 1-6 out of 6
Completed processing for very_detailed detail level
Combined results saved to merge_test_batch_merged.csv

Completed: Annotation merging (all levels) (125.43s = 2.09min)

✓ Annotation merging completed
  - Result shape: (6, 15)
  - Columns added: ['Merged_Grouping_1', 'Merged_Grouping_2', 'Merged_Grouping_3']

================================================================================
VALIDATING MERGED ANNOTATIONS
================================================================================
Validating merged annotation columns:
  ✓ Found Merged_Grouping_1: 6/6 non-null values
  ✓ Found Merged_Grouping_2: 6/6 non-null values
  ✓ Found Merged_Grouping_3: 6/6 non-null values

Analyzing grouping hierarchy:
  - Merged_Grouping_1: 4 unique groups
  - Merged_Grouping_2: 5 unique groups
  - Merged_Grouping_3: 6 unique groups
  ✓ Hierarchy validation passed: Broad → Detailed → Very Detailed

Example groupings by level:
--------------------------------------------------------------------------------

Cluster: monocyte
  Predicted: Classical Monocyte
  Broad (G1): Myeloid cells
  Detailed (G2): Monocytes
  Very Detailed (G3): Classical monocytes

Cluster: plasma cell
  Predicted: Plasma Cell
  Broad (G1): B cells
  Detailed (G2): Plasma cells
  Very Detailed (G3): Plasma cells

Cluster: cd8_positive_alpha_beta_t_cell
  Predicted: Cytotoxic T Cell
  Broad (G1): T cells
  Detailed (G2): CD8 T cells
  Very Detailed (G3): Cytotoxic CD8+ T cells

================================================================================
ARCHIVING RESULTS
================================================================================
  ✓ Archived merged annotations: 20251007_144523_merged_annotations.csv
  ✓ Cleaned up prerequisite files
  ✓ Saved test summary: 20251007_144523_merge_summary.json

================================================================================
MERGING RESULTS SUMMARY
================================================================================
Clusters processed: 6
Merging levels completed: 3
Unique groups per level:
  - Broad (Grouping_1): 4
  - Detailed (Grouping_2): 5
  - Very Detailed (Grouping_3): 6

Test completed successfully (260.64s = 4.34min)

================================================================================
✅ ANNOTATION MERGING TEST COMPLETED SUCCESSFULLY
📁 Results directory: results/
📄 Log file: results/20251007_144523_test_log.txt
================================================================================
```

### Generated Files

```
test/04_merge_annotations/
│
├── config.json                                    # Test configuration
├── test_merge_annotations.py                      # Test script
├── README.md                                      # This file
│
└── results/                                       # Timestamped results
    ├── 20251007_144523_merged_annotations.csv     # Merged annotations
    ├── 20251007_144523_merge_summary.json        # Test metadata
    └── 20251007_144523_test_log.txt              # Detailed log
```

### Merged Annotations CSV Structure

The output CSV includes original columns plus three new merged grouping columns:

| True Cell Type | Predicted Main Cell Type | Merged_Grouping_1 | Merged_Grouping_2 | Merged_Grouping_3 |
|----------------|--------------------------|-------------------|-------------------|-------------------|
| monocyte | Classical Monocyte | Myeloid cells | Monocytes | Classical monocytes |
| plasma cell | Plasma Cell | B cells | Plasma cells | Plasma cells |
| cd8+ T cell | Cytotoxic T Cell | T cells | CD8 T cells | Cytotoxic CD8+ T cells |
| transit amplifying cell | Transit Amplifying Cell | Epithelial cells | Intestinal epithelial cells | Transit amplifying cells |
| enteroendocrine cell | Enteroendocrine Cell | Epithelial cells | Enteroendocrine cells | Enteroendocrine cells |
| intestinal stem cell | Intestinal Stem Cell | Stem cells | Intestinal stem cells | Lgr5+ intestinal stem cells |

## Merging Workflow

### How It Works

1. **Input Validation**: Verifies required columns exist in batch results
2. **Parallel Processing**: Runs 3 LLM calls simultaneously (one per level)
3. **LLM Grouping**: Each level uses tailored prompts for appropriate granularity
4. **Result Combination**: Merges all three grouping columns into single DataFrame
5. **Output**: Returns/saves DataFrame with original + merged columns

### Prompting Strategy

Each detail level uses a different prompt strategy:

**Broad (Grouping_1)**:
```
Prompt: "Create general cell lineage categories"
Example: "macrophage, inflammatory macrophage" → "Myeloid cells"
Strategy: Maximize consolidation into major lineages
```

**Detailed (Grouping_2)**:
```
Prompt: "Create intermediate-level cell groupings"
Example: "macrophage, inflammatory macrophage" → "Macrophages"
Strategy: Balance specificity and generality
```

**Very Detailed (Grouping_3)**:
```
Prompt: "Normalize and preserve specific annotations"
Example: "macrophage, inflammatory macrophage" → "Inflammatory macrophages"
Strategy: Standardize nomenclature while preserving detail
```

### Parallel Execution

The function uses `ThreadPoolExecutor` with 3 workers to run all three levels concurrently:

```python
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    # Submit all 3 levels at once
    futures = {
        executor.submit(merge_annotations, detail_level="broad"): "broad",
        executor.submit(merge_annotations, detail_level="detailed"): "detailed",
        executor.submit(merge_annotations, detail_level="very_detailed"): "very_detailed"
    }
    # Collect results as they complete
```

This reduces runtime from ~6-9 minutes (sequential) to ~3-5 minutes (parallel).

## Validation Checks

### 1. Column Existence

- ✅ `Merged_Grouping_1` created
- ✅ `Merged_Grouping_2` created
- ✅ `Merged_Grouping_3` created
- ✅ All merged columns have non-null values

### 2. Hierarchy Validation

- ✅ Number of unique groups: Grouping_1 ≤ Grouping_2 ≤ Grouping_3
- ✅ Grouping_1 has fewest unique values (most general)
- ✅ Grouping_3 has most unique values (most specific)

### 3. Biological Validity

- ✅ Groupings make biological sense
- ✅ Related cell types grouped appropriately
- ✅ Terminology is consistent and standardized
- ✅ Functional/activation states preserved where appropriate

## Sample Results

### Merge Summary JSON

```json
{
  "total_clusters": 6,
  "merged_columns": [
    "Merged_Grouping_1",
    "Merged_Grouping_2",
    "Merged_Grouping_3"
  ],
  "unique_groups": {
    "Merged_Grouping_1": 4,
    "Merged_Grouping_2": 5,
    "Merged_Grouping_3": 6
  }
}
```

### Grouping Statistics

| Level | Unique Groups | Granularity | Use Case |
|-------|---------------|-------------|----------|
| Grouping_1 (Broad) | 4 | Lineage-level | Tissue composition |
| Grouping_2 (Detailed) | 5 | Cell type-level | Cross-study comparison |
| Grouping_3 (Very Detailed) | 6 | Subtype-level | Fine-grained analysis |

## Troubleshooting

### 1. Missing Required Columns

**Error**: `Required columns not found: Predicted Main Cell Type`

**Cause**: Batch results missing necessary columns

**Solution**: Ensure batch annotation creates both `Predicted Main Cell Type` and `Predicted Sub Cell Types`:
```python
CASSIA.runCASSIA_batch(...)  # Should output both columns
```

### 2. Hierarchy Validation Failed

**Warning**: `Hierarchy may not follow expected pattern`

**Possible Causes**:
- LLM not following instructions properly
- Insufficient differentiation in prompts
- Small dataset (all cells same type)

**Solutions**:
- Try different model
- Increase temperature slightly (0.3 → 0.5)
- Use larger/more diverse dataset
- Check if warning is spurious (some datasets naturally have similar counts)

### 3. Parallel Processing Errors

**Error**: `All parallel processing tasks failed`

**Solutions**:
- Check API rate limits (3 simultaneous calls)
- Verify API key is valid
- Try sequential processing (modify code to loop instead of parallel)
- Check network connectivity

### 4. Empty Merged Columns

**Issue**: Merged_Grouping columns exist but contain empty/null values

**Cause**: LLM response parsing failure

**Solutions**:
```json
{
  "temperature": 0.1,  // Lower temperature
  "batch_size": 10     // Smaller batches
}
```

Or check logs for LLM response format issues.

### 5. Prerequisite Batch Results Missing

**Error**: `Batch results file not found`

**Solution**:
```json
{
  "prerequisites": {
    "requires_batch_results": true  // Auto-create
  }
}
```

Or run Test 01 first:
```bash
cd ../01_runCASSIA_batch/
python test_batch.py
# Copy results to ../04_merge_annotations/
```

### 6. Inconsistent Groupings

**Issue**: Same cell type gets different groupings across runs

**Cause**: High temperature or model variability

**Solutions**:
```json
{
  "temperature": 0.0,  // Deterministic
  "model": "google/gemini-2.5-flash-preview"  // More consistent model
}
```

## Advanced Usage

### Use Pre-existing Batch Results

Skip prerequisite creation if you already have batch results:

```json
{
  "prerequisites": {
    "requires_batch_results": false,
    "batch_results_file": "/path/to/existing_batch_results_full.csv"
  }
}
```

### Customize Context for Better Groupings

Add domain-specific context:

```json
{
  "method_params": {
    "additional_context": "human pancreatic islets, focus on beta cell subtypes and their maturation states"
  }
}
```

### Adjust Batch Size

For API efficiency or rate limit management:

```json
{
  "method_params": {
    "batch_size": 10  // Smaller (safer) or 30 (faster)
  }
}
```

### Test Single Level Only

Modify test script to test one level:

```python
# In test_merge_annotations.py
from CASSIA.merging_annotation import merge_annotations  # Not merge_annotations_all

merged_df = merge_annotations(
    csv_path=batch_csv,
    detail_level="broad",  # or "detailed" or "very_detailed"
    ...
)
```

### Compare Different Models

Test how different models handle grouping:

```python
models = [
    "google/gemini-2.5-flash-preview",
    "anthropic/claude-3.5-sonnet",
    "openai/gpt-4o"
]

for model in models:
    config['model'] = model
    # Run test and compare results
```

## Performance Benchmarks

### Expected Performance (gemini-2.5-flash-preview)

| Component | Time | Notes |
|-----------|------|-------|
| Prerequisite Batch | 2-3 min | 6 clusters |
| Merging (parallel) | 2-3 min | 3 levels simultaneously |
| Validation | <5 sec | Column checks |
| **Total** | **3-5 min** | **Complete test** |

### Sequential vs Parallel

| Approach | Time | Speed-up |
|----------|------|----------|
| Sequential (3 levels) | 6-9 min | 1.0x |
| Parallel (3 workers) | 3-5 min | 2.0x |

### Resource Usage

- **Memory**: 200-400 MB
- **Disk Space**: ~5-10 MB per test run
- **Network**: ~30-50 KB per LLM call (3 calls total)
- **CPU**: Multi-threaded (3 parallel workers)

### Optimization Tips

1. **Reuse batch results** across multiple merge tests
2. **Increase batch_size** to reduce number of LLM calls (20-30 recommended)
3. **Use faster model** if quality is acceptable
4. **Cache results** for development/testing

## Interpreting Results

### Good Merging Results

Signs of successful merging:
- Clear hierarchy (Broad < Detailed < Very Detailed unique counts)
- Biologically sensible groupings at each level
- Consistent terminology within each level
- Related cell types grouped together appropriately

### Example Good Hierarchy

```
Monocyte cluster:
  Grouping_1: "Myeloid cells"        (broad lineage)
  Grouping_2: "Monocytes"            (specific cell type)
  Grouping_3: "Classical monocytes"  (detailed subtype)

CD4 T cell cluster:
  Grouping_1: "T cells"              (broad lineage)
  Grouping_2: "CD4 T cells"          (specific cell type)
  Grouping_3: "Naive CD4+ T cells"   (detailed subtype)
```

### Poor Merging Results

Warning signs:
- No hierarchy (all levels have same unique counts)
- Inconsistent terminology (e.g., "T cells" and "T lymphocytes" in same column)
- Over-specific Grouping_1 (e.g., each cluster gets unique group)
- Over-general Grouping_3 (e.g., everything grouped into one category)

## Use Cases

### When to Use Merging

1. **Cross-dataset comparison** at standardized granularity
2. **Multi-resolution analysis** of tissue composition
3. **Annotation harmonization** across different naming conventions
4. **Hierarchical clustering** with biological categories
5. **Publication figures** requiring different detail levels

### Downstream Analysis Examples

**Broad Level (Grouping_1)**:
```python
# Tissue composition pie chart
composition = merged_df['Merged_Grouping_1'].value_counts()
plt.pie(composition, labels=composition.index)
```

**Detailed Level (Grouping_2)**:
```python
# Cross-study comparison
study1_groups = study1_df['Merged_Grouping_2'].value_counts(normalize=True)
study2_groups = study2_df['Merged_Grouping_2'].value_counts(normalize=True)
correlation = pearson(study1_groups, study2_groups)
```

**Very Detailed Level (Grouping_3)**:
```python
# Fine-grained differential abundance
abundance_matrix = create_abundance_matrix(merged_df, level='Merged_Grouping_3')
differential_test(abundance_matrix, conditions=['healthy', 'disease'])
```

## Related Tests

- **01_runCASSIA_batch** - Creates prerequisite batch results
- **02_runCASSIA_pipeline** - Uses merging as integrated stage
- **07_celltype_comparison** - Compares merged annotations across datasets
- **11_report_generation** - Visualizes merged annotations

## References

- [CASSIA Documentation](../../data/CASSIA_Package_Documentation.md)
- [merge_annotations_all Source](../../merging_annotation.py:299)
- [Testing Framework Plan](../../../../dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md)
- [Cell Type Ontology](http://www.obofoundry.org/ontology/cl.html) - Standard nomenclature

---

**Last Updated**: 2025-10-07
**Version**: 1.0
**Status**: ✅ Active
