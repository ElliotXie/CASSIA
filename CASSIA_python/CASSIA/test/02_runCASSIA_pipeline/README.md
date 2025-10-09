# Test: runCASSIA_pipeline

## Overview

**Method**: `CASSIA.runCASSIA_pipeline()`
**Purpose**: Test full end-to-end CASSIA pipeline
**Category**: Core Annotation
**Priority**: High
**Expected Runtime**: 10-15 minutes

## Description

This test validates the **complete CASSIA pipeline** from raw marker data to final annotated results with quality scoring, annotation boosting for uncertain clusters, and merged annotations. It's the most comprehensive test in the framework, exercising all major CASSIA components in an integrated workflow.

### Pipeline Stages Tested

1. **Batch Annotation** - Parallel annotation of all clusters
2. **Quality Scoring** - LLM-based quality assessment of annotations
3. **Annotation Boost** - Deep-dive analysis for low-scoring clusters
4. **Annotation Merging** - Consolidation at different granularity levels
5. **Report Generation** - HTML reports with visualizations

### What This Test Validates

- ✅ End-to-end pipeline execution
- ✅ Stage integration (annotation → scoring → boost → merge)
- ✅ Organized output directory structure
- ✅ Quality-based conditional boost logic
- ✅ HTML report generation
- ✅ Final results consolidation
- ✅ Error handling across stages
- ✅ Multi-model coordination

## Prerequisites

### Required

- Python 3.8+
- CASSIA package installed
- OpenRouter API key
- ~15 minutes execution time
- ~500 MB disk space for results

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
cd test/02_runCASSIA_pipeline/

# Run test
python test_pipeline.py
```

## Configuration

The test is configured via `config.json`:

```json
{
  "test_name": "test_runCASSIA_pipeline",
  "annotation_model": "google/gemini-2.5-flash-preview",
  "score_model": "google/gemini-2.5-flash-preview",
  "annotationboost_model": "google/gemini-2.5-flash-preview",
  "merge_model": "google/gemini-2.5-flash-preview",
  "annotation_provider": "openrouter",
  "method_params": {
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4,
    "score_threshold": 75,
    "merge_annotations": true
  }
}
```

### Key Parameters

- **score_threshold** (75): Clusters scoring below this trigger annotation boost
- **merge_annotations** (true): Enable final annotation merging
- **max_workers** (4): Parallel processing threads
- **conversation_history_mode** ("final"): History tracking for boost analysis

### Customization

Edit `config.json` to:

- **Use different models** per stage (annotation, scoring, boost, merge)
- **Adjust score threshold** to control boost sensitivity
- **Disable merging** by setting `merge_annotations: false`
- **Change tissue/species** for different contexts

## Expected Outputs

### Console Output

```
================================================================================
CASSIA Test: runCASSIA_pipeline (Full End-to-End)
================================================================================
✓ API key validated for provider: openrouter
✓ Loaded 6 rows from processed.csv
  - Number of clusters: 6

================================================================================
STARTING FULL CASSIA PIPELINE
================================================================================

Pipeline stages:
  1. Batch annotation
  2. Quality scoring
  3. Annotation boost (for low-scoring clusters)
  4. Annotation merging
  5. Report generation

=== Starting cell type analysis ===
Analyzing monocyte...
Analyzing plasma cell...
...
✓ Cell type analysis completed

=== Starting scoring process ===
Processed row 1: Score = 85
Processed row 2: Score = 92
...
✓ Scoring process completed

=== Creating final combined results ===
✓ Final combined results saved

=== Generating main reports ===
✓ Main reports generated

=== Analyzing low-scoring clusters ===
Found 1 clusters with scores below 75:
['transit amplifying cell of large intestine']

=== Starting boost annotation for low-scoring clusters ===
Processing low score cluster: transit amplifying cell
...
✓ Boost annotation completed

Completed: Complete CASSIA pipeline (645.23s = 10.75min)

================================================================================
VALIDATING PIPELINE OUTPUT
================================================================================

✓ Found pipeline output directory: CASSIA_large_intestine_human_20251007_153045

Validating directory structure:
  ✓ Found directory: 01_annotation_results
  ✓ Found directory: 02_reports
  ✓ Found directory: 03_boost_analysis
  ✓ Found final results: pipeline_test_results_FINAL_RESULTS.csv

Validating final results:
  ✓ Final results shape: (6, 15)
  ✓ All expected columns present
  ✓ Score range: 72.0 - 95.0

Validating HTML reports:
  ✓ Found 7 HTML report(s)
    - report_monocyte.html
    - report_plasma_cell.html
    - ...
    - pipeline_test_results.html

Validating boost analysis:
  ✓ Found boost analysis for 1 cluster(s)
    - transit_amplifying_cell_of_large_intestine

Archiving results to test directory:
  ✓ Archived final results: results/20251007_153045_pipeline_final.csv
  ✓ Saved pipeline summary: results/20251007_153045_pipeline_summary.json

================================================================================
PIPELINE RESULTS SUMMARY
================================================================================
Output directory: CASSIA_large_intestine_human_20251007_153045
Total clusters processed: 6
Average quality score: 84.17
Clusters below threshold (75): 1
✓ Merged annotations generated

Sample pipeline results:
  True Cell Type: monocyte | Predicted Main Cell Type: Classical Monocyte | Score: 85 | Merged_Annotation_Broad: Monocyte
  True Cell Type: plasma cell | Predicted Main Cell Type: Plasma Cell | Score: 92 | Merged_Annotation_Broad: Plasma Cell

================================================================================
✅ PIPELINE TEST COMPLETED SUCCESSFULLY
📁 Pipeline output: CASSIA_large_intestine_human_20251007_153045
📁 Test results: results/
📄 Log file: results/20251007_153045_test_log.txt
================================================================================
```

### Generated Directory Structure

```
CASSIA_large_intestine_human_[timestamp]/
│
├── 01_annotation_results/
│   ├── pipeline_test_results_FINAL_RESULTS.csv    # Main output
│   └── intermediate_files/
│       ├── pipeline_test_results_full.csv
│       ├── pipeline_test_results_summary.csv
│       ├── pipeline_test_results_scored.csv
│       └── pipeline_test_results_merged.csv
│
├── 02_reports/
│   ├── pipeline_test_results.html                 # Index page
│   ├── report_monocyte.html
│   ├── report_plasma_cell.html
│   ├── report_cd8_positive_alpha_beta_t_cell.html
│   ├── report_transit_amplifying_cell.html
│   ├── report_intestinal_enteroendocrine_cell.html
│   └── report_intestinal_crypt_stem_cell.html
│
└── 03_boost_analysis/
    └── transit_amplifying_cell_of_large_intestine/
        ├── pipeline_test_results_transit_amplifying_cell_boosted_conversation.json
        ├── pipeline_test_results_transit_amplifying_cell_boosted_summary.html
        └── pipeline_test_results_transit_amplifying_cell_boosted_raw.txt
```

### Test Results Directory (`results/`)

```
results/
├── 20251007_153045_pipeline_final.csv        # Archived final results
├── 20251007_153045_pipeline_summary.json     # Pipeline metadata
└── 20251007_153045_test_log.txt             # Detailed test log
```

## Pipeline Workflow

### Stage 1: Batch Annotation

```python
# Annotates all clusters in parallel
runCASSIA_batch(
    marker=data,
    model="google/gemini-2.5-flash-preview",
    provider="openrouter",
    tissue="large intestine",
    species="human"
)
```

**Output**: `*_full.csv`, `*_summary.csv`

### Stage 2: Quality Scoring

```python
# Scores each annotation for quality
runCASSIA_score_batch(
    input_file="*_full.csv",
    model="google/gemini-2.5-flash-preview",
    provider="openrouter"
)
```

**Output**: `*_scored.csv` with Score column (0-100)

### Stage 3: Annotation Boost (Conditional)

```python
# Only runs for clusters with Score < threshold
if cluster_score < 75:
    runCASSIA_annotationboost(
        cluster_name=cluster,
        num_iterations=5,
        model="google/gemini-2.5-flash-preview"
    )
```

**Output**: Boost analysis folder with conversation history

### Stage 4: Annotation Merging

```python
# Merges annotations at different granularities
merge_annotations_all(
    csv_path="*_scored.csv",
    model="google/gemini-2.5-flash-preview",
    provider="openrouter"
)
```

**Output**: `*_merged.csv` with Merged_Annotation columns

### Stage 5: Report Generation

```python
# Generates HTML reports
runCASSIA_generate_score_report(
    csv_path="*_scored.csv",
    index_name="results_summary"
)
```

**Output**: HTML reports in `02_reports/`

## Validation Checks

### 1. Directory Structure
- ✅ `01_annotation_results/` exists with CSV files
- ✅ `02_reports/` exists with HTML files
- ✅ `03_boost_analysis/` exists if low-scoring clusters present
- ✅ `intermediate_files/` subfolder created

### 2. Final Results File
- ✅ `*_FINAL_RESULTS.csv` exists
- ✅ Contains all expected columns
- ✅ Row count matches input clusters
- ✅ No null values in critical columns

### 3. Quality Scores
- ✅ Score column present
- ✅ Scores in valid range [0-100]
- ✅ Average score calculated
- ✅ Low-scoring clusters identified

### 4. Boost Analysis
- ✅ Triggered for scores below threshold
- ✅ Conversation history saved
- ✅ HTML summary generated
- ✅ Raw text conversation saved

### 5. Reports
- ✅ HTML index page created
- ✅ Individual cluster reports generated
- ✅ Reports contain formatted results
- ✅ Styling and structure correct

## Sample Results

### Final Results CSV

| True Cell Type | Predicted Main Cell Type | Score | Merged_Annotation_Broad | Merged_Annotation_Granular |
|----------------|--------------------------|-------|-------------------------|----------------------------|
| monocyte | Classical Monocyte | 85 | Monocyte | Classical Monocyte |
| plasma cell | Plasma Cell | 92 | Plasma Cell | IgA+ Plasma Cell |
| cd8+ T cell | Cytotoxic T Cell | 88 | T Cell | CD8+ Cytotoxic T Cell |
| transit amplifying cell | Transit Amplifying Cell | 72 | Epithelial Cell | Transit Amplifying Cell |
| enteroendocrine cell | Enteroendocrine Cell | 90 | Enteroendocrine Cell | GLP-1+ Enteroendocrine Cell |
| intestinal stem cell | Intestinal Stem Cell | 95 | Stem Cell | Lgr5+ Intestinal Stem Cell |

### Pipeline Summary JSON

```json
{
  "output_directory": "CASSIA_large_intestine_human_20251007_153045",
  "final_results_file": "...FINAL_RESULTS.csv",
  "total_clusters": 6,
  "score_range": "72.0 - 95.0",
  "html_reports": 7,
  "boost_clusters": 1
}
```

## Troubleshooting

### 1. Pipeline Timeout

**Issue**: Test takes longer than expected

**Solutions**:
```json
{
  "method_params": {
    "max_workers": 2,        // Reduce parallelism
    "score_threshold": 90     // Reduce boost triggers
  }
}
```

### 2. Output Directory Not Found

**Error**: `Pipeline output directory not found`

**Cause**: Directory name mismatch or creation failure

**Solution**:
```bash
# Check for any CASSIA_* directories
ls -d CASSIA_*

# Manually specify output name without spaces
```

### 3. Boost Analysis Fails

**Error**: `Error in pipeline processing cluster X`

**Solutions**:
- Check API rate limits
- Verify marker data quality
- Check conversation_history_mode setting
- Reduce num_iterations in boost

### 4. Merging Fails

**Error**: `Error during annotation merging`

**Solutions**:
```json
{
  "method_params": {
    "merge_annotations": false  // Disable merging
  }
}
```

Or check merge_model/provider settings

### 5. Missing HTML Reports

**Issue**: No HTML files in reports/

**Solutions**:
- Check scoring completed successfully
- Verify write permissions
- Check for errors in report generation stage

### 6. Low Scores Across All Clusters

**Issue**: All clusters scoring below threshold

**Possible Causes**:
- Poor marker quality
- Wrong tissue/species context
- Scoring model too strict

**Solutions**:
- Use better marker data (unprocessed.csv)
- Adjust tissue/species parameters
- Lower score_threshold
- Try different scoring model

## Advanced Usage

### Use Different Models Per Stage

```json
{
  "annotation_model": "google/gemini-2.5-flash-preview",
  "score_model": "anthropic/claude-3.5-sonnet",
  "annotationboost_model": "google/gemini-2.5-flash-preview",
  "merge_model": "openai/gpt-4o"
}
```

### Adjust Boost Sensitivity

```json
{
  "method_params": {
    "score_threshold": 85,  // Higher = more boost triggers
    "conversation_history_mode": "full",  // More context
    "report_style": "total_summary"  // Different report format
  }
}
```

### Run on Larger Dataset

```json
{
  "data_file": "unprocessed.csv",  // Larger dataset
  "method_params": {
    "max_workers": 8,  // More parallelism
    "score_threshold": 80
  }
}
```

### Disable Specific Stages

```python
# In test_pipeline.py, modify pipeline call:
CASSIA.runCASSIA_pipeline(
    marker=marker_data,
    merge_annotations=False,  # Skip merging
    score_threshold=100,      # Skip boost (nothing below 100)
    ...
)
```

## Performance Benchmarks

### Expected Performance (gemini-2.5-flash-preview)

| Stage | Time | Operations |
|-------|------|------------|
| Batch Annotation | 3-5 min | 6 clusters |
| Quality Scoring | 2-3 min | 6 evaluations |
| Annotation Boost | 3-5 min | 1-2 clusters (if triggered) |
| Annotation Merging | 1-2 min | 2 merge levels |
| Report Generation | <1 min | 6+ HTML files |
| **Total** | **10-15 min** | **Full pipeline** |

### Resource Usage

- **Memory**: 500 MB - 1 GB
- **Disk Space**: ~50-100 MB per run
- **Network**: ~100-200 KB per API call
- **CPU**: Multi-threaded (4 workers default)

### Optimization Tips

1. **Reduce max_workers** for stability (fewer parallel calls)
2. **Increase score_threshold** to skip boost for more clusters
3. **Disable merging** if not needed
4. **Use faster models** for non-critical stages
5. **Pre-process markers** to reduce n_genes

## Related Tests

- **01_runCASSIA_batch** - Just the annotation stage
- **03_annotation_boost** - Deep-dive on boost analysis
- **04_merge_annotations** - Just the merging stage
- **05_uncertainty_quantification** - Multiple pipeline runs for UQ

## References

- [CASSIA Documentation](../../data/CASSIA_Package_Documentation.md)
- [runCASSIA_pipeline Source](../../tools_function.py:1602)
- [Testing Framework Plan](../../../../dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md)
- [Test 01 README](../01_runCASSIA_batch/README.md) - Simpler batch-only test

---

**Last Updated**: 2025-10-07
**Version**: 1.0
**Status**: ✅ Active
