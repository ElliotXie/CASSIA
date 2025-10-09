# Test: runCASSIA_batch

## Overview

**Method**: `CASSIA.runCASSIA_batch()`
**Purpose**: Test basic batch annotation functionality
**Category**: Core Annotation
**Priority**: High
**Expected Runtime**: 3-5 minutes

## Description

This test validates the core batch annotation functionality of CASSIA. It processes multiple cell clusters in parallel, annotating each one based on marker gene expression patterns. The test uses the `processed.csv` dataset which contains 6 clean, well-characterized cell types from human large intestine tissue.

### What This Test Validates

- ✅ Basic batch processing workflow
- ✅ Parallel annotation of multiple clusters
- ✅ Output file generation (full + summary CSVs)
- ✅ Data format compliance
- ✅ Error handling and retry mechanisms
- ✅ Model integration (gemini-2.5-flash-preview via OpenRouter)

## Prerequisites

### Required

- Python 3.8+
- CASSIA package installed
- OpenRouter API key

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
cd test/01_runCASSIA_batch/

# Run test
python test_batch.py
```

## Configuration

The test is configured via `config.json`:

```json
{
  "test_name": "test_runCASSIA_batch",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "processed.csv",
  "method_params": {
    "n_genes": 50,
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4,
    "max_retries": 1
  }
}
```

### Customization

Edit `config.json` to customize:

- **model**: Change LLM model
- **temperature**: Adjust sampling randomness (0-1)
- **n_genes**: Number of top marker genes per cluster
- **max_workers**: Parallel processing threads
- **data_file**: Use different sample dataset

## Expected Outputs

### Console Output

```
================================================================================
Starting test: test_runCASSIA_batch
Model: google/gemini-2.5-flash-preview
Provider: openrouter
Data: processed.csv
================================================================================
✓ API key validated for provider: openrouter
✓ Loaded 6 rows from processed.csv
  - Number of clusters: 6
  - Cluster names: ['monocyte', 'plasma cell', 'cd8-positive, alpha-beta t cell', ...]

Starting CASSIA batch annotation...
Analyzing monocyte...
Analyzing plasma cell...
...
Completed: CASSIA batch annotation (125.34s = 2.09min)

Validating output format...
✓ Validation passed

Saving timestamped results...
✓ Saved full results: results/20251007_143022_batch_full.csv
✓ Saved summary results: results/20251007_143022_batch_summary.csv

================================================================================
TEST RESULTS SUMMARY
================================================================================
Total clusters processed: 6
Successful annotations: 6
Average iterations: 1.00

Sample annotations:
  monocyte → Classical Monocyte
  plasma cell → Plasma Cell
  cd8-positive, alpha-beta t cell → Cytotoxic T Cell

================================================================================
✅ TEST COMPLETED SUCCESSFULLY
📁 Results saved to: results/
📄 Log file: results/20251007_143022_test_log.txt
================================================================================
```

### Generated Files

#### Results Directory (`results/`)

1. **`[timestamp]_batch_full.csv`** - Complete annotation results
   ```csv
   True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types,Possible Mixed Cell Types,Marker Number,Marker List,Iterations,Model,Provider,Tissue,Species,Additional Info,Conversation History
   monocyte,Classical Monocyte,"CD14+ Monocyte, CD16- Monocyte",,50,"NKAIN3,SPP1,CDH19,...",1,google/gemini-2.5-flash-preview,openrouter,large intestine,human,...
   ```

2. **`[timestamp]_batch_summary.csv`** - Condensed results
   ```csv
   True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types,Possible Mixed Cell Types,Marker List,Iterations,Model,Provider,Tissue,Species
   monocyte,Classical Monocyte,"CD14+ Monocyte, CD16- Monocyte",,NKAIN3,SPP1,...,1,google/gemini-2.5-flash-preview,openrouter,large intestine,human
   ```

3. **`[timestamp]_test_log.txt`** - Detailed execution log
   ```
   2025-10-07 14:30:22 - test_runCASSIA_batch - INFO - Starting test: test_runCASSIA_batch
   2025-10-07 14:30:22 - test_runCASSIA_batch - INFO - Model: google/gemini-2.5-flash-preview
   ...
   ```

## Validation Checks

The test performs the following validations:

### 1. Output Format
- ✅ Required columns present: `True Cell Type`, `Predicted Main Cell Type`, `Predicted Sub Cell Types`
- ✅ Correct data types
- ✅ No missing column headers

### 2. Data Completeness
- ✅ Minimum 6 rows (one per cluster)
- ✅ No null values in `Predicted Main Cell Type` column
- ✅ All clusters processed successfully

### 3. Execution Metrics
- ✅ Completes within expected time (3-5 minutes)
- ✅ No API errors or authentication issues
- ✅ Proper error handling and retries

## Sample Results

### Input Data (processed.csv)

| Cluster | Cell Type | Top Markers |
|---------|-----------|-------------|
| 1 | monocyte | NKAIN3, SPP1, CDH19, ... |
| 2 | plasma cell | IGLL5, IGLV6-57, JCHAIN, ... |
| 3 | cd8+ T cell | KLRC2, XCL2, LINC02446, ... |
| 4 | transit amplifying cell | HJURP, UBE2C, CENPA, ... |
| 5 | enteroendocrine cell | GCG, PCSK1, AC090679.1, ... |
| 6 | intestinal crypt stem cell | AL031284.1, AC012085.1, IGFBPL1, ... |

### Expected Annotations

| True Cell Type | Predicted Main Cell Type | Confidence |
|----------------|--------------------------|------------|
| monocyte | Classical Monocyte | High |
| plasma cell | Plasma Cell | High |
| cd8-positive, alpha-beta t cell | Cytotoxic T Cell | High |
| transit amplifying cell | Transit Amplifying Cell | High |
| intestinal enteroendocrine cell | Enteroendocrine Cell | High |
| intestinal crypt stem cell | Intestinal Stem Cell | High |

## Troubleshooting

### API Key Error

**Error**:
```
EnvironmentError: OPENROUTER_API_KEY not found in environment
```

**Solution**:
```bash
export OPENROUTER_API_KEY='your-api-key-here'
```

Or set it programmatically:
```python
import os
os.environ['OPENROUTER_API_KEY'] = 'your-key'
```

### Data File Not Found

**Error**:
```
FileNotFoundError: Data file not found: .../data/processed.csv
```

**Solution**:
- Ensure you're running from the test directory: `cd test/01_runCASSIA_batch/`
- Check that `../../data/processed.csv` exists
- Verify project structure is intact

### Import Errors

**Error**:
```
ModuleNotFoundError: No module named 'CASSIA'
```

**Solution**:
```bash
# Install CASSIA in development mode
cd CASSIA_python/
pip install -e .
```

### Slow Performance

**Issue**: Test taking longer than 5 minutes

**Solutions**:
1. Reduce `max_workers` in config (less parallel, but more stable)
2. Reduce `n_genes` (fewer markers to analyze)
3. Check network connectivity to OpenRouter API
4. Try a different time (API may be busy)

### Authentication Errors (401)

**Error**:
```
401 Client Error: Unauthorized
```

**Solutions**:
1. Verify API key is correct
2. Check API key has sufficient credits
3. Ensure OpenRouter account is active
4. Test API key with curl:
   ```bash
   curl -H "Authorization: Bearer $OPENROUTER_API_KEY" \
        https://openrouter.ai/api/v1/models
   ```

### Empty or Incomplete Results

**Issue**: Some clusters not annotated

**Possible Causes**:
1. API rate limiting - add delays between requests
2. Network timeout - increase timeout in config
3. Model failure - check model status at OpenRouter

**Solution**:
```json
{
  "method_params": {
    "max_retries": 3,
    "max_workers": 2
  }
}
```

## Advanced Usage

### Run with Different Model

```json
{
  "model": "anthropic/claude-3.5-sonnet",
  "provider": "openrouter"
}
```

### Use Subset of Data

```python
# Modify test_batch.py
from sample_data import SampleDataLoader

loader = SampleDataLoader()
marker_data = loader.load_subset(
    dataset="processed",
    clusters=["monocyte", "plasma cell"],  # Only test 2 clusters
    n_clusters=None
)
```

### Custom Validation

```json
{
  "validation": {
    "check_output_format": true,
    "expected_columns": ["True Cell Type", "Predicted Main Cell Type"],
    "min_rows": 2,
    "max_nulls": {
      "Predicted Main Cell Type": 0.0,
      "Predicted Sub Cell Types": 0.2
    }
  }
}
```

## Performance Benchmarks

### Expected Performance (gemini-2.5-flash-preview)

| Metric | Expected | Typical |
|--------|----------|---------|
| Total Runtime | 3-5 min | 3.5 min |
| Per-cluster Time | 30-50 sec | 35 sec |
| API Calls | 6-12 | 6 |
| Retries | 0-2 | 0 |

### Resource Usage

- **Memory**: ~200-500 MB
- **CPU**: Multi-threaded (4 workers)
- **Network**: ~10-50 KB per request

## Related Tests

- **02_runCASSIA_pipeline** - Full end-to-end pipeline including batch annotation
- **03_annotation_boost** - Deep-dive analysis for ambiguous clusters
- **05_uncertainty_quantification** - Multiple runs for annotation stability

## References

- [CASSIA Documentation](../../data/CASSIA_Package_Documentation.md)
- [runCASSIA_batch Source](../../tools_function.py:367)
- [Testing Framework Plan](../../../../dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md)
- [Master Test README](../README.md)

---

**Last Updated**: 2025-10-07
**Version**: 1.0
**Status**: ✅ Active
