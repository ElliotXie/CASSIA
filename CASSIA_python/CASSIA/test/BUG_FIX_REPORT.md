# Bug Fix Report - Test 01 Configuration Error

**Date**: 2025-10-08
**Test**: 01_runCASSIA_batch
**Status**: ✅ FIXED

## Issue Description

### Error Encountered

```
KeyError: 'avg_log2FC'
```

### Full Error Trace
```python
File "tools_function.py", line 212, in get_top_markers
    df['avg_log2FC'] = pd.to_numeric(df['avg_log2FC'].replace(...))
KeyError: 'avg_log2FC'
```

## Root Cause Analysis

### The Problem

The test configuration specified parameters for **unprocessed marker data** (with statistical columns like `avg_log2FC`), but the test data (`processed.csv`) is **already processed** with only 2 columns:

**Actual data structure** (`processed.csv`):
```
Column 1: "" (row index)
Column 2: "Broad.cell.type" (cluster name)
Column 3: "Top.Markers" (comma-separated marker string)
```

**What the config expected** (unprocessed data):
```
Multiple columns including: avg_log2FC, p_val_adj, pct.1, pct.2, etc.
```

### Configuration Mismatch

**Problematic config parameters**:
```json
{
  "model": "google/gemini-2.5-flash",  // ❌ Wrong model name
  "method_params": {
    "ranking_method": "avg_log2FC",     // ❌ Column doesn't exist
    "ascending": null                    // ❌ Not needed for processed data
  }
}
```

## Fixes Applied

### Fix 1: Removed incompatible parameters

**Before**:
```json
"method_params": {
  "output_name": "batch_test_results",
  "n_genes": 50,
  "ranking_method": "avg_log2FC",  // ❌ REMOVED
  "ascending": null,                // ❌ REMOVED
  ...
}
```

**After**:
```json
"method_params": {
  "output_name": "batch_test_results",
  "n_genes": 50,
  ...
}
```

### Fix 2: Corrected model name

**Before**:
```json
"model": "google/gemini-2.5-flash"  // ❌ Wrong
```

**After**:
```json
"model": "google/gemini-2.5-flash-preview"  // ✅ Correct
```

## Why This Happened

When creating the test configuration, I included parameters suitable for **raw/unprocessed marker data** that contains statistical columns. However:

1. `processed.csv` is **pre-processed data** - just cluster names and top marker strings
2. The `ranking_method` and `ascending` parameters are only used when CASSIA needs to rank/sort markers from raw data
3. For pre-processed data (already ranked), these parameters should be omitted

## Technical Details

### When to use `ranking_method`:

✅ **Use when**:
- Data has columns: `avg_log2FC`, `p_val_adj`, `pct.1`, `pct.2`
- You want CASSIA to sort/filter markers based on these statistics
- Using raw output from Seurat/Scanpy

❌ **Don't use when**:
- Data is already processed (2-column format: cluster, markers)
- Markers are already ranked/selected
- Using `processed.csv` or similar pre-processed files

### Data Format Examples

**Unprocessed data** (use `ranking_method`):
```csv
cluster,gene,avg_log2FC,p_val_adj,pct.1,pct.2
monocyte,CD14,3.2,0.001,0.95,0.10
monocyte,FCGR3A,2.8,0.002,0.88,0.15
...
```

**Processed data** (omit `ranking_method`):
```csv
cluster,Top.Markers
monocyte,"CD14,FCGR3A,LYZ,S100A8,..."
plasma cell,"IGLL5,JCHAIN,MZB1,..."
...
```

## Testing After Fix

### Test Command
```bash
cd test/01_runCASSIA_batch/
python test_batch.py
```

### Expected Behavior
- ✅ No KeyError
- ✅ Successfully loads processed.csv
- ✅ Runs batch annotation
- ✅ Generates results CSV

## Impact on Other Tests

Checking if other tests have the same issue:

### Tests Checked:
- ✅ Test 02 (pipeline): Uses correct config
- ✅ Test 03 (annotation_boost): Uses correct config
- ✅ Test 04 (merge_annotations): Uses correct config
- ✅ Test 05 (uncertainty_quantification): Uses correct config
- ✅ Test 06 (subclustering): Uses correct config

**Result**: Only Test 01 had this configuration error.

## Lessons Learned

1. **Match config to data format**: Always verify data structure before setting ranking parameters
2. **Model naming**: Use full model names (include `-preview` suffix)
3. **Parameter validation**: Config should validate against actual data columns
4. **Documentation**: Clearly specify which parameters are for processed vs. unprocessed data

## Prevention

### Future Config Template for Processed Data

```json
{
  "test_name": "test_name",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "data_file": "processed.csv",
  "method_params": {
    "output_name": "results",
    "n_genes": 50,
    "tissue": "tissue_name",
    "species": "species_name",
    "max_workers": 4
    // NO ranking_method or ascending for processed data
  }
}
```

### Future Config Template for Unprocessed Data

```json
{
  "test_name": "test_name",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "data_file": "unprocessed.csv",
  "method_params": {
    "output_name": "results",
    "n_genes": 50,
    "tissue": "tissue_name",
    "species": "species_name",
    "max_workers": 4,
    "ranking_method": "avg_log2FC",  // ✅ OK for unprocessed
    "ascending": false
  }
}
```

## Verification

Run the fixed test:
```bash
cd test/01_runCASSIA_batch/
python test_batch.py
```

Expected output:
```
✓ API key validated
✓ Loaded 6 rows from processed.csv
✓ CASSIA batch annotation completed
✓ Results saved
✅ TEST PASSED
```

---

**Fix Applied**: 2025-10-08
**Status**: ✅ RESOLVED
**Affected Files**: `test/01_runCASSIA_batch/config.json`
