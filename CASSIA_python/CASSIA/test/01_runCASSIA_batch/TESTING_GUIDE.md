# Test 01: runCASSIA_batch - Complete Testing Guide

## Overview

Test 01 now supports **both data formats** commonly used with CASSIA:

1. **Processed Data** - Simple 2-column format (cluster, markers)
2. **Unprocessed Data** - Full Seurat FindAllMarkers output with statistics

## Quick Start

### Test Processed Data Only (Default)
```bash
cd test/01_runCASSIA_batch/
python test_batch.py
```
**Runtime**: 3-5 minutes
**Data**: `processed.csv` (6 clusters, pre-selected markers)

### Test Unprocessed Data Only
```bash
python test_batch.py  # Use config_unprocessed.json manually
```

### Test BOTH Formats
```bash
python test_both.py
```
**Runtime**: 6-10 minutes
**Tests**: Both processed and unprocessed data formats

## Data Format Comparison

### Format 1: Processed Data

**File**: `data/processed.csv`

**Structure**:
```csv
"","Broad.cell.type","Top.Markers"
"1","monocyte","NKAIN3,SPP1,CDH19,TSPAN11,PLP1,..."
"2","plasma cell","IGLL5,IGLV6-57,JCHAIN,FAM92B,..."
```

**Characteristics**:
- ✅ Simple 2-column format
- ✅ Markers already ranked/selected
- ✅ Fast to process
- ✅ Good for quick testing

**Config Parameters**:
```json
{
  "data_file": "processed.csv",
  "method_params": {
    "n_genes": 50
    // NO ranking_method or ascending
  }
}
```

### Format 2: Unprocessed Data

**File**: `data/unprocessed.csv`

**Structure**:
```csv
"","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene"
"CDH19",0,8.866,0.996,0.014,0,"monocyte","CDH19"
"NRXN1",0,8.417,0.974,0.013,0,"monocyte","NRXN1"
```

**Characteristics**:
- ✅ Full Seurat/Scanpy output
- ✅ Contains statistical information
- ✅ Allows custom ranking/filtering
- ✅ Realistic production scenario

**Config Parameters**:
```json
{
  "data_file": "unprocessed.csv",
  "method_params": {
    "n_genes": 50,
    "ranking_method": "avg_log2FC",  // Use log2FC for ranking
    "ascending": false               // Highest values first
  }
}
```

## Available Test Scripts

### 1. test_batch.py (Standard Test)

**Purpose**: Test with processed data (default config)

**Usage**:
```bash
python test_batch.py
```

**Config**: Uses `config.json` (processed data)

**What it tests**:
- Loading processed.csv
- Basic batch annotation
- Result validation
- Output archiving

**Expected output**:
```
✓ API key validated
✓ Loaded 6 rows from processed.csv
✓ CASSIA batch annotation completed (120.5s)
✓ Results validated
✓ Archived results
✅ TEST PASSED
```

### 2. test_both.py (Comprehensive Test)

**Purpose**: Test BOTH data formats in one run

**Usage**:
```bash
python test_both.py
```

**Configs**: Uses both `config.json` AND `config_unprocessed.json`

**What it tests**:
- Processed data (2-column format)
- Unprocessed data (FindAllMarkers format)
- Config parameter differences
- Both ranking methods

**Expected output**:
```
================================================================================
TEST 1/2: PROCESSED DATA (2-column format)
================================================================================
✓ Loaded 6 rows from processed.csv
✓ CASSIA batch (processed) completed (125.3s)
✓ All expected columns present
✓ Archived: processed_results.csv

================================================================================
TEST 2/2: UNPROCESSED DATA (FindAllMarkers format)
================================================================================
✓ Loaded 3700 rows from unprocessed.csv
✓ Found avg_log2FC column (unprocessed data)
✓ CASSIA batch (unprocessed) completed (132.1s)
✓ All expected columns present
✓ Archived: unprocessed_results.csv

================================================================================
FINAL SUMMARY
================================================================================
Total tests: 2
  ✓ Passed: 2
  ✗ Failed: 0

✅ ALL TESTS PASSED
```

## Configuration Files

### config.json (Processed Data)

```json
{
  "test_name": "test_runCASSIA_batch",
  "description": "Test basic batch annotation functionality using sample data",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "processed.csv",
  "method_params": {
    "output_name": "batch_test_results",
    "n_genes": 50,
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4,
    "max_retries": 1,
    "validator_involvement": "v1"
    // Note: NO ranking_method for processed data
  }
}
```

### config_unprocessed.json (Unprocessed Data)

```json
{
  "test_name": "test_runCASSIA_batch_unprocessed",
  "description": "Test batch annotation with unprocessed FindAllMarkers output",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "unprocessed.csv",
  "method_params": {
    "output_name": "batch_test_unprocessed_results",
    "n_genes": 50,
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4,
    "max_retries": 1,
    "ranking_method": "avg_log2FC",  // ✅ Required for unprocessed
    "ascending": false,               // ✅ Required for unprocessed
    "validator_involvement": "v1"
  }
}
```

## Key Differences

| Feature | Processed Data | Unprocessed Data |
|---------|----------------|------------------|
| **File** | processed.csv | unprocessed.csv |
| **Rows** | 6 | ~3,700 |
| **Columns** | 3 | 7 |
| **Markers** | Pre-selected | All markers |
| **ranking_method** | ❌ Not used | ✅ Required (`avg_log2FC`) |
| **ascending** | ❌ Not used | ✅ Required (`false`) |
| **Speed** | Fast (3-5 min) | Slower (5-7 min) |
| **Use Case** | Quick testing | Production-like |

## Common Parameters

Both configs share these settings:

```json
{
  "n_genes": 50,           // Top N genes per cluster
  "tissue": "large intestine",
  "species": "human",
  "max_workers": 4,        // Parallel processing
  "max_retries": 1,        // API retry count
  "validator_involvement": "v1"
}
```

## Troubleshooting

### Error: KeyError 'avg_log2FC'

**Cause**: Using `ranking_method` with processed data

**Fix**: Remove `ranking_method` and `ascending` from config
```json
// ❌ Wrong for processed data
{
  "ranking_method": "avg_log2FC",
  "ascending": false
}

// ✅ Correct for processed data
{
  // No ranking parameters
}
```

### Error: Column not found

**Cause**: Wrong data file for the config

**Fix**: Match config to data format
- `processed.csv` → Use `config.json`
- `unprocessed.csv` → Use `config_unprocessed.json`

### Test runs but no results

**Cause**: API key issues or model errors

**Check**:
```bash
echo $OPENROUTER_API_KEY  # Should not be empty
```

## Expected Results

### Processed Data Output

**File**: `results/YYYYMMDD_HHMMSS_processed_results.csv`

**Structure**:
```csv
True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types
monocyte,Classical Monocyte,CD14+ Classical Monocyte
plasma cell,Plasma Cell,IgA+ Plasma Cell
...
```

### Unprocessed Data Output

**File**: `results/YYYYMMDD_HHMMSS_unprocessed_results.csv`

**Structure**: Same as processed (CASSIA normalizes output format)

```csv
True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types
monocyte,Classical Monocyte,CD14+ Classical Monocyte
plasma cell,Plasma Cell,IgA+ Plasma Cell
...
```

## Performance Benchmarks

| Test | Data | Rows | Runtime | API Calls |
|------|------|------|---------|-----------|
| test_batch.py | Processed | 6 | 3-5 min | 6 |
| test_both.py (part 1) | Processed | 6 | 3-5 min | 6 |
| test_both.py (part 2) | Unprocessed | 6 | 5-7 min | 6 |
| **Total (both)** | Both | 6 | **8-12 min** | **12** |

*Note: Unprocessed is slower due to marker processing overhead*

## When to Use Which Test

### Use test_batch.py (Processed) When:
- ✅ Quick validation of CASSIA installation
- ✅ Testing after code changes
- ✅ CI/CD pipeline (fast tests)
- ✅ Learning CASSIA basics

### Use test_both.py (Both Formats) When:
- ✅ Comprehensive validation before release
- ✅ Testing different data formats
- ✅ Validating ranking logic
- ✅ Production readiness testing

### Use Custom Config When:
- ✅ Testing your own data
- ✅ Different ranking methods
- ✅ Custom tissue/species
- ✅ Different model/provider

## Advanced Usage

### Create Custom Config for Your Data

```json
{
  "test_name": "my_custom_test",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "data_file": "my_data.csv",
  "method_params": {
    "output_name": "my_results",
    "n_genes": 100,                    // Adjust as needed
    "tissue": "your_tissue",
    "species": "your_species",
    "ranking_method": "avg_log2FC",    // If unprocessed
    "ascending": false                 // If unprocessed
  }
}
```

### Run with Custom Config

```bash
# Edit test_batch.py line 51
config = load_config(test_dir / "my_config.json")

python test_batch.py
```

## Summary

✅ **Two data formats supported**
✅ **Three test scripts available**
✅ **Clear configuration examples**
✅ **Comprehensive documentation**

**Recommended first test**: `python test_batch.py` (fastest)
**Recommended full test**: `python test_both.py` (comprehensive)

---

**Last Updated**: 2025-10-08
**Test Version**: 1.1
**Status**: ✅ Both formats validated
