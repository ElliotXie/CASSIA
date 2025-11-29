# Test 02: Batch Annotation

## Purpose
Tests the `runCASSIA_batch()` function by annotating all 6 cell clusters in parallel.

## What it Tests
- Batch processing functionality
- Parallel execution with multiple workers
- CSV output generation
- All clusters processed successfully

## Test Clusters (6 total)
1. monocyte
2. plasma cell
3. cd8-positive, alpha-beta t cell
4. transit amplifying cell of large intestine
5. intestinal enteroendocrine cell
6. intestinal crypt stem cell

## Expected Output
- `batch_results_full.csv`: Full annotation results
- `batch_results_summary.csv`: Summary of results
- All 6 clusters should be annotated

## Running the Test

### Python
```bash
python test_batch_annotation.py
```

### R
```bash
Rscript test_batch_annotation.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `batch_results_full.csv`: Full results
- `batch_results_summary.csv`: Summary
