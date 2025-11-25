# Test 01: Single Cluster Annotation

## Purpose
Tests the core `runCASSIA()` function by annotating a single cell cluster.

## What it Tests
- Basic annotation functionality
- Model/provider configuration
- Result structure validation

## Test Cluster
- **Cluster**: plasma cell
- **Data source**: `data/markers/processed.csv`

## Expected Output
- `main_cell_type`: Should contain "plasma" or similar
- `sub_cell_types`: Optional list of subtypes

## Running the Test

### Python
```bash
python test_single_annotation.py
```

### R
```bash
Rscript test_single_annotation.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Annotation results
