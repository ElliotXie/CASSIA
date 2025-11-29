# Test 03: Validator Comparison

## Purpose
Compares the behavior of v0 (strict) and v1 (moderate) validators.

## What it Tests
- Both validator modes work correctly
- Both validators can process all clusters
- Comparison of strictness levels

## Validators
- **v0 (strict)**: More rigorous validation
- **v1 (moderate)**: Balanced validation (default)

## Expected Output
- Two sets of results (one per validator)
- Both should successfully annotate all 6 clusters
- May show differences in annotation confidence

## Running the Test

### Python
```bash
python test_validator_comparison.py
```

### R
```bash
Rscript test_validator_comparison.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results_v0_full.csv`: v0 validator results
- `results_v1_full.csv`: v1 validator results
- `results.json`: Comparison summary
