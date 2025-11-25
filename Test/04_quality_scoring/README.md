# Test 04: Quality Scoring

## Purpose
Tests the `runCASSIA_score_batch()` function to score annotation quality.

## What it Tests
- Quality scoring functionality
- Score generation for each cluster
- Score values in expected range (0-100)

## Prerequisites
- Requires batch annotation results (will run batch annotation if none exist)

## Expected Output
- `scored_results.csv`: Original results with added score column
- Scores should be integers between 0-100
- Average, min, and max scores reported

## Running the Test

### Python
```bash
python test_quality_scoring.py
```

### R
```bash
Rscript test_quality_scoring.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `scored_results.csv`: Results with quality scores
- `results.json`: Scoring statistics
