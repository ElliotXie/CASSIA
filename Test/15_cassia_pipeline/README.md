# Test 15: CASSIA Pipeline

## Purpose
Tests the `runCASSIA_pipeline()` function, which is the complete end-to-end cell type annotation orchestrator.

## What it Tests
- Full pipeline execution from markers to final results
- Output directory structure creation
- Initial annotation via `runCASSIA_batch`
- Annotation merging (optional)
- Quality scoring
- HTML report generation
- Annotation boosting for low-scoring clusters
- Result organization into structured directories

## Pipeline Steps
1. Initial LLM-based cell type annotation
2. Optional annotation merging
3. Quality scoring (0-100)
4. HTML report generation
5. Annotation boosting for low-scoring clusters (below threshold)
6. Final result consolidation

## Expected Output Structure
```
CASSIA_[tissue]_[species]_[timestamp]/
├── 01_annotation_results/
│   ├── FINAL_RESULTS.csv
│   └── intermediate_files/
├── 02_reports/
│   └── *.html
└── 03_boost_analysis/
```

## Running the Test

### Python
```bash
python test_cassia_pipeline.py
```

### R
```bash
Rscript test_cassia_pipeline.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- Pipeline output directory with all generated files
