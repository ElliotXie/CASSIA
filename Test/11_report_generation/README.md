# Test 11: Report Generation

## Purpose
Tests the `generate_batch_html_report()` function for generating interactive HTML reports from CASSIA batch analysis results, without needing to re-run the full batch analysis.

## What it Tests
- HTML report generation from existing CSV results
- Report generation with custom title
- Report generation from in-memory data
- HTML file creation and validation
- Interactive features (search, filters, modals)

## Use Case
After running `runCASSIA_batch()`, you may want to regenerate the HTML report with a different title, or generate a report from previously saved CSV results without re-running the entire analysis.

## Functions Tested

### generate_batch_html_report()
Main function for report generation from CSV:
- Takes path to `batch_results_full.csv`
- Optionally specify output path and title
- Returns path to generated HTML file

### generate_batch_html_report_from_data()
Alternative function for report generation from data:
- Takes list of dictionaries (rows)
- Useful when data is already in memory

## Test Parameters
- **Input**: Existing batch results CSV (from Test 02) or sample data
- **Output**: Interactive HTML report
- **Custom title**: "CASSIA Test Report"

## Expected Output
```
Report Generation Results:
  Input CSV: batch_results_full.csv (6 rows)
  Output HTML: test_report.html
  File size: ~150KB

  Report features:
    - Cluster cards: 6
    - Search functionality: Yes
    - Filter dropdowns: Yes
    - Modal popups: Yes
```

## Running the Test

### Python
```bash
python test_report_generation.py
```

### R
```bash
Rscript test_report_generation.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Report generation results
- `test_report.html`: Generated HTML report

## Notes
- Test uses existing batch results from Test 02 if available
- Falls back to creating sample data if no batch results exist
- Does NOT re-run batch analysis (that's the point!)
- Generated HTML is self-contained with embedded CSS/JS
