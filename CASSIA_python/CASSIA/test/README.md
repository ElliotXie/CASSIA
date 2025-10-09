# CASSIA Testing Framework

**Version**: 1.0.0
**Last Updated**: 2025-10-07
**Author**: Elliot Yixuan Xie & Claude Code

---

## Overview

This is the **modular testing framework** for CASSIA (Cell-type Annotation using Single-cell RNA-sequencing and Artificial Intelligence). Each CASSIA method has its own dedicated test folder with standardized structure, configuration, and documentation.

### Framework Features

- ✅ **Modular Design** - Each method in its own folder
- ✅ **Standardized Structure** - Consistent layout across all tests
- ✅ **Real Sample Data** - Uses actual CASSIA datasets
- ✅ **Configurable** - JSON-based configuration
- ✅ **Automated** - Master test runner for all tests
- ✅ **Documented** - Each test has comprehensive README
- ✅ **Version Controlled** - Results timestamped and traceable

---

## Quick Start

### Prerequisites

1. **Python 3.8+** installed
2. **CASSIA package** installed in development mode
3. **API Key** set in environment:
   ```bash
   export OPENROUTER_API_KEY='your-api-key-here'
   ```

### Run a Single Test

```bash
# Navigate to test directory
cd test/01_runCASSIA_batch/

# Run the test
python test_batch.py
```

### Run All Tests

```bash
# From test/ directory
python run_all_tests.py
```

---

## Directory Structure

```
test/
├── README.md                    # This file
├── run_all_tests.py            # Master test runner
│
├── shared/                     # Shared utilities
│   ├── __init__.py
│   ├── test_config.py         # Configuration loader
│   ├── test_utils.py          # Logging, timing, validation
│   └── sample_data.py         # Data loading helpers
│
├── 01_runCASSIA_batch/        # Core batch annotation
├── 02_runCASSIA_pipeline/     # Full pipeline
├── 03_annotation_boost/       # Iterative deep-dive
├── 04_merge_annotations/      # Annotation merging
├── 05_uncertainty_quantification/  # UQ analysis
├── 06_subclustering/          # Subcluster annotation
├── 07_celltype_comparison/    # Cell type comparison
├── 08_symphony_compare/       # Multi-agent discussion
├── 09_llm_utils/              # LLM utilities
├── 10_model_settings/         # Model configuration
└── 11_report_generation/      # Report generation
```

### Standard Test Directory

Each test folder follows this structure:

```
XX_test_name/
├── README.md              # Test documentation
├── test_[name].py         # Main test script
├── config.json            # Test configuration
├── results/               # Test outputs (CSV, JSON, TXT)
│   ├── .gitkeep
│   └── [timestamp]_*.csv
└── reports/               # HTML reports (if applicable)
    ├── .gitkeep
    └── [timestamp]_*.html
```

---

## Available Tests

### Core Annotation Tests (Priority: High)

| # | Test Name | Method | Runtime | Description |
|---|-----------|--------|---------|-------------|
| 01 | runCASSIA_batch | `runCASSIA_batch()` | 3-5 min | Basic batch annotation |
| 02 | runCASSIA_pipeline | `runCASSIA_pipeline()` | 10-15 min | Full end-to-end pipeline |
| 03 | annotation_boost | `runCASSIA_annotationboost()` | 5-10 min | Iterative marker analysis |

### Post-Processing Tests (Priority: High)

| # | Test Name | Method | Runtime | Description |
|---|-----------|--------|---------|-------------|
| 04 | merge_annotations | `merge_annotations()` | 2-3 min | Annotation merging |

### Advanced Analysis Tests (Priority: Medium)

| # | Test Name | Method | Runtime | Description |
|---|-----------|--------|---------|-------------|
| 05 | uncertainty_quantification | `runCASSIA_batch_n_times()` | 15-20 min | Multiple runs for UQ |
| 06 | subclustering | `runCASSIA_subclusters()` | 5-8 min | Hierarchical annotation |

### Comparison Tests (Priority: Medium)

| # | Test Name | Method | Runtime | Description |
|---|-----------|--------|---------|-------------|
| 07 | celltype_comparison | `compareCelltypes()` | 3-5 min | Cell type comparison |
| 08 | symphony_compare | `symphonyCompare()` | 8-12 min | Multi-agent debate |

### Utility Tests (Priority: Low)

| # | Test Name | Method | Runtime | Description |
|---|-----------|--------|---------|-------------|
| 09 | llm_utils | `call_llm()` | 1-2 min | LLM utilities |
| 10 | model_settings | `resolve_model_name()` | <1 min | Model configuration |
| 11 | report_generation | `generate_html_report()` | 1-2 min | Report generation |

---

## Configuration

### Default Test Settings

All tests use these default settings (customizable in `config.json`):

```json
{
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "processed.csv",
  "method_params": {
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4,
    "n_genes": 50
  }
}
```

### Sample Data

Tests use real sample data from `../data/`:

- **processed.csv** - 6 clean clusters (standard testing)
- **unprocessed.csv** - Large dataset (stress testing)
- **subcluster_results.csv** - Subclustering data

### Environment Variables

Required:
- `OPENROUTER_API_KEY` - Your OpenRouter API key

Optional:
- `OPENAI_API_KEY` - For OpenAI-specific tests
- `ANTHROPIC_API_KEY` - For Anthropic-specific tests

---

## Running Tests

### Individual Test

```bash
# Navigate to test directory
cd test/01_runCASSIA_batch/

# Run test
python test_batch.py
```

**Output**:
```
================================================================================
Starting test: test_runCASSIA_batch
Model: google/gemini-2.5-flash-preview
Provider: openrouter
Data: processed.csv
================================================================================
Loading data: processed.csv
Loaded 6 rows
Starting: CASSIA method execution
...
Completed: CASSIA method execution (125.34s = 2.09min)
Validating output format...
✓ Validation passed
Saving results...
✓ Results saved: results/20251007_143022_results.csv
================================================================================
Test Status: ✓ SUCCESS
Elapsed time: 130.45 seconds (2.17 minutes)
================================================================================
```

### All Tests (Sequential)

```bash
cd test/
python run_all_tests.py
```

### All Tests (Parallel - Future)

```bash
python run_all_tests.py --parallel --workers 3
```

### Specific Tests

```bash
# Run only core tests
python run_all_tests.py --filter "01,02,03"

# Run only utility tests
python run_all_tests.py --filter "09,10,11"
```

---

## Test Results

### Result Files

All results are timestamped and saved in `results/` directories:

```
results/
├── 20251007_143022_results.csv      # Main results
├── 20251007_143022_summary.csv      # Summary (if applicable)
└── 20251007_143022_test_log.txt     # Execution log
```

### Naming Convention

Format: `[YYYYMMDD]_[HHMMSS]_[description].[ext]`

Examples:
- `20251007_143022_batch_full.csv`
- `20251007_143022_pipeline_merged.csv`
- `20251007_143022_boost_conversation.json`

### Reports

HTML reports (when generated) are saved in `reports/`:

```
reports/
└── 20251007_143022_report.html
```

---

## Adding a New Test

### Step 1: Create Test Directory

```bash
cd test/
mkdir 12_new_test_name
cd 12_new_test_name
mkdir results reports
touch results/.gitkeep reports/.gitkeep
```

### Step 2: Create config.json

```json
{
  "test_name": "test_new_method",
  "description": "Test new CASSIA method functionality",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "processed.csv",
  "method_params": {
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4
  },
  "validation": {
    "check_output_format": true,
    "expected_columns": ["cluster", "annotation"],
    "min_rows": 1
  }
}
```

### Step 3: Create test_[method].py

Use the template from `shared/` or copy from existing test and modify.

### Step 4: Create README.md

Document:
- Test purpose
- Expected runtime
- Sample outputs
- Troubleshooting

### Step 5: Test and Validate

```bash
python test_[method].py
```

---

## Shared Utilities

### test_config.py

```python
from test_config import load_config, get_data_path

# Load configuration
config = load_config(Path("config.json"))

# Get data path
data_path = get_data_path("processed.csv")
```

### test_utils.py

```python
from test_utils import (
    setup_logging,
    log_test_start,
    log_test_end,
    save_results,
    validate_output,
    Timer
)

# Setup logging
logger = setup_logging("my_test")

# Time operations
with Timer("Data loading", logger):
    data = load_data()

# Save results
save_results(df, "test_name", results_dir)
```

### sample_data.py

```python
from sample_data import SampleDataLoader, load_sample_data

# Quick load
data = load_sample_data("processed")

# Advanced loading
loader = SampleDataLoader()
subset = loader.load_subset(
    dataset="processed",
    clusters=["monocyte", "plasma cell"]
)
markers = loader.get_marker_list("monocyte", max_genes=50)
```

---

## Validation

Each test validates:

1. **Output Format** - Correct columns and data types
2. **Row Count** - Minimum number of results
3. **Null Values** - No unexpected nulls in critical columns
4. **Business Logic** - Method-specific checks

Example validation:

```python
validate_output(
    result_df,
    expected_columns=["cluster", "annotation", "confidence"],
    min_rows=6,
    max_nulls={"annotation": 0.0}  # No nulls allowed
)
```

---

## Troubleshooting

### API Key Not Set

**Error**:
```
EnvironmentError: OPENROUTER_API_KEY not found in environment
```

**Solution**:
```bash
export OPENROUTER_API_KEY='your-api-key-here'
```

### Data File Not Found

**Error**:
```
FileNotFoundError: Data file not found: .../data/processed.csv
```

**Solution**:
- Ensure you're running from the test directory
- Check that `../data/processed.csv` exists
- Verify data directory structure

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

### Permission Errors

**Error**:
```
PermissionError: [Errno 13] Permission denied: 'results/'
```

**Solution**:
- Check directory permissions
- Ensure results/ directory exists
- Run with appropriate user permissions

---

## Best Practices

### 1. Always Use Shared Utilities

❌ **Don't**:
```python
import json
config = json.load(open("config.json"))
```

✅ **Do**:
```python
from test_config import load_config
config = load_config(Path("config.json"))
```

### 2. Timestamp All Results

❌ **Don't**:
```python
df.to_csv("results.csv")
```

✅ **Do**:
```python
from test_utils import save_results
save_results(df, test_name, results_dir, suffix="results")
```

### 3. Log Everything

```python
logger = setup_logging("test_name")
logger.info("Starting test...")

with Timer("Method execution", logger):
    result = CASSIA.method()

logger.info("✓ Test completed")
```

### 4. Validate Outputs

```python
# Always validate before saving
validate_output(result, expected_columns, min_rows)
save_results(result, test_name, results_dir)
```

### 5. Handle Errors Gracefully

```python
try:
    result = CASSIA.method()
    log_test_end(logger, start_time, success=True)
except Exception as e:
    logger.error(f"Test failed: {str(e)}", exc_info=True)
    log_test_end(logger, start_time, success=False)
    return 1
```

---

## CI/CD Integration (Future)

### GitHub Actions Example

```yaml
name: CASSIA Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.9
      - name: Install dependencies
        run: |
          pip install -e .
      - name: Run tests
        env:
          OPENROUTER_API_KEY: ${{ secrets.OPENROUTER_API_KEY }}
        run: |
          cd test/
          python run_all_tests.py
```

---

## Performance Benchmarks

### Expected Runtimes (gemini-2.5-flash-preview)

| Test | Clusters | Expected | Actual (Avg) |
|------|----------|----------|--------------|
| 01_runCASSIA_batch | 6 | 3-5 min | TBD |
| 02_runCASSIA_pipeline | 6 | 10-15 min | TBD |
| 03_annotation_boost | 1 | 5-10 min | TBD |
| 04_merge_annotations | 6 | 2-3 min | TBD |
| 05_uncertainty_quantification | 3 | 15-20 min | TBD |
| 06_subclustering | 3 | 5-8 min | TBD |
| 07_celltype_comparison | 3 | 3-5 min | TBD |
| 08_symphony_compare | 1 | 8-12 min | TBD |

**Total Suite Runtime**: ~60-90 minutes

---

## Migration from Old Tests

Old test files from `test_code/` have been archived to `test_code_legacy/`.

### Mapping Old → New

| Old File | New Location |
|----------|--------------|
| test.py | 01_runCASSIA_batch/, 02_runCASSIA_pipeline/ |
| test_batch_analysis.py | 01_runCASSIA_batch/ |
| run_batch_simple.py | 01_runCASSIA_batch/ |
| quick_test.py | 09_llm_utils/ |
| test_import.py | (Deprecated) |
| test_model_settings.ipynb | 10_model_settings/ |

---

## Contributing

### Adding New Tests

1. Create test directory following naming convention
2. Use shared utilities consistently
3. Follow standard structure (config.json, test script, README)
4. Document expected outputs and runtime
5. Add validation checks
6. Update this master README

### Code Style

- Follow PEP 8
- Type hints where appropriate
- Comprehensive docstrings
- Clear variable names
- Error handling

### Commit Messages

- `test: Add test for [method]`
- `test: Fix [test] validation`
- `test: Update [test] configuration`

---

## Resources

### Documentation

- [CASSIA Package Documentation](../data/CASSIA_Package_Documentation.md)
- [Testing Framework Plan](../../../dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md)
- [Individual Test READMEs](./01_runCASSIA_batch/README.md)

### Support

- **Issues**: https://github.com/elliotxe/CASSIA/issues
- **Discussions**: https://github.com/elliotxe/CASSIA/discussions

---

## Changelog

### Version 1.0.0 (2025-10-07)

- ✅ Initial framework implementation
- ✅ 11 test modules created
- ✅ Shared utilities implemented
- ✅ Sample data integration
- ✅ Master test runner (pending)
- ✅ Comprehensive documentation

---

## License

Same as CASSIA package - see LICENSE file in repository root.

---

**Happy Testing! 🧪🧬**
