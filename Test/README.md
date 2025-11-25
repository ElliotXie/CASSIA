# CASSIA Test Suite

A modular test suite for CASSIA (Python and R) with organized test folders and comprehensive result management.

## Configuration

### Default Settings
- **LLM Provider**: OpenRouter
- **Model**: `google/gemini-2.5-flash`
- **Marker Data**: `data/markers/processed.csv` (6 intestinal cell clusters)

### API Keys Setup

1. Copy the example file:
   ```bash
   cp config/api_keys.env.example config/api_keys.env
   ```

2. Edit `config/api_keys.env` and add your API keys:
   ```bash
   OPENROUTER_API_KEY=your-actual-key-here
   ```

## Quick Start

### Run All Tests (Python)
```bash
python run_all_tests.py
```

### Run All Tests (R)
```bash
Rscript run_all_tests.R
```

### Run Individual Test
```bash
# Python
python run_test.py 01                    # By number
python run_test.py batch                 # By name substring
python run_test.py 03_validator          # By partial name

# R
Rscript run_test.R 01
Rscript run_test.R batch
```

## Test Suite Structure

```
test/
├── config/
│   ├── test_config.yaml          # Main configuration
│   └── api_keys.env.example      # API key template
│
├── shared/
│   ├── python/                   # Shared Python utilities
│   │   ├── fixtures.py           # Data loading
│   │   ├── test_utils.py         # Test utilities
│   │   └── result_manager.py     # Results management
│   │
│   └── r/                        # Shared R utilities
│       ├── fixtures.R
│       ├── test_utils.R
│       └── result_manager.R
│
├── data/markers/
│   └── processed.csv             # Test marker data
│
├── 01_single_annotation/         # Test 01
├── 02_batch_annotation/          # Test 02
├── 03_validator_comparison/      # Test 03
├── 04_quality_scoring/           # Test 04
├── 05_model_settings/            # Test 05
├── 06_annotation_boost/          # Test 06
├── 07_merging_annotation/        # Test 07
├── 08_uncertainty_quantification/ # Test 08
├── 09_subclustering/             # Test 09
├── 10_symphony_compare/          # Test 10
│
├── run_test.py                   # Individual test runner (Python)
├── run_test.R                    # Individual test runner (R)
├── run_all_tests.py              # Full test runner (Python)
└── run_all_tests.R               # Full test runner (R)
```

## Available Tests

| # | Test | Description |
|---|------|-------------|
| 01 | Single Annotation | Test `runCASSIA()` on a single cluster |
| 02 | Batch Annotation | Test `runCASSIA_batch()` on all 6 clusters |
| 03 | Validator Comparison | Compare v0 (strict) vs v1 (moderate) validators |
| 04 | Quality Scoring | Test `runCASSIA_score_batch()` scoring |
| 05 | Model Settings | Test model name resolution and shortcuts |
| 06 | Annotation Boost | Test `runCASSIA_annotationboost()` iterative deep analysis |
| 07 | Merging Annotation | Test `merge_annotations()` grouping at different detail levels |
| 08 | Uncertainty Quantification | Test `runCASSIA_n_times_similarity_score()` robustness analysis |
| 09 | Subclustering | Test `runCASSIA_subclusters()` subcluster annotation |
| 10 | Symphony Compare | Test `symphonyCompare()` multi-model consensus building |

## Test Data

The test suite uses 6 intestinal cell clusters from `processed.csv`:
1. monocyte
2. plasma cell
3. cd8-positive, alpha-beta t cell
4. transit amplifying cell of large intestine
5. intestinal enteroendocrine cell
6. intestinal crypt stem cell

## Results Management

Each test run creates a timestamped results folder:

```
01_single_annotation/results/
└── 20241124_143022/
    ├── test_metadata.json    # Config, duration, status
    ├── results.json          # Test-specific results
    └── output_*.csv          # Generated files
```

## Command Line Options

### run_all_tests.py
```bash
# Run all tests
python run_all_tests.py

# Skip specific tests
python run_all_tests.py --skip 03,04

# Run only specific tests
python run_all_tests.py --only 01,02

# Don't save JSON report
python run_all_tests.py --no-report
```

### run_test.py
```bash
# List available tests
python run_test.py --list

# Run by number
python run_test.py 01

# Run by name
python run_test.py single_annotation
```

## Configuration File

Edit `config/test_config.yaml` to customize settings:

```yaml
llm:
  provider: "openrouter"
  model: "google/gemini-2.5-flash"
  temperature: 0.3
  max_workers: 3

data:
  marker_file: "data/markers/processed.csv"
  tissue: "large intestine"
  species: "human"
  n_genes: 30

results:
  keep_last_n: 10

validator:
  default: "v1"
```

## Adding New Tests

1. Create a new numbered folder (e.g., `06_new_test/`)
2. Add `test_new_test.py` and `test_new_test.R`
3. Add `README.md` describing the test
4. Create `results/` subfolder for outputs

## Troubleshooting

### API Key Issues
- Ensure `OPENROUTER_API_KEY` is set in `config/api_keys.env`
- Verify the key is valid at openrouter.ai

### Import Errors
- Ensure CASSIA Python is installed or accessible
- For R, install CASSIA R package first

### Test Timeouts
- Default timeout is 10 minutes per test
- Reduce `max_workers` in config for slow connections
