# CASSIA Test Framework - Quick Start Guide

**🎯 Get Started in 2 Minutes**

---

## What's Been Built

A **modular, standardized testing framework** for CASSIA where each method has its own test folder with:
- Configuration files (JSON)
- Test scripts (Python)
- Documentation (README)
- Results tracking (timestamped outputs)

### Current Status
✅ **Infrastructure Complete** - All shared utilities ready
✅ **First Test Ready** - `01_runCASSIA_batch` fully functional
📁 **10 More Tests** - Structure ready, implementation pending

---

## Quick Test Run

### 1. Set API Key
```bash
export OPENROUTER_API_KEY='your-key-here'
```

### 2. Navigate to Test
```bash
cd CASSIA_python/CASSIA/test/01_runCASSIA_batch/
```

### 3. Run Test
```bash
python test_batch.py
```

### 4. Check Results
```bash
ls results/
# See: 20251007_143022_batch_full.csv
#      20251007_143022_batch_summary.csv
#      20251007_143022_test_log.txt
```

---

## File Locations

### 📁 Main Documentation
```
dev_docs/
├── MODULAR_TEST_FRAMEWORK_PLAN.md         # Master plan (550 lines)
├── TEST_FRAMEWORK_IMPLEMENTATION_SUMMARY.md  # What's done (400 lines)
└── TEST_FRAMEWORK_QUICK_START.md          # This file
```

### 📁 Test Framework
```
CASSIA_python/CASSIA/test/
├── README.md                    # Master test documentation (800 lines)
│
├── shared/                      # Shared utilities
│   ├── test_config.py          # Configuration loading
│   ├── test_utils.py           # Logging, timing, validation
│   └── sample_data.py          # Data loading
│
└── 01_runCASSIA_batch/         # First test (COMPLETE)
    ├── README.md               # Test documentation
    ├── config.json             # Test configuration
    ├── test_batch.py           # Test script
    └── results/                # Test outputs
```

---

## Using Shared Utilities

### Load Configuration
```python
from test_config import load_config

config = load_config(Path("config.json"))
model = config['model']  # 'google/gemini-2.5-flash-preview'
```

### Load Sample Data
```python
from sample_data import SampleDataLoader

loader = SampleDataLoader()
data = loader.load_processed()  # 6 clean clusters
subset = loader.load_subset(clusters=["monocyte", "T cell"])
markers = loader.get_marker_list("monocyte", max_genes=50)
```

### Setup Logging
```python
from test_utils import setup_logging, log_test_start, log_test_end

logger = setup_logging("my_test")
start = log_test_start(logger, config.to_dict())
# ... do work ...
log_test_end(logger, start, success=True)
```

### Time Operations
```python
from test_utils import Timer

with Timer("Data processing", logger):
    result = process_data()
    # Automatically logs: "Completed: Data processing (2.5s)"
```

### Save Results
```python
from test_utils import save_results

# Saves with timestamp: 20251007_143022_results.csv
save_results(df, "test_name", Path("results"))
```

---

## Standard Test Structure

### config.json
```json
{
  "test_name": "test_my_method",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "processed.csv",
  "method_params": {
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4
  }
}
```

### test_[method].py
```python
#!/usr/bin/env python3
import sys
from pathlib import Path
import CASSIA
from test_config import load_config, validate_api_key
from test_utils import setup_logging, log_test_start, Timer
from sample_data import load_sample_data

def main():
    config = load_config(Path("config.json"))
    logger = setup_logging(config['test_name'])

    api_key = validate_api_key(config['provider'])
    CASSIA.set_api_key(api_key, provider=config['provider'])

    data = load_sample_data("processed")

    start = log_test_start(logger, config.to_dict())
    try:
        with Timer("CASSIA method", logger):
            result = CASSIA.method_name(**config['method_params'])

        save_results(result, config['test_name'], Path("results"))
        log_test_end(logger, start, success=True)
        return 0
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        log_test_end(logger, start, success=False)
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

---

## Next Steps for You

### Option 1: Run Existing Test ✅
```bash
cd test/01_runCASSIA_batch/
python test_batch.py
```

### Option 2: Create New Test 📝
```bash
# Use the pattern from 01_runCASSIA_batch/
cd test/02_runCASSIA_pipeline/

# Copy template
cp ../01_runCASSIA_batch/config.json .
cp ../01_runCASSIA_batch/test_batch.py test_pipeline.py
cp ../01_runCASSIA_batch/README.md .

# Modify for pipeline test
# - Update config.json
# - Change CASSIA.runCASSIA_batch → CASSIA.runCASSIA_pipeline
# - Update documentation
```

### Option 3: Continue Implementation 🚀
Follow the roadmap in `MODULAR_TEST_FRAMEWORK_PLAN.md`:
- Day 3-5: Core tests (pipeline, annotation_boost, merge)
- Day 6-8: Advanced tests (UQ, subclustering)
- Day 9-10: Comparison & utility tests
- Day 11: Master test runner
- Day 12: Migration & cleanup

---

## Key Features

### ✅ Standardization
- Every test follows same structure
- Consistent naming conventions
- JSON-based configuration
- Timestamped outputs

### ✅ Reusability
- Shared utilities eliminate duplication
- Data loaders handle all formats
- Validation framework for outputs
- Logging infrastructure ready

### ✅ Documentation
- Master README (800 lines)
- Test-specific READMEs
- Comprehensive plan document
- Implementation summary
- This quick start guide

### ✅ Real Data Integration
- Uses actual CASSIA sample datasets
- processed.csv - 6 clean clusters
- unprocessed.csv - Large raw dataset
- subcluster_results.csv - Subclustering data

---

## Expected Test Results

### Test 01: runCASSIA_batch
- **Input**: 6 cell type clusters
- **Output**: Annotated cell types
- **Runtime**: 3-5 minutes
- **Files**: 3 (full CSV, summary CSV, log)

**Sample Output**:
```csv
True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types
monocyte,Classical Monocyte,"CD14+ Monocyte, CD16- Monocyte"
plasma cell,Plasma Cell,"IgA+ Plasma Cell, IgG+ Plasma Cell"
```

---

## Troubleshooting

### API Key Not Set
```bash
export OPENROUTER_API_KEY='sk-or-v1-...'
```

### Import Errors
```bash
cd CASSIA_python/
pip install -e .
```

### Data Not Found
```bash
# Check data directory
ls CASSIA_python/CASSIA/data/
# Should see: processed.csv, unprocessed.csv, etc.
```

### Permission Errors
```bash
# Make test scripts executable
chmod +x test/01_runCASSIA_batch/test_batch.py
```

---

## Resources

### 📚 Documentation
- [Master Test README](../CASSIA_python/CASSIA/test/README.md)
- [Framework Plan](MODULAR_TEST_FRAMEWORK_PLAN.md)
- [Implementation Summary](TEST_FRAMEWORK_IMPLEMENTATION_SUMMARY.md)
- [Test 01 README](../CASSIA_python/CASSIA/test/01_runCASSIA_batch/README.md)

### 🔗 Links
- Main docs in `dev_docs/`
- Test docs in `test/*/README.md`
- Shared code in `test/shared/`
- Sample data in `CASSIA/data/`

---

## Summary

**What You Have**:
- ✅ Complete test framework infrastructure
- ✅ Shared utilities for all tests
- ✅ Comprehensive documentation
- ✅ First working test (runCASSIA_batch)
- ✅ Templates for 10 more tests

**What To Do Next**:
1. Run the existing test to validate setup
2. Read the framework plan for full details
3. Implement remaining tests using templates
4. Create master test runner
5. Archive old test files

**Time to Full Implementation**: ~8-10 days following the roadmap

---

**Ready to Test?** 🧪
```bash
cd CASSIA_python/CASSIA/test/01_runCASSIA_batch/
export OPENROUTER_API_KEY='your-key'
python test_batch.py
```

---

**Last Updated**: 2025-10-07
**Status**: Phase 1 Complete ✅
