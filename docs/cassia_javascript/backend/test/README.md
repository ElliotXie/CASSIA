# CASSIA JavaScript Tests

This directory contains comprehensive tests for the CASSIA JavaScript implementation.

## Test Files

- **`test_import.js`** - Basic import verification (no API calls)
- **`test_simple.js`** - Simple functionality test with LLM calls  
- **`test_quick.js`** - Quick test for basic functionality ✅ **Default test**
- **`test_basic.js`** - Comprehensive functionality test
- **`test_comprehensive.js`** - Full parameter compatibility test
- **`test_scoring.js`** - Scoring system functionality test ⭐ **New**
- **`test_reports.js`** - HTML report generation test ⭐ **New**
- **`test_annotation_boost.js`** - Annotation boost iterative analysis test 🚀 **New**
- **`test_sample_data.csv`** - Sample data for testing (3 cell types)

## Test Organization

### 📁 **test_results/** 
All test output files are automatically saved to this directory:
- CSV result files (`*_full.csv`, `*_summary.csv`)
- Analysis reports
- Generated outputs from batch processing

### 🧪 **Running Tests**

```bash
# Quick test (recommended)
npm test

# Specific tests
npm run test-import      # Import verification only
npm run test-simple      # Simple LLM functionality  
npm run test-basic       # Comprehensive functionality
npm run test-comprehensive  # Full parameter testing
npm run test-scoring     # Scoring system testing
npm run test-reports     # HTML report generation
npm run test-annotation-boost  # Annotation boost testing

# Clean up results
npm run clean-results
```

### ⚙️ **Test Configuration**

All tests use the API key configured in the test files. To run tests:

1. **Set API Key**: Edit the test files or set environment variable
2. **Choose Test**: Run appropriate test based on your needs
3. **Check Results**: Review outputs in `test_results/` directory

### 📊 **Test Coverage**

- ✅ **LLM Provider Integration** (OpenAI, Anthropic, OpenRouter, Custom)
- ✅ **Single Cluster Analysis** (`runCASSIA`)
- ✅ **Batch Processing** (`runCASSIABatch`) 
- ✅ **Scoring System** (`scoreSingleAnalysis`, `runCASSIAScoreBatch`) ⭐ **New**
- ✅ **Report Generation** (`runCASSIAGenerateScoreReport`) ⭐ **New**
- ✅ **Annotation Boost** (`runCASSIAAnnotationboost`, iterative analysis) 🚀 **New**
- ✅ **All Validator Versions** (v0, v1, v2, tissue-blind)
- ✅ **All Ranking Methods** (avg_log2FC, p_val_adj, pct_diff, Score)
- ✅ **Complete Parameter Set** (100% Python compatibility)
- ✅ **Error Handling** (retry logic, authentication errors)
- ✅ **Output Generation** (CSV files, conversation history)

### 📋 **Sample Test Results**

The tests will generate files like:
```
test_results/
├── quick_test_full.csv          # Complete results with conversation history
├── quick_test_summary.csv       # Summary results without conversation 
├── ranking_avg_log2FC_full.csv  # Ranking method tests
├── ranking_p_val_adj_full.csv   # Different ranking tests
└── full_params_test_full.csv    # Comprehensive parameter tests
```

### 🎯 **Expected Results**

- **Single Analysis**: Should identify T cells correctly from CD3D, CD3E, CD8A markers
- **Batch Analysis**: Should process 3 cell types (monocytes, T cells, B cells) 
- **All Tests Passing**: Indicates 100% compatibility with Python implementation