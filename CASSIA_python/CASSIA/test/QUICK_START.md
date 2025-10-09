# CASSIA Test Framework - Quick Start Guide

## 🚀 Get Started in 2 Minutes

### 1. Set Your API Key

```bash
export OPENROUTER_API_KEY='your-key-here'
```

### 2. Run All Tests

```bash
cd test/
python run_all_tests.py
```

### 3. Run Quick Tests (Fast Validation)

```bash
python run_all_tests.py --quick
```

## 📋 Available Tests

**Fully Implemented (Ready to Use):**
- ✅ Test 01: `runCASSIA_batch` - Basic batch annotation (3-5 min)
- ✅ Test 02: `runCASSIA_pipeline` - Full pipeline (10-15 min)
- ✅ Test 03: `annotation_boost` - Deep-dive analysis (5-10 min)
- ✅ Test 04: `merge_annotations` - Multi-level merging (3-5 min)
- ✅ Test 05: `uncertainty_quantification` - Stability testing (15-20 min)
- ✅ Test 06: `subclustering` - Hierarchical annotation (3-5 min)

**Placeholders (To Be Expanded):**
- 🚧 Tests 07-11: Basic structure in place, full implementation pending

## 🎯 Common Usage Patterns

### Run Specific Tests

```bash
# Run tests 01, 04, and 06 only
python run_all_tests.py 01 04 06
```

### List All Available Tests

```bash
python run_all_tests.py --list
```

### Verbose Output

```bash
python run_all_tests.py --verbose
```

## 📁 Test Results

All test results are automatically saved to timestamped files:
- CSV results: `test/XX_test_name/results/YYYYMMDD_HHMMSS_*.csv`
- Logs: `test/XX_test_name/results/YYYYMMDD_HHMMSS_test_log.txt`
- JSON summaries: `test/XX_test_name/results/YYYYMMDD_HHMMSS_*_summary.json`

## 🔧 Run Individual Tests

```bash
# Navigate to specific test
cd test/01_runCASSIA_batch/

# Run test directly
python test_batch.py
```

## 📖 Documentation

- **Master README**: `test/README.md` - Comprehensive guide
- **Per-test READMEs**: Each test folder has detailed documentation
- **Implementation Plan**: `dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md`
- **Final Summary**: `dev_docs/TEST_FRAMEWORK_FINAL_SUMMARY.md`

## ⚡ Quick Test Mode (Fastest)

Run only the fastest tests for quick validation:

```bash
python run_all_tests.py --quick
# Runs: Test 01, 04, 06 (total ~10-15 min)
```

## 🐛 Troubleshooting

**API Key Issues:**
```bash
# Verify key is set
echo $OPENROUTER_API_KEY

# Check it's valid
python -c "import os; print('Set' if os.getenv('OPENROUTER_API_KEY') else 'NOT SET')"
```

**Test Failures:**
1. Check log files in `test/XX_test_name/results/`
2. Review per-test README for troubleshooting section
3. Try running test individually with verbose output

**Missing Data:**
Tests automatically use sample data from `CASSIA/data/` directory. If files are missing, some tests create synthetic data.

## 💡 Tips

- **First Run**: Start with `python run_all_tests.py 01` to validate setup
- **Development**: Use `--quick` mode for frequent validation
- **Full Validation**: Run all tests before major releases
- **Debugging**: Use verbose mode and check individual test logs

## 🎓 Next Steps

1. **Customize Tests**: Edit `config.json` files in each test folder
2. **Add Tests**: Expand placeholders (Tests 07-11) as needed
3. **Integrate CI/CD**: Use `run_all_tests.py` in automated workflows

## 📞 Help

- Master README: `test/README.md`
- Per-test docs: `test/XX_test_name/README.md`
- Framework plan: `dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md`

---

**Framework Version**: 1.0
**Last Updated**: 2025-10-07
**Status**: ✅ Ready to Use
