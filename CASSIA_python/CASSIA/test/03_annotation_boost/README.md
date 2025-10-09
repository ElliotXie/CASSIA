# Test: runCASSIA_annotationboost

## Overview

**Method**: `CASSIA.runCASSIA_annotationboost()`
**Purpose**: Test iterative marker analysis for deep-dive annotation
**Category**: Advanced Annotation
**Priority**: High
**Expected Runtime**: 5-10 minutes

## Description

This test validates the **annotation boost functionality** which performs iterative marker analysis for ambiguous or uncertain clusters. The LLM can request additional markers iteratively to refine its annotation, making this particularly useful for:

- Clusters with low initial annotation scores
- Ambiguous cell types requiring deeper investigation
- Cases where standard batch annotation needs refinement

### What This Test Validates

- ✅ Iterative marker request mechanism
- ✅ Conversation history tracking
- ✅ Multiple iteration support (up to configured max)
- ✅ Gene checking tool usage by LLM
- ✅ Conversation JSON structure
- ✅ HTML summary report generation
- ✅ Raw conversation text export
- ✅ Integration with batch annotation results

## Prerequisites

### Required

- Python 3.8+
- CASSIA package installed
- OpenRouter API key
- Prerequisite batch annotation results
- ~5-10 minutes execution time

### Environment Setup

```bash
# Set API key
export OPENROUTER_API_KEY='your-api-key-here'

# Verify CASSIA is installed
python -c "import CASSIA; print(CASSIA.__version__)"
```

### Data Dependencies

**IMPORTANT**: This test requires batch annotation results as input. The test will automatically create prerequisite batch results if configured to do so.

## Quick Start

```bash
# Navigate to test directory
cd test/03_annotation_boost/

# Run test (will create prerequisites automatically)
python test_annotation_boost.py
```

## Configuration

The test is configured via `config.json`:

```json
{
  "test_name": "test_annotation_boost",
  "description": "Test iterative marker analysis for deep-dive annotation",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0,
  "data_file": "processed.csv",
  "test_cluster": "monocyte",
  "method_params": {
    "cluster_name": "monocyte",
    "major_cluster_info": "human large intestine",
    "output_name": "boost_test_monocyte",
    "num_iterations": 5,
    "conversation_history_mode": "final",
    "report_style": "per_iteration"
  },
  "prerequisites": {
    "requires_batch_results": true,
    "batch_results_file": "boost_test_batch_full.csv"
  },
  "validation": {
    "check_conversation_history": true,
    "min_iterations": 1,
    "max_iterations": 10,
    "check_html_report": true
  }
}
```

### Key Parameters

- **cluster_name**: Which cluster to perform deep-dive analysis on
- **major_cluster_info**: Tissue/species context (e.g., "human large intestine")
- **num_iterations**: Maximum iterations allowed for marker requests
- **conversation_history_mode**: How to track conversation
  - `"final"`: Only final state
  - `"full"`: Complete history
- **report_style**: Report format
  - `"per_iteration"`: Show each iteration separately
  - `"total_summary"`: Consolidated summary

### Customization

Edit `config.json` to:

- **Test different clusters** by changing `cluster_name` and `test_cluster`
- **Adjust iteration limit** with `num_iterations` (1-10 recommended)
- **Change context** via `major_cluster_info` for different tissues/species
- **Use existing batch results** by setting `requires_batch_results: false`

## Expected Outputs

### Console Output

```
================================================================================
CASSIA Test: runCASSIA_annotationboost (Iterative Deep-Dive)
================================================================================
✓ API key validated for provider: openrouter
✓ Loaded 6 rows from processed.csv
  - Available clusters: ['monocyte', 'plasma cell', ...]
  - Test cluster: monocyte

================================================================================
Test configuration:
  Test: test_annotation_boost
  Model: google/gemini-2.5-flash-preview
  Provider: openrouter
  Cluster: monocyte
  Max iterations: 5
================================================================================

================================================================================
CREATING PREREQUISITE BATCH RESULTS
================================================================================
Running batch annotation to create prerequisite data...
Output name: boost_test_batch

=== Starting cell type analysis ===
Analyzing monocyte...
Analyzing plasma cell...
...
✓ Cell type analysis completed

Starting: Prerequisite batch annotation
Completed: Prerequisite batch annotation (125.43s = 2.09min)
✓ Created prerequisite batch results: boost_test_batch_full.csv

================================================================================
RUNNING ANNOTATION BOOST
================================================================================

Boost configuration:
  - Cluster: monocyte
  - Context: human large intestine
  - Max iterations: 5
  - Conversation mode: final
  - Report style: per_iteration

Starting: Annotation boost analysis

[LLM requests markers iteratively...]

Completed: Annotation boost analysis (245.67s = 4.09min)

✓ Annotation boost completed
  - Result length: 1523 characters
  - Messages count: 8

================================================================================
VALIDATING BOOST OUTPUTS
================================================================================
  ✓ Found conversation_json: boost_test_monocyte_conversation.json
  ✓ Found html_summary: boost_test_monocyte_summary.html
  ✓ Found raw_text: boost_test_monocyte_raw_conversation.txt

Validating conversation history:

Conversation history structure:
  - Type: <class 'list'>
  - Messages: 8
  - Iterations with gene checks: 3

Sample conversation messages:
  [0] user: You are analyzing markers for cluster 'monocyte'...
  [1] assistant: I'll analyze these markers. Let me check additional genes...
  [2] user: Here are the markers for genes: ['CD14', 'FCGR3A']...

================================================================================
ARCHIVING RESULTS
================================================================================
  ✓ Archived conversation_json: 20251007_143052_conversation_json.json
  ✓ Archived html_summary: 20251007_143052_html_summary.html
  ✓ Archived raw_text: 20251007_143052_raw_text.txt
  ✓ Cleaned up prerequisite files
  ✓ Saved test summary: 20251007_143052_boost_summary.json

================================================================================
BOOST ANALYSIS SUMMARY
================================================================================
Cluster analyzed: monocyte
Boost outputs created: 3
Files archived: 4

Test completed successfully (371.10s = 6.19min)

================================================================================
✅ ANNOTATION BOOST TEST COMPLETED SUCCESSFULLY
📁 Results directory: results/
📄 Log file: results/20251007_143052_test_log.txt
================================================================================
```

### Generated Files

```
test/03_annotation_boost/
│
├── config.json                                    # Test configuration
├── test_annotation_boost.py                       # Test script
├── README.md                                      # This file
│
└── results/                                       # Timestamped results
    ├── 20251007_143052_conversation_json.json     # Conversation history
    ├── 20251007_143052_html_summary.html          # HTML report
    ├── 20251007_143052_raw_text.txt              # Raw conversation
    ├── 20251007_143052_boost_summary.json        # Test metadata
    └── 20251007_143052_test_log.txt              # Detailed log
```

### Conversation JSON Structure

The conversation JSON contains the full message history:

```json
[
  {
    "role": "user",
    "content": "You are analyzing markers for cluster 'monocyte' in human large intestine..."
  },
  {
    "role": "assistant",
    "content": "Based on the markers, I'll perform initial analysis...\n\n<check_genes>CD14,FCGR3A,S100A8</check_genes>"
  },
  {
    "role": "user",
    "content": "Here are the expression values for requested genes:\nCD14: 8.5\nFCGR3A: 6.2\nS100A8: 7.1"
  },
  {
    "role": "assistant",
    "content": "With these additional markers, I can now refine the annotation..."
  }
]
```

### HTML Summary Report

The HTML summary includes:
- Cluster information
- Iteration-by-iteration analysis
- Final annotation conclusion
- Markers checked at each step
- Confidence assessment

## Annotation Boost Workflow

### How It Works

1. **Initial Context**: LLM receives cluster name, tissue/species, and initial markers
2. **Analysis Phase**: LLM analyzes available markers
3. **Gene Request**: LLM can request additional genes using `<check_genes>` tag
4. **Gene Lookup**: System retrieves expression values for requested genes
5. **Iteration**: Steps 3-4 repeat up to `num_iterations` times
6. **Final Annotation**: LLM provides refined annotation based on all information

### Example Iteration Flow

```
Iteration 0:
  User: "Here are the top markers for monocyte: [LYZ, CD68, CST3...]"
  LLM: "These suggest myeloid lineage. Let me check: <check_genes>CD14,FCGR3A</check_genes>"

Iteration 1:
  User: "CD14=8.5, FCGR3A=6.2"
  LLM: "High CD14 and moderate FCGR3A suggest classical monocytes. Let me verify: <check_genes>CD16,CCR2</check_genes>"

Iteration 2:
  User: "CD16=2.1, CCR2=7.8"
  LLM: "Low CD16 and high CCR2 confirm classical monocyte phenotype. Final annotation: Classical Monocyte (CD14+ CCR2+)"
```

### Gene Check Mechanism

The LLM uses XML-style tags to request genes:

```xml
<check_genes>GENE1,GENE2,GENE3</check_genes>
```

The system parses these tags and provides expression values from the marker data.

## Validation Checks

### 1. Conversation History

- ✅ JSON structure is list of message dicts
- ✅ Each message has 'role' and 'content'
- ✅ At least one iteration with gene checks
- ✅ Iteration count within configured bounds

### 2. Output Files

- ✅ Conversation JSON created
- ✅ HTML summary generated
- ✅ Raw text conversation saved
- ✅ Files properly formatted and non-empty

### 3. Iteration Behavior

- ✅ Gene check tags properly parsed
- ✅ Requested genes found in data
- ✅ Iterations respect configured maximum
- ✅ Conversation progresses logically

## Sample Results

### Boost Summary JSON

```json
{
  "test_cluster": "monocyte",
  "boost_outputs": ["conversation_json", "html_summary", "raw_text"],
  "conversation_messages": 8,
  "result_length": 1523
}
```

### Conversation Statistics

| Metric | Value |
|--------|-------|
| Total Messages | 8 |
| User Messages | 4 |
| Assistant Messages | 4 |
| Iterations with Gene Checks | 3 |
| Unique Genes Checked | 12 |
| Final Annotation | Classical Monocyte |

## Troubleshooting

### 1. No Gene Checks Performed

**Issue**: LLM doesn't use `<check_genes>` tags

**Possible Causes**:
- Initial markers too comprehensive
- LLM confidence already high
- Model not following instructions

**Solutions**:
```json
{
  "method_params": {
    "num_iterations": 10,  // Allow more chances
    "conversation_history_mode": "full"  // More context
  }
}
```

Or try a different model that better follows structured output instructions.

### 2. Prerequisite Batch Results Missing

**Error**: `Batch results file not found: boost_test_batch_full.csv`

**Solution**:
```json
{
  "prerequisites": {
    "requires_batch_results": true  // Auto-create
  }
}
```

Or manually create batch results first:
```bash
cd test/01_runCASSIA_batch/
python test_batch.py
# Copy results to test/03_annotation_boost/
```

### 3. Too Many/Too Few Iterations

**Issue**: Hits max iterations or stops after 1

**Solutions**:

For more iterations:
```json
{
  "method_params": {
    "num_iterations": 10,
    "conversation_history_mode": "full"
  },
  "validation": {
    "min_iterations": 3  // Expect at least 3
  }
}
```

For fewer iterations (faster testing):
```json
{
  "method_params": {
    "num_iterations": 2
  }
}
```

### 4. Conversation JSON Parse Error

**Error**: `JSON decode error in conversation history`

**Cause**: Malformed conversation JSON

**Solution**:
- Check raw conversation text file for clues
- Verify model output is valid JSON
- Try different model if persistent

### 5. Cluster Not Found in Data

**Error**: `Test cluster 'X' not found in data`

**Solution**: Check available clusters:
```python
import pandas as pd
data = pd.read_csv('../../data/processed.csv')
print(data.iloc[:, 1].unique())  # Shows cluster names
```

Update config with valid cluster name.

### 6. HTML Report Not Generated

**Issue**: Missing `*_summary.html`

**Solutions**:
- Check CASSIA version supports HTML generation
- Verify report_style parameter is valid
- Check for errors in boost execution
- Ensure write permissions in test directory

## Advanced Usage

### Test Multiple Clusters

Run boost analysis on multiple clusters sequentially:

```python
# In test_annotation_boost.py, add loop:
test_clusters = ['monocyte', 'plasma cell', 'cd8_positive_alpha_beta_t_cell']

for cluster in test_clusters:
    config['test_cluster'] = cluster
    config['method_params']['cluster_name'] = cluster
    config['method_params']['output_name'] = f"boost_test_{cluster}"

    # Run boost for this cluster
    result, messages = CASSIA.runCASSIA_annotationboost(...)
```

### Compare Conversation Modes

Test different conversation history modes:

```json
// Test with "final" mode
{
  "method_params": {
    "conversation_history_mode": "final"
  }
}

// vs "full" mode
{
  "method_params": {
    "conversation_history_mode": "full"
  }
}
```

### Use Pre-existing Batch Results

Skip prerequisite creation if you already have batch results:

```json
{
  "prerequisites": {
    "requires_batch_results": false,
    "batch_results_file": "/path/to/existing_batch_results_full.csv"
  }
}
```

### Custom Iteration Limits

Test extreme cases:

```json
// Minimal iterations (fast)
{
  "method_params": {
    "num_iterations": 1
  }
}

// Maximum iterations (thorough)
{
  "method_params": {
    "num_iterations": 15
  }
}
```

### Different Report Styles

```json
// Per-iteration detail
{
  "method_params": {
    "report_style": "per_iteration"
  }
}

// Consolidated summary
{
  "method_params": {
    "report_style": "total_summary"
  }
}
```

## Performance Benchmarks

### Expected Performance (gemini-2.5-flash-preview)

| Component | Time | Notes |
|-----------|------|-------|
| Prerequisite Batch | 2-3 min | 6 clusters |
| Boost Analysis | 3-5 min | 3-5 iterations |
| Validation | <10 sec | Output checks |
| **Total** | **5-10 min** | **Complete test** |

### Iteration Timing

| Iterations | Approximate Time |
|------------|------------------|
| 1 | 1-2 min |
| 3 | 3-4 min |
| 5 | 4-6 min |
| 10 | 8-12 min |

### Resource Usage

- **Memory**: 200-500 MB
- **Disk Space**: ~10-20 MB per test run
- **Network**: ~50-100 KB per iteration
- **CPU**: Single-threaded LLM calls

### Optimization Tips

1. **Reduce prerequisite time** by using smaller dataset:
   ```python
   loader = SampleDataLoader()
   marker_data = loader.load_subset(n_clusters=3)  # Only 3 clusters
   ```

2. **Limit iterations** for faster testing:
   ```json
   {"num_iterations": 2}
   ```

3. **Reuse batch results** across multiple boost tests:
   ```json
   {"requires_batch_results": false}
   ```

4. **Use faster model** for prerequisites:
   ```python
   CASSIA.runCASSIA_batch(..., model="google/gemini-2.5-flash-preview")
   ```

## Interpreting Results

### Good Annotation Boost

Signs of successful boost:
- 2-5 iterations with meaningful gene checks
- Progressive refinement of annotation
- Final annotation more specific than initial
- High confidence in final conclusion

### Poor Annotation Boost

Warning signs:
- No gene checks performed (0 iterations)
- Repetitive gene requests (same genes multiple times)
- No convergence to final annotation
- Contradictory conclusions across iterations

### Example Good Progression

```
Iteration 0: "Likely myeloid cell based on LYZ, CD68"
Iteration 1: "Monocyte subtype with high CD14"
Iteration 2: "Classical monocyte confirmed by CCR2+ CD16-"
Final: "Classical Monocyte (CD14+ CCR2+ CD16-)"
```

### Example Poor Progression

```
Iteration 0: "This is a monocyte"
Iteration 1: "Let me check CD14... yes, monocyte"
Iteration 2: "Let me check CD14 again... still monocyte"
Final: "Monocyte"
```

## Related Tests

- **01_runCASSIA_batch** - Creates prerequisite batch results
- **02_runCASSIA_pipeline** - Uses boost as integrated stage
- **04_merge_annotations** - Downstream analysis of boosted annotations
- **05_uncertainty_quantification** - Multiple runs for confidence

## Use Cases

### When to Use Annotation Boost

1. **Low-scoring clusters** from batch annotation (Score < 75)
2. **Ambiguous cell types** requiring expert-level analysis
3. **Novel cell populations** not well-represented in training
4. **Quality control** to verify batch annotations
5. **Research mode** for exploratory analysis

### When NOT to Use Annotation Boost

1. **High-confidence annotations** (Score > 90)
2. **Well-characterized cell types** (e.g., plasma cells)
3. **Time-sensitive** batch processing
4. **Large-scale** annotation (use batch instead)

## References

- [CASSIA Documentation](../../data/CASSIA_Package_Documentation.md)
- [runCASSIA_annotationboost Source](../../tools_function.py:1410)
- [Testing Framework Plan](../../../../dev_docs/MODULAR_TEST_FRAMEWORK_PLAN.md)
- [Test 02 README](../02_runCASSIA_pipeline/README.md) - Pipeline includes boost
- [Conversation History Docs](../../data/CASSIA_Package_Documentation.md#conversation-modes)

---

**Last Updated**: 2025-10-07
**Version**: 1.0
**Status**: ✅ Active
