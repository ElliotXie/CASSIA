# Test 06: Annotation Boost

## Purpose
Tests the `runCASSIA_annotationboost()` function for iterative deep analysis of cell type annotations.

## What it Tests
- Annotation boost iterative analysis functionality
- Conversation history summarization
- Gene checking and validation
- Report generation (HTML summary + raw conversation)

## How Annotation Boost Works
Annotation boost performs deeper analysis by:
1. Taking existing batch annotation results
2. Loading conversation history from prior analysis
3. Iteratively checking additional marker genes
4. Generating hypotheses and validating them
5. Producing detailed reports with reasoning

## Prerequisites
- Requires batch annotation results (from Test 02)
- If none exist, the test will run batch annotation first

## Test Parameters
- **Cluster tested**: plasma cell (or first available)
- **Iterations**: 3 (reduced for testing)
- **Search strategy**: breadth (tests multiple hypotheses)
- **Report style**: per_iteration

## Expected Output
- `boost_<cluster>_summary.html`: Formatted summary report
- `boost_<cluster>_raw_conversation.txt`: Full conversation history
- Status: success with valid analysis text

## Search Strategies
- **breadth**: Tests multiple hypotheses per iteration (default)
- **depth**: Focuses on one hypothesis at a time, goes deeper

## Running the Test

### Python
```bash
python test_annotation_boost.py
```

### R
```bash
Rscript test_annotation_boost.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Boost results summary
- `boost_*_summary.html`: Generated summary report
- `boost_*_raw_conversation.txt`: Full conversation log

## Notes
- This test takes longer than others due to iterative LLM calls
- Reduced iterations (3) for faster testing
- Full analysis typically uses 5-10 iterations
