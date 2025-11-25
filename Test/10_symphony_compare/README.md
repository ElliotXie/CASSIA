# Test 10: Symphony Compare

## Purpose
Tests the `symphonyCompare()` function for multi-model cell type comparison with AI consensus building. Symphony Compare orchestrates multiple AI models to analyze marker expression patterns and reach consensus through structured discussion.

## What it Tests
- Multi-model parallel analysis
- Consensus detection across models
- Discussion rounds when models disagree
- Interactive HTML report generation
- Structured CSV output

## How Symphony Compare Works

1. **Initial Analysis**: Multiple AI models analyze markers independently
2. **Consensus Check**: Calculate agreement level across models
3. **Discussion Rounds**: If no consensus, models review each other's analyses
4. **Final Consensus**: Aggregate results and determine winner

## Model Presets

| Preset | Models | Description |
|--------|--------|-------------|
| symphony | Claude, GPT-4, Gemini Pro | High-performance ensemble |
| quartet | Claude, GPT-4, Gemini Pro, Llama | Balanced 4-model ensemble |
| budget | Gemini Flash, DeepSeek, Grok | Cost-effective models |
| custom | User-defined | Specify custom_models list |

## Test Parameters
- **Cell types compared**: Plasma cell, B cell, T cell
- **Marker set**: Top 15 markers from plasma cell cluster
- **Model preset**: budget (for cost-effective testing)
- **Discussion**: Disabled (for faster testing)
- **Consensus threshold**: 60%

## Expected Output
```
Symphony Compare Results:
  Consensus: Plasma cell
  Confidence: 85.0%
  CSV file: symphony_test.csv
  HTML report: symphony_test_report.html

  Summary:
    Models used: 3
    Total rounds: 1
    Consensus reached: True

  Cell Type Scores:
    Plasma cell: mean=85.0, range=[80.0-90.0]
    B cell: mean=45.0, range=[40.0-50.0]
    T cell: mean=15.0, range=[10.0-20.0]
```

## Running the Test

### Python
```bash
python test_symphony_compare.py
```

### R
```bash
Rscript test_symphony_compare.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Symphony compare results summary
- `symphony_test.csv`: Detailed model responses and scores
- `symphony_test_report.html`: Interactive visual report

## Output Files

### CSV Format
| Column | Description |
|--------|-------------|
| model | AI model used |
| researcher | Persona name (e.g., "Dr. Ada Lovelace") |
| round | Analysis round (initial, discussion_1, etc.) |
| {celltype}_score | Match score (0-100) for each cell type |
| {celltype}_reasoning | Detailed reasoning for score |
| discussion | Model's critique of other analyses (if applicable) |

### HTML Report
Interactive report with:
- Score progression charts
- Model-by-model analysis
- Discussion round summaries
- Final consensus visualization

## Notes
- Test uses "budget" preset for cost-effective testing
- Discussion rounds disabled to reduce API calls
- Full analysis (symphony preset + discussion) provides best results
- Consensus threshold of 0.8 recommended for production
