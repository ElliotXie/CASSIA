---
title: Quality Scoring
---


Quality scoring helps evaluate the reliability of cell type annotations. CASSIA provides automated scoring functionality through the `runCASSIA_score_batch` function, which analyzes the reasoning and evidence behind each annotation.

## Running Quality Scoring

### Basic Usage
```R
runCASSIA_score_batch(
    input_file = "my_annotation_full.csv",
    output_file = "my_annotation_scored.csv",
    max_workers = 4,
    model = "google/gemini-2.5-flash-preview",
    provider = "openrouter"
)
```

### Parameter Details

- **Input/Output Files**:
   - `input_file`: Path to the full annotation results (from `runCASSIA_batch`)
   - `output_file`: Where to save the scored results
   
- **Processing Parameters**:
   - `max_workers`: Number of parallel scoring threads
   - Recommended: Use fewer workers than annotation step to avoid API limits if provider set to anthropic

- **Model Configuration**:
   - Recommended model: `google/gemini-2.5-flash-preview`
   - Recommended provider: `openrouter`

### API Provider Considerations

#### OpenRouter
- **Advantages**:
  - Higher rate limits
  - Easy to switch models
- **Setup**:
  ```R
  provider <- "openrouter"
  model <- "google/gemini-2.5-flash-preview"
  ```

#### Anthropic Direct
- **Considerations**:
  - New users have usage limits
  - May need to reduce `max_workers`
  - Better for smaller datasets
- **Setup**:
  ```R
  provider <- "anthropic"
  model <- "claude-3-5-sonnet-20241022"
  ```

### Output Format
The scored output file contains:
- Original annotation data
- Quality scores (0-100)
- Confidence metrics
- Detailed reasoning for scores

### Interpreting Scores

- **90-100**: High confidence, strong evidence
- **76-89**: Good confidence, adequate evidence
- **<75**: Low confidence, need to run through Annotation Boost Agent and Compare Agent

# Report Generation

Generate detailed reports from your analysis. This step typically follows after quality scoring.

The score report includes all outputs from CASSIA, including structured outputs, conversation histories, and quality scores.

### Batch Reports from Scored Results

```R
runCASSIA_generate_score_report(
  csv_path = "my_annotation_scored.csv",
  output_name = "CASSIA_reports_summary"
)
```

_Generates individual reports and an index page from `scored_results.csv`._
