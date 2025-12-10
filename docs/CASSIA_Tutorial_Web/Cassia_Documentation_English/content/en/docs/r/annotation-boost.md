---
title: Annotation Boost (Optional)
---


Annotation Boost is an advanced validation tool that enhances annotation confidence through multiple iterations of analysis. It's particularly useful for:
- Validating low-confidence annotations
- Getting detailed insights into specific cell clusters
- Resolving ambiguous cell type assignments
- Generating comprehensive validation reports

### Required Components

1. **Input Data**:
   - Full results CSV from CASSIA batch analysis
   - Original marker gene file (***Note: The marker file should not be filtered!***)
   - Cluster context information
   - Specific cluster identifier

2. **Model Configuration**:
   - Use a model at least more advanced than `gemini-2.5-flash`
   - Recommended: `anthropic/claude-sonnet-4.5` via OpenRouter or Anthropic

### Running Annotation Boost

```R
# Setup parameters
validation_config <- list(
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)

# Define cluster information
cluster_info <- "Human PBMC"

# Specify the cluster you want to validate
target_cluster = "CD4+ T cell"

# Run validation
runCASSIA_annotationboost(
    # Required parameters
    full_result_path = "batch_results_summary.csv",
    marker = marker_data,
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    output_name = "Cluster1_report",

    # Optional parameters
    num_iterations = 5,             # Number of validation rounds
    model = validation_config$model,
    provider = validation_config$provider,

    # Conversation history (from batch annotation)
    conversations_json_path = "batch_results_conversations.json",
    conversation_history_mode = "full",  # "full", "final" (summarize), or "none"

    # Advanced options
    search_strategy = "breadth",    # "breadth" or "depth"
    report_style = "per_iteration", # "per_iteration" or "total_summary"
    validator_involvement = "v1"    # "v0" (high) or "v1" (moderate)
)
```

### Parameter Details

- **`full_result_path`**: Path to the CASSIA results CSV file (`_summary.csv`).
- **`marker`**: Marker gene data (data frame or path). **Important**: Use the same marker data as the initial analysis (do not filter).
- **`cluster_name`**: Exact name of the target cluster to validate.
- **`major_cluster_info`**: Context about the dataset (e.g., "Human PBMC", "Mouse Brain").
- **`output_name`**: Base name for the output validation report.
- **`num_iterations`**: Number of validation rounds (default: 5).
- **`model`**: LLM model to use for validation.
- **`provider`**: API provider for the model.
- **`conversations_json_path`**: Path to the conversations JSON file from batch annotation (`_conversations.json`). This provides the prior annotation context for deeper analysis.
- **`conversation_history_mode`**: How to use prior conversation history.
    - `"full"` (default): Use the complete prior conversation history as context.
    - `"final"`: Summarize the history using an LLM before using it as context.
    - `"none"`: Don't use any prior conversation history.
- **`search_strategy`**: Strategy for exploring hypotheses.
    - `"breadth"` (default): Tests multiple hypotheses in parallel.
    - `"depth"`: Focuses deeply on verifying a single hypothesis.
- **`report_style`**: Format of the final report.
    - `"per_iteration"`: Detailed breakdown of each round.
    - `"total_summary"`: Concise overview.
- **`validator_involvement`**: Level of validation strictness.
    - `"v1"` (default): Moderate involvement.
    - `"v0"`: High involvement (stricter checks).

### Output Files

The analysis generates the following output files:
- `{output_name}_summary.html`: HTML report with detailed analysis results and visualizations.
- `{output_name}_raw_conversation.txt`: Raw conversation text containing the full analysis dialogue.

### Troubleshooting

1. **Low Confidence Results**:
   - Increase `num_iterations` to gather more evidence.
   - Review marker gene quality.
   - Check for potential issues like doublets, batch effects, or background contamination.

2. **Inconsistent Results**:
   - Verify input data quality and consistency.
   - Consider biological variability or mixed populations (subclustering might be needed).
