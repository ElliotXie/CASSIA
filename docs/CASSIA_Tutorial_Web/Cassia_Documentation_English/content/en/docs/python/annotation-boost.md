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

The monocyte cluster is sometimes annotated as mixed population of immune cell and neuron/glia cells. Here we use annotation boost agent to test these hypotheses in more detail.

```python
# Run validation plus for the high mitochondrial content cluster
CASSIA.runCASSIA_annotationboost(
    full_result_path = output_name + "_summary.csv",
    marker = unprocessed_markers,
    output_name = "monocyte_annotationboost",
    cluster_name = "monocyte",
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter",
    conversations_json_path = output_name + "_conversations.json",  # Use prior conversation history
    conversation_history_mode = "full",  # "full", "final" (summarize), or "none"
    reasoning = "low"  # Optional: "low", "medium", "high"
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
- **`conversation_history_mode`**: How to use prior conversation history:
  - `"full"` (default): Use the complete prior conversation history as context.
  - `"final"`: Summarize the history using an LLM before using it as context.
  - `"none"`: Don't use any prior conversation history.
- **`search_strategy`**: Strategy for exploring hypotheses ("breadth" or "depth").
- **`report_style`**: Format of the final report ("per_iteration" or "total_summary").
- **`validator_involvement`**: Level of validation strictness ("v1" or "v0").
- **`reasoning`**: (Optional) Reasoning effort level ("low", "medium", "high"). Controls how much the model "thinks" before responding. Only supported by OpenAI GPT-5 series models (e.g., `gpt-5.1`). Via OpenRouter, no additional verification needed. Via direct OpenAI API, identity verification (KYC) is required.

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
