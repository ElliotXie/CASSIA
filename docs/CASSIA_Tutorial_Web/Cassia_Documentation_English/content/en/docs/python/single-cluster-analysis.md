---
title: Single Cluster Analysis
---

The `runCASSIA` function analyzes a single cluster of marker genes to identify the cell type.
Note that CASSIA is designed to handle multiple clusters at once. This function is specifically designed for users who only have one cluster to analyze.

### Example

For detailed information on model settings and recommendations, please refer to the **[How to Select Models and Providers](setting-up-cassia.md#how-to-select-models-and-providers)** section.

#### Example Code

```python
import CASSIA

# Parameters
model = "anthropic/claude-sonnet-4.5"  # Model to use
temperature = 0
marker_list = ["CD3D", "CD3E", "CD2", "TRAC"]
tissue = "blood"
species = "human"
additional_info = None
provider = "openrouter"  # or "openai", "anthropic"

# Run the analysis
result, conversation_history, _ = CASSIA.runCASSIA(
    model=model,
    temperature=temperature,
    marker_list=marker_list,
    tissue=tissue,
    species=species,
    additional_info=additional_info,
    provider=provider,
    validator_involvement="v1",
    reasoning="medium"  # Optional: "high", "medium", "low" for compatible models
)

# View structured output
print(result['main_cell_type'])
print(result['sub_cell_types'])

# View conversation history
print(conversation_history)
```

### Parameter Details

- **`model`**: The LLM model to use for analysis. See [Setting Up CASSIA](setting-up-cassia.md) for options.
- **`temperature`**: Controls the randomness of the model's output (0 = deterministic, 1 = creative). Default is 0.
- **`marker_list`**: A list of marker gene names for the single cluster.
- **`tissue`**: The tissue of origin for the sample.
- **`species`**: The species of the sample (e.g., "human", "mouse").
- **`additional_info`**: (Optional) Any extra context about the experiment or sample.
- **`provider`**: The API provider to use ("openrouter", "openai", "anthropic").
- **`validator_involvement`**: The level of validation strictness ("v1" for moderate, "v0" for high).
- **`reasoning`**: (Optional) Controls reasoning depth for compatible models ("high", "medium", "low"). Python also accepts dict format: `{"effort": "high"}`. Omit for standard mode. See [Reasoning Effort Parameter](setting-up-cassia.md#reasoning-effort-parameter) for details.

For more detailed parameter explanations, see [Batch Processing Parameter Details](batch-processing.md#parameter-details).

### Return Value

The function returns a tuple: `(result, conversation_history, _)`

- **`result`**: A dictionary containing:
  - `main_cell_type`: The primary cell type prediction
  - `sub_cell_types`: List of possible subtypes
- **`conversation_history`**: The full conversation with the model for transparency

_Note:_ When using OpenRouter, specify the complete model ID.
