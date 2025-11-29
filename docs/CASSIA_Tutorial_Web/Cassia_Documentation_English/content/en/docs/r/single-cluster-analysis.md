---
title: Single Cluster Analysis
---

The `runCASSIA` function analyzes a single cluster of marker genes to identify the cell type.
Note that CASSIA is designed to handle multiple clusters at once, this function is specifically designed for users who only have one cluster to analyze.

### Example

For detailed information on model settings and recommendations, please refer to the **[How to Select Models and Providers](setting-up-cassia.md#how-to-select-models-and-providers)** section.

#### Example Code

```R
# Parameters
model <- "anthropic/claude-4.5-sonnet"  # Model to use
temperature <- 0
marker_list <- c("CD3D", "CD3E", "CD2", "TRAC")
tissue <- "blood"
species <- "human"
additional_info <- NULL
provider <- "openrouter"  # or "openai", "anthropic"

# Run the analysis
result <- runCASSIA(
  model = model,
  temperature = temperature,
  marker_list = marker_list,
  tissue = tissue,
  species = species,
  additional_info = additional_info,
  provider = provider,
  validator_involvement = "v1"
)

# View structured output
print(result$structured_output)

# View conversation history
print(result$conversation_history)
```

### Parameter Details

- **`model`**: The LLM model to use for analysis. See [Setting Up CASSIA](setting-up-cassia.md) for options.
- **`temperature`**: Controls the randomness of the model's output (0 = deterministic, 1 = creative). Default is 0.
- **`marker_list`**: A character vector of marker gene names for the single cluster.
- **`tissue`**: The tissue of origin for the sample.
- **`species`**: The species of the sample (e.g., "human", "mouse").
- **`additional_info`**: (Optional) Any extra context about the experiment or sample.
- **`provider`**: The API provider to use ("openrouter", "openai", "anthropic").
- **`validator_involvement`**: The level of validation strictness ("v1" for moderate, "v0" for high).

_Note:_ When using OpenRouter, specify the complete model ID.
