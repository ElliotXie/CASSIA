---
title: Single Cluster Analysis
---

The `runCASSIA` function analyzes a single cluster of marker genes to identify the cell type.
Note that CASSIA is designed to handle multiple clusters at once, this function is specifically designed for users who only have one cluster to analyze.

### Example

If you're using OpenRouter as your provider, you can specify models like `"openai/gpt-4o-2024-11-20"` or `"anthropic/claude-3.5-sonnet"`. Here are some model recommendations:

- **Gemini 2.5 Flash** (Best performance, cost-effective)
    - Model ID: `"google/gemini-2.5-flash-preview"`
- **Deepseek v3** (Open source, excellent performance, cost-effective)
    - Model ID: `"deepseek/deepseek-chat-v3-0324"`
- **Deepseek v3 Free** (Free tier, slower)
    - Model ID: `"deepseek/deepseek-chat-v3-0324:free"`

#### Example Code

```R
# Parameters
model <- "google/gemini-2.5-flash-preview"  # Model ID when using OpenRouter
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
  provider = provider
)

# View structured output
print(result$structured_output)

# View conversation history
print(result$conversation_history)
```

_Note:_ When using OpenRouter, specify the complete model ID.
