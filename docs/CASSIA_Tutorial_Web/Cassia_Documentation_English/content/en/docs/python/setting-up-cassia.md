---
title: Setting Up CASSIA
---

First, install the CASSIA package using pip:

```python
!pip install CASSIA
```

## Import CASSIA

```python
import CASSIA
```

## Setting API Keys

To use LLMs like OpenAI's GPT-4, Anthropic's Claude, or models via OpenRouter, you will first need to get your API keys from the provider and then set your API keys using the `set_api_key` function.

**Note: You must set at least one API key to use CASSIA.**

```python
# Set API keys
CASSIA.set_api_key("your-openai-key", provider="openai")
CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")
```

- Replace `"your-key"` with your actual API key.
- Set `provider` to `"openai"`, `"anthropic"`, or `"openrouter"` depending on your provider.

## How to Select Models and Providers

There are three providers to choose from: `openrouter`, `openai`, and `anthropic`. Each provider has its own models and pricing.
**Note that the model name must be set exactly as shown below or the model will not be found.**

### OpenRouter

OpenRouter is a platform that offers access to almost all the models supported by major providers. It is recommended to use OpenRouter due to its higher rate limits and access to a variety of models including open source options.

- `anthropic/claude-sonnet-4.5`: Current default high-performance model.
- `openai/gpt-5.1`: Balanced option.
- `google/gemini-2.5-flash`: One of the best-performed low-cost models.

### OpenAI

- `gpt-5.1`: High-performance model.

### Anthropic

- `claude-sonnet-4.5`: High-performance model.

## Smart Model Settings (Recommended)

CASSIA includes a smart model selection system that allows you to use simple aliases or "tiers" instead of remembering exact model version strings. This makes your code more robust to model version updates.

### Tier Shortcuts
You can use these shortcuts with any provider to get the appropriate model for your needs:

- `"best"`: Selects the highest performing model (e.g., `gpt-5.1`, `claude-opus-4.5`)
- `"balanced"`: Selects a good balance of performance and cost (e.g., `gpt-4o`, `claude-sonnet-4.5`)
- `"fast"`: Selects the fastest/cheapest model (e.g., `gpt-5-mini`, `claude-haiku-4.5`)

Example:
```python
# This will automatically select the best model for OpenAI
CASSIA.runCASSIA_pipeline(..., model = "best", provider = "openai")
```

### Fuzzy Matching & Aliases
You can also use common names, and CASSIA will resolve them to the correct version:

- `"gpt"` -> resolves to `gpt-5.1` (for OpenAI)
- `"claude"` -> resolves to `claude-sonnet-4.5` (for Anthropic)
- `"gemini"` -> resolves to `google/gemini-2.5-flash` (for OpenRouter)

## Loading Example Markers

You can load example marker data provided with the package:

```python
processed_markers = CASSIA.loadmarker(marker_type="processed")
unprocessed_markers = CASSIA.loadmarker(marker_type="unprocessed")
subcluster_results = CASSIA.loadmarker(marker_type="subcluster_results")

# List available marker sets
available_markers = CASSIA.list_available_markers()
print(available_markers) 
```
