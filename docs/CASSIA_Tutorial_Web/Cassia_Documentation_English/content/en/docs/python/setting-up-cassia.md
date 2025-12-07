---
title: Setting Up CASSIA
---

First, install the CASSIA package using pip:

```bash
pip install CASSIA
```

## Import CASSIA

```python
import CASSIA
```

## Setting API Keys

To use LLMs like OpenAI's GPT-4, Anthropic's Claude, or models via OpenRouter, you will first need to get your API keys from the provider and then set your API keys using the `set_api_key` function.

**Note: You must set at least one API key to use CASSIA.**

**You only need to choose one provider.** OpenRouter is recommended as it provides access to multiple models.

```python
# Set API key (choose one provider)
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")  # Recommended
# CASSIA.set_api_key("your-openai-key", provider="openai")
# CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
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

## Reasoning Effort Parameter

The `reasoning` parameter controls the model's reasoning depth for compatible models.

### Syntax

```python
# Simple string format
reasoning="high"    # Maximum reasoning depth
reasoning="medium"  # Balanced (recommended for GPT-5.1)
reasoning="low"     # Minimal reasoning

# Dictionary format (Python only)
reasoning={"effort": "high"}

# Standard mode (no extended reasoning)
reasoning=None  # or omit the parameter
```

### Provider Notes

- **OpenAI via OpenRouter**: Full control with reasoning parameter. **Recommended** to avoid identity verification required by direct OpenAI API.
- **OpenAI direct**: Requires identity verification for reasoning models.
- **Anthropic Claude**: Default uses highest reasoning effort automatically.
- **Gemini**: Dynamic thinking - model decides when and how much to think.

### Recommendations

- **GPT-5.1**: Use `reasoning="medium"` - highest effort may take a very long time
- **GPT-4o**: Still excellent performance without reasoning parameter
- **Claude**: No need to set - uses optimal reasoning automatically
- **General**: Higher effort = longer processing time + higher cost

Example:
```python
# Using reasoning with GPT-5.1 via OpenRouter
CASSIA.runCASSIA_batch(
    marker=markers,
    output_name="results",
    model="openai/gpt-5.1",
    provider="openrouter",
    reasoning="medium",  # Recommended for balanced speed/quality
    tissue="brain",
    species="human"
)
```

## Custom API Providers

CASSIA supports any OpenAI-compatible API endpoint, allowing you to use custom providers like DeepSeek, local LLM servers, or other third-party services.

### Setting Up Custom Providers

To use a custom API provider, specify the full URL as the `provider` parameter:

```python
# Set API key for custom provider
CASSIA.set_api_key("your-api-key", provider="https://api.deepseek.com")

# Use in analysis
CASSIA.runCASSIA_batch(
    marker=markers,
    output_name="results",
    provider="https://api.deepseek.com",
    model="deepseek-chat",
    tissue="brain",
    species="human"
)
```

### DeepSeek Example

DeepSeek offers high-performance models at competitive prices:

1. Get your API key from [DeepSeek Platform](https://platform.deepseek.com/)
2. Set up in CASSIA:

```python
CASSIA.set_api_key("your-deepseek-key", provider="https://api.deepseek.com")

CASSIA.runCASSIA_pipeline(
    output_file_name="analysis",
    marker=markers,
    annotation_provider="https://api.deepseek.com",
    annotation_model="deepseek-chat",
    tissue="brain",
    species="human"
)
```

### Compatible Providers

Any API that follows the OpenAI chat completions format can be used, including:
- DeepSeek (`https://api.deepseek.com`)
- Local LLM servers (e.g., Ollama, vLLM)
- Other OpenAI-compatible services
