---
title: Setting Up CASSIA
---

First, ensure you have the reticulate package and the devtools package installed.

```r
install.packages("reticulate")
install.packages("devtools")
```

Next, you can install the CASSIA package from github:

```r
# Install the CASSIA package
library(devtools)
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

## Setting Up the Python Environment

CASSIA relies on Python for some of its backend processing. When you load the CASSIA package, it attempts to set up the required Python environment automatically. However, if you encounter issues, you can use the `setup_cassia_env()` function to create and configure the necessary Python environment automatically.

```r
library(CASSIA)

# Automatically set up the Python environment if needed
setup_cassia_env()
```

This function will:

- Create a new Python virtual environment named `cassia_env` (default). If virtual environment creation fails, it will attempt to create a Conda environment instead.
- Install the required Python packages: `openai`, `pandas`, `numpy`, `requests`, `anthropic`, `matplotlib`, and `seaborn`.

## Setting API Keys


To use LLMs like OpenAI's GPT-4, Anthropic's Claude, or models via OpenRouter, you will first need to get your API keys from the provider and then set your API keys using the `setLLMApiKey()` function. Normally it will take about 3 minutes to get the API key from the provider. We recommend starting with OpenRouter, as it offers access to a variety of models including open source options.

**Note: You must set at least one API key to use CASSIA.**

```r

# For OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)

# For OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# For Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)
```

- Replace `"your_api_key"` with your actual API key.
- Set `provider` to `"openai"`, `"anthropic"`, or `"openrouter"` depending on your provider.
- Setting `persist = TRUE` saves the key in your `.Renviron` file for future sessions.

## Validating API Keys

You can verify that your API keys are working correctly before running analyses:

```r
# Validate all configured providers
validate_api_keys(force_revalidate = TRUE)

# Validate a specific provider
validate_api_keys("openai", force_revalidate = TRUE)
```

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

## Smart Model Settings

CASSIA includes a smart model selection system that allows you to use simple aliases or "tiers" instead of remembering exact model version strings.

### Tier Shortcuts
You can use these shortcuts with any provider to get the appropriate model for your needs:

- `"best"`: Selects the highest performing model (e.g., `gpt-5.1`, `claude-opus-4.5`)
- `"balanced"`: Selects a good balance of performance and cost (e.g., `gpt-4o`, `claude-sonnet-4.5`)
- `"fast"`: Selects the fastest/cheapest model (e.g., `gpt-5-mini`, `claude-haiku-4.5`)

Example:
```r
# This will automatically select the best model for OpenAI (gpt-5.1)
runCASSIA_pipeline(..., model = "best", provider = "openai")
```

### Fuzzy Matching & Aliases
You can also use common names, and CASSIA will resolve them to the correct version:

- `"gpt"` -> resolves to `gpt-5.1` (for OpenAI)
- `"claude"` -> resolves to `claude-sonnet-4.5` (for Anthropic)
- `"gemini"` -> resolves to `google/gemini-2.5-flash` (for OpenRouter)

This makes your code more robust to model version updates.

## Reasoning Effort Parameter

The `reasoning` parameter controls the model's reasoning depth for compatible models.

### Syntax

```r
reasoning = "high"    # Maximum reasoning depth
reasoning = "medium"  # Balanced (recommended for GPT-5.1)
reasoning = "low"     # Minimal reasoning
reasoning = NULL      # Standard mode (or omit the parameter)
```

### Provider Notes

- **OpenAI via OpenRouter**: Full control with reasoning parameter. **Recommended** to avoid identity verification required by direct OpenAI API.
- **OpenAI direct**: Requires identity verification for reasoning models.
- **Anthropic Claude**: Default uses highest reasoning effort automatically.
- **Gemini**: Dynamic thinking - model decides when and how much to think.

### Recommendations

- **GPT-5.1**: Use `reasoning = "medium"` - highest effort may take a very long time
- **GPT-4o**: Still excellent performance without reasoning parameter
- **Claude**: No need to set - uses optimal reasoning automatically
- **General**: Higher effort = longer processing time + higher cost

Example:
```r
# Using reasoning with GPT-5.1 via OpenRouter
runCASSIA_batch(
    marker = markers,
    output_name = "results",
    model = "openai/gpt-5.1",
    provider = "openrouter",
    reasoning = "medium",  # Recommended for balanced speed/quality
    tissue = "brain",
    species = "human"
)
```

## Custom API Providers

CASSIA supports any OpenAI-compatible API endpoint, allowing you to use custom providers like DeepSeek, local LLM servers, or other third-party services.

### Setting Up Custom Providers

To use a custom API provider, specify the full URL as the `provider` parameter:

```r
# Set API key for custom provider
setLLMApiKey("your-api-key", provider = "https://api.deepseek.com", persist = TRUE)

# Use in analysis
runCASSIA_batch(
    marker = markers,
    output_name = "results",
    provider = "https://api.deepseek.com",
    model = "deepseek-chat",
    tissue = "brain",
    species = "human"
)
```

### DeepSeek Example

DeepSeek offers high-performance models at competitive prices:

1. Get your API key from [DeepSeek Platform](https://platform.deepseek.com/)
2. Set up in CASSIA:

```r
setLLMApiKey("your-deepseek-key", provider = "https://api.deepseek.com", persist = TRUE)

runCASSIA_pipeline(
    output_file_name = "analysis",
    marker = markers,
    annotation_provider = "https://api.deepseek.com",
    annotation_model = "deepseek-chat",
    tissue = "brain",
    species = "human"
)
```

### Compatible Providers

Any API that follows the OpenAI chat completions format can be used, including:
- DeepSeek (`https://api.deepseek.com`)
- Local LLM servers (e.g., Ollama, LM Studio, vLLM)
- Other OpenAI-compatible services

### Local LLMs (Ollama, LM Studio)

For complete data privacy and zero API costs, you can run LLMs locally. CASSIA supports any OpenAI-compatible local server.

**No API key required for localhost URLs.**

#### Ollama Setup

1. Install Ollama from [ollama.ai](https://ollama.ai)
2. Pull a model: `ollama pull gpt-oss:20b`
3. Ollama runs automatically on `http://localhost:11434`

#### Usage

```r
runCASSIA_batch(
    marker = markers,
    output_name = "results",
    provider = "http://localhost:11434/v1",
    model = "gpt-oss:20b",
    tissue = "brain",
    species = "human"
)
```

#### LM Studio Setup

1. Download LM Studio from [lmstudio.ai](https://lmstudio.ai)
2. Load a model and start the local server
3. Default port: `http://localhost:1234/v1`
