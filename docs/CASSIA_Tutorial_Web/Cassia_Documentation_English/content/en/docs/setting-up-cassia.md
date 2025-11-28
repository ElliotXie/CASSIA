---
title: Setting Up CASSIA
---

First, ensure you have the reticulate package and the devtools package installed.

\`\`\`r
install.packages("reticulate")
install.packages("devtools")
\`\`\`

Next, you can install the CASSIA package from github:

\`\`\`r
# Install the CASSIA package
library(devtools)
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
\`\`\`

## Setting Up the Python Environment

CASSIA relies on Python for some of its backend processing. When you load the CASSIA package, it attempts to set up the required Python environment automatically. However, if you encounter issues, you can use the `setup_cassia_env()` function to create and configure the necessary Python environment automatically.

\`\`\`r
library(CASSIA)

# Automatically set up the Python environment if needed
setup_cassia_env()
\`\`\`

This function will:

- Create a new Conda environment named `cassia_env` if it doesn't already exist.
- Install the required Python packages: `openai`, `pandas`, `numpy`, `scikit-learn`, `requests`, and `anthropic`.

## Setting API Keys


To use LLMs like OpenAI's GPT-4, Anthropic's Claude, or models via OpenRouter, you will first need to get your API keys from the provider and then set your API keys using the `setLLMApiKey()` function. Normally it will take about 3 minutes to get the API key from the provider. We recommend starting with OpenRouter, as it offers access to a variety of models including open source options.

**Note: You must set at least one API key to use CASSIA.**

\`\`\`r

# For OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)

# For OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# For Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)
\`\`\`

- Replace `"your_api_key"` with your actual API key.
- Set `provider` to `"openai"`, `"anthropic"`, or `"openrouter"` depending on your provider.
- Setting `persist = TRUE` saves the key in your `.Renviron` file for future sessions.



## How to Select Models and Providers

There are three providers to choose from: `openrouter`, `openai`, and `anthropic`. Each provider has its own models and pricing.
**Note that the model name must be set exactly as shown below or the model will not be found.**

### OpenRouter

OpenRouter is a platform that offers access to almost all the models supported by major providers. It is recommended to use OpenRouter due to its higher rate limits and access to a variety of models including open source options.

- `google/gemini-2.5-flash-preview`: One of the best-performed low-cost models, comparable with models like gpt-4o (Most recommended)
- `deepseek/deepseek-chat-v3-0324`: One of the best-performed open-source models, turns to gives very detailed annotations (Recommended)
- `deepseek/deepseek-chat-v3-0324:free`: Free but slower
- `anthropic/claude-3.7-sonnet`
- `openai/gpt-4o`

### OpenAI

- `gpt-4o`: Used in the benchmark

### Anthropic

- `claude-3-7-sonnet-latest`: The latest High-performance model
