"""
API Key Configuration File
--------------------------
This file sets API keys as environment variables for the current Python session.
DO NOT commit this file to git - it contains sensitive credentials.

Usage:
    import set_api_keys  # Run this at the start of your script/notebook
"""

import os

# =============================================================================
# ADD YOUR API KEYS BELOW
# =============================================================================

# Anthropic API Key (for Claude models)
ANTHROPIC_API_KEY = "your-anthropic-api-key-here"

# OpenAI API Key
OPENAI_API_KEY = "your-openai-api-key-here"

# OpenRouter API Key
OPENROUTER_API_KEY = "your-openrouter-api-key-here"

# Custom API Key (if needed)
CUSTOMIZED_API_KEY = "your-custom-key-here"

# =============================================================================
# SET ENVIRONMENT VARIABLES (runs automatically on import)
# =============================================================================

def set_keys():
    """Set all API keys as environment variables."""
    keys = {
        "ANTHROPIC_API_KEY": ANTHROPIC_API_KEY,
        "OPENAI_API_KEY": OPENAI_API_KEY,
        "OPENROUTER_API_KEY": OPENROUTER_API_KEY,
        "CUSTOMIZED_API_KEY": CUSTOMIZED_API_KEY,
    }

    for key_name, key_value in keys.items():
        if key_value and not key_value.startswith("your-"):
            os.environ[key_name] = key_value
            print(f"Set {key_name}")

# Auto-run when imported
set_keys()


# Optionally print confirmation (uncomment for debugging)
# print("API keys configured:", list(k for k, v in os.environ.items() if k.endswith("_API_KEY")))