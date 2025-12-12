import os
import json
import requests
from typing import Dict, Any, Optional

# Import CASSIA logger for actionable error messages
try:
    from .logging_config import get_logger
except ImportError:
    from logging_config import get_logger

logger = get_logger(__name__)

# Import model settings for automatic model name resolution
try:
    from .model_settings import resolve_model_name
except ImportError:
    # Fallback function if model_settings is not available
    def resolve_model_name(model_name: str, provider: str = None):
        return model_name, provider or "openrouter"


def _handle_api_error(exc: Exception, provider: str, model: str) -> None:
    """
    Log actionable error messages for common API failures.

    Args:
        exc: The exception that was raised
        provider: The LLM provider name
        model: The model name being used
    """
    error_str = str(exc).lower()

    # Authentication errors (401)
    if "401" in str(exc) or "unauthorized" in error_str or "invalid api key" in error_str:
        logger.error(
            f"Authentication failed for {provider}. "
            f"Please check your API key is valid. "
            f"Set it with: CASSIA.set_api_key('{provider}', 'your-key')"
        )

    # Rate limit errors (429)
    elif "429" in str(exc) or "rate limit" in error_str or "too many requests" in error_str:
        logger.error(
            f"Rate limit exceeded for {provider}. "
            f"Wait a few minutes and try again, or use a different model. "
            f"Consider reducing max_workers in batch processing."
        )

    # Timeout errors
    elif "timeout" in error_str or "timed out" in error_str:
        logger.error(
            f"Request timed out for {provider}. "
            f"The model may be overloaded. Try again or use a faster model like 'gemini-flash'."
        )

    # Model not found errors (404)
    elif "404" in str(exc) or "not found" in error_str or "does not exist" in error_str:
        logger.error(
            f"Model '{model}' not found for {provider}. "
            f"Run CASSIA.print_available_models('{provider}') to see available models. "
            f"Or try a common model like 'gpt-4o' or 'claude-sonnet-4-5'."
        )

    # Insufficient quota/credits
    elif "quota" in error_str or "insufficient" in error_str or "credit" in error_str:
        logger.error(
            f"Insufficient API credits for {provider}. "
            f"Please check your account balance and billing settings."
        )

    # Context length exceeded
    elif "context" in error_str and ("length" in error_str or "limit" in error_str or "too long" in error_str):
        logger.error(
            f"Input too long for model '{model}'. "
            f"Try reducing n_genes parameter or use a model with larger context window."
        )

    # Generic error with details
    else:
        logger.error(
            f"API call to {provider} failed with model '{model}'. "
            f"Error: {exc}"
        )

def call_llm(
    prompt: str,
    provider: str = "openai",
    model: str = None,
    api_key: Optional[str] = None,
    temperature: float = 0.7,
    max_tokens: int = 4096,
    system_prompt: Optional[str] = None,
    additional_params: Optional[Dict[str, Any]] = None,
    reasoning: Optional[Dict[str, Any]] = None
) -> str:
    """
    Call an LLM from various providers and return the generated text.

    Args:
        prompt: The user prompt to send to the LLM
        provider: One of "openai", "anthropic", or "openrouter"
        model: Specific model from the provider to use (e.g., "gpt-4" for OpenAI)
        api_key: API key for the provider (if None, gets from environment)
        temperature: Sampling temperature (0-1)
        max_tokens: Maximum tokens to generate
        system_prompt: Optional system prompt for providers that support it
        additional_params: Additional parameters to pass to the provider's API
        reasoning: Optional reasoning configuration for models that support it.
            Controls how much the model "thinks" before responding.
            Options:
            - effort: "high", "medium", "low" (OpenAI/Anthropic/OpenRouter)
            Example: {"effort": "high"} or {"effort": "medium"}

            Provider-specific behavior:
            - OpenAI: Uses Responses API with reasoning parameter (GPT-5 series)
            - Anthropic: Uses beta.messages.create with effort parameter (Claude Opus 4.5)
            - OpenRouter: Passes reasoning to chat completions endpoint

    Returns:
        str: The generated text response
    """
    provider = provider.lower()
    additional_params = additional_params or {}
    
    # Resolve model name using model settings if model is provided
    if model:
        try:
            resolved_model, resolved_provider = resolve_model_name(model, provider)
            # Use the resolved model (provider should stay the same)
            model = resolved_model
        except Exception as e:
            # If resolution fails, warn and continue with original model name
            print(f"Warning: Model resolution failed for '{model}': {str(e)}. Using original name.")
    
    # Default models for each provider if not specified
    default_models = {
        "openai": "gpt-4o",
        "anthropic": "claude-3-5-sonnet-latest",
        "openrouter": "google/gemini-2.5-flash",
    }
    
    # Use default model if not specified
    if not model:
        model = default_models.get(provider)
        if not model:
            raise ValueError(f"No model specified and no default available for provider: {provider}")
    
    # Get API key from environment if not provided
    if not api_key:
        env_var_names = {
            "openai": "OPENAI_API_KEY",
            "anthropic": "ANTHROPIC_API_KEY",
            "openrouter": "OPENROUTER_API_KEY",
        }
        env_var = env_var_names.get(provider)
        if env_var:
            api_key = os.environ.get(env_var)
            if not api_key:
                raise ValueError(f"API key not provided and {env_var} not found in environment")
    
    # Prepare messages format
    messages = []
    if system_prompt:
        messages.append({"role": "system", "content": system_prompt})
    
    messages.append({"role": "user", "content": prompt})
    
    # OpenAI API call
    if provider == "openai":
        try:
            import openai
        except ImportError:
            raise ImportError("Please install openai package: pip install openai")

        client = openai.OpenAI(api_key=api_key)

        # Handle message history from additional_params (for conversation history)
        params_copy = additional_params.copy() if additional_params else {}
        api_messages = messages.copy()

        if 'messages' in params_copy:
            # Use the full conversation history from additional_params
            api_messages = params_copy.pop('messages')
            # Only add system prompt if not already in history
            if system_prompt and not any(msg.get('role') == 'system' for msg in api_messages):
                api_messages.insert(0, {"role": "system", "content": system_prompt})

        # Use Responses API when reasoning is specified (for GPT-5 reasoning models)
        if reasoning:
            try:
                # Convert messages to input format for Responses API
                # Note: Responses API uses "input" not "messages", and "developer" not "system"
                input_messages = []
                for msg in api_messages:
                    role = "developer" if msg["role"] == "system" else msg["role"]
                    input_messages.append({"role": role, "content": msg["content"]})

                response = client.responses.create(
                    model=model,
                    input=input_messages,
                    reasoning=reasoning,
                    **params_copy
                )
                return response.output_text
            except Exception as e:
                _handle_api_error(e, provider, model)
                raise

        # Standard Chat Completions API (no reasoning)
        # GPT-5 models and o-series require max_completion_tokens instead of max_tokens
        model_lower = model.lower() if model else ""
        uses_max_completion_tokens = any(m in model_lower for m in ["gpt-5", "gpt5", "o1", "o3", "o4"])

        try:
            if uses_max_completion_tokens:
                response = client.chat.completions.create(
                    model=model,
                    messages=api_messages,
                    temperature=temperature,
                    max_completion_tokens=max_tokens,
                    **params_copy
                )
            else:
                response = client.chat.completions.create(
                    model=model,
                    messages=api_messages,
                    temperature=temperature,
                    max_tokens=max_tokens,
                    **params_copy
                )
            return response.choices[0].message.content
        except Exception as e:
            _handle_api_error(e, provider, model)
            raise
    
    # Custom OpenAI-compatible API call (base_url as provider)
    elif provider.startswith("http"):
        try:
            import openai
        except ImportError:
            raise ImportError("Please install openai package: pip install openai")
        custom_api_key = api_key or os.environ.get("CUSTOMIZED_API_KEY")

        # For localhost URLs, API key is optional (local LLMs like Ollama don't need auth)
        is_localhost = any(x in provider.lower() for x in ["localhost", "127.0.0.1"])
        if not custom_api_key:
            if is_localhost:
                custom_api_key = "ollama"  # Placeholder for local LLMs
            else:
                raise ValueError("API key not provided and CUSTOMIZED_API_KEY not found in environment")

        client = openai.OpenAI(api_key=custom_api_key, base_url=provider)

        # Handle message history properly
        api_messages = messages.copy()

        # If additional_params contains message history, merge it properly
        if 'messages' in additional_params:
            # Use the full conversation history from additional_params instead
            history_messages = additional_params.pop('messages')
            api_messages = history_messages

            # Only add system prompt if it's not already in the history
            if system_prompt and not any(msg.get('role') == 'system' for msg in api_messages):
                api_messages.insert(0, {"role": "system", "content": system_prompt})

        # Call the API with the proper message history
        try:
            response = client.chat.completions.create(
                model=model,
                messages=api_messages,
                temperature=temperature,
                max_tokens=max_tokens,
                **additional_params
            )
            return response.choices[0].message.content
        except Exception as e:
            _handle_api_error(e, "custom", model)
            raise
    
    # Anthropic API call
    elif provider == "anthropic":
        try:
            import anthropic
        except ImportError:
            raise ImportError("Please install anthropic package: pip install anthropic")

        client = anthropic.Anthropic(api_key=api_key)

        # Handle message history from additional_params (for conversation history)
        params_copy = additional_params.copy() if additional_params else {}

        if 'messages' in params_copy:
            # Use the full conversation history from additional_params
            history_messages = params_copy.pop('messages')

            # Anthropic doesn't accept "system" role in messages - filter it out
            # and use the system_prompt parameter instead
            api_messages = [
                msg for msg in history_messages
                if msg.get('role') != 'system'
            ]
        else:
            # No history, just use the single prompt
            api_messages = [{"role": "user", "content": prompt}]

        # Create the message params
        message_params = {
            "model": model,
            "max_tokens": max_tokens,
            "temperature": temperature,
            "messages": api_messages
        }

        # Add system prompt if provided (Anthropic uses separate system parameter)
        if system_prompt:
            message_params["system"] = system_prompt

        # Add any remaining additional parameters (skip model to prevent override)
        for key, value in params_copy.items():
            if key != "model":
                message_params[key] = value

        # Use beta API when effort/reasoning is specified (for Claude Opus 4.5)
        if reasoning and reasoning.get("effort"):
            try:
                # Use beta.messages.create with effort parameter
                response = client.beta.messages.create(
                    betas=["effort-2025-11-24"],
                    output_config={"effort": reasoning["effort"]},
                    **message_params
                )

                # Extract the text content from the response
                if hasattr(response, 'content') and len(response.content) > 0:
                    content_block = response.content[0]
                    if hasattr(content_block, 'text'):
                        return content_block.text
                    elif isinstance(content_block, dict) and 'text' in content_block:
                        return content_block['text']
                    else:
                        return str(response.content)
                else:
                    return "No content returned from Anthropic API"
            except Exception as e:
                _handle_api_error(e, provider, model)
                raise

        # Standard API call (no effort/reasoning)
        try:
            response = client.messages.create(**message_params)

            # Extract the text content from the response
            if hasattr(response, 'content') and len(response.content) > 0:
                content_block = response.content[0]
                if hasattr(content_block, 'text'):
                    return content_block.text
                elif isinstance(content_block, dict) and 'text' in content_block:
                    return content_block['text']
                else:
                    return str(response.content)
            else:
                return "No content returned from Anthropic API"
        except Exception as e:
            _handle_api_error(e, provider, model)
            raise
    
    # OpenRouter API call
    elif provider == "openrouter":
        url = "https://openrouter.ai/api/v1/chat/completions"

        headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }

        # Handle message history properly (similar to custom HTTP provider)
        api_messages = messages.copy()
        params_copy = additional_params.copy() if additional_params else {}

        # If additional_params contains message history, use it instead of the single prompt
        if 'messages' in params_copy:
            api_messages = params_copy.pop('messages')
            # Only add system prompt if not already in history
            if system_prompt and not any(msg.get('role') == 'system' for msg in api_messages):
                api_messages.insert(0, {"role": "system", "content": system_prompt})

        # GPT-5 and reasoning models require max_completion_tokens instead of max_tokens
        model_lower = model.lower() if model else ""
        uses_max_completion_tokens = any(m in model_lower for m in ["gpt-5", "gpt5", "o1", "o3", "o4"])

        # Auto-default reasoning to "medium" for GPT-5 series models if not explicitly set
        # Covers: gpt-5, gpt5, gpt-5.1, openai/gpt-5.1, etc.
        if reasoning is None and ("gpt-5" in model_lower or "gpt5" in model_lower):
            reasoning = {"effort": "medium"}

        data = {
            **params_copy,
            "model": model,
            "messages": api_messages,
            "temperature": temperature,
        }

        # Use appropriate token parameter based on model
        if uses_max_completion_tokens:
            data["max_completion_tokens"] = max_tokens
        else:
            data["max_tokens"] = max_tokens

        # Add reasoning configuration if provided (for models that support it)
        if reasoning:
            data["reasoning"] = reasoning

        try:
            response = requests.post(url, headers=headers, data=json.dumps(data), timeout=180)
            response.raise_for_status()
            return response.json()["choices"][0]["message"]["content"]
        except Exception as e:
            _handle_api_error(e, provider, model)
            raise

    else:
        raise ValueError(f"Unsupported provider: {provider}") 