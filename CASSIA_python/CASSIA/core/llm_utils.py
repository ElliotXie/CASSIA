import os
import json
import requests
import time
from typing import Dict, Any, Optional

# Import CASSIA logger for actionable error messages
try:
    from .logging_config import get_logger
except ImportError:
    from logging_config import get_logger

logger = get_logger(__name__)

_LLM_USAGE_LOG = []


def reset_llm_usage_log() -> None:
    """Clear the in-process LLM usage/cost log."""
    _LLM_USAGE_LOG.clear()


def get_llm_usage_log():
    """Return a copy of the in-process LLM usage/cost log."""
    return [entry.copy() for entry in _LLM_USAGE_LOG]


def get_llm_usage_summary(reset: bool = False) -> Dict[str, Any]:
    """Summarize recorded LLM usage for the current Python process.

    Cost units are provider-native. OpenRouter reports cost in credits in the
    response ``usage`` object; OpenAI/Anthropic SDK responses usually expose
    token counts but not a normalized dollar cost.
    """
    summary = {
        "requests": len(_LLM_USAGE_LOG),
        "prompt_tokens": 0,
        "completion_tokens": 0,
        "reasoning_tokens": 0,
        "cached_tokens": 0,
        "total_tokens": 0,
        "cost": 0.0,
        "by_model": {},
    }

    for entry in _LLM_USAGE_LOG:
        model_key = f"{entry.get('provider', 'unknown')}::{entry.get('model', 'unknown')}"
        model_summary = summary["by_model"].setdefault(
            model_key,
            {
                "requests": 0,
                "prompt_tokens": 0,
                "completion_tokens": 0,
                "reasoning_tokens": 0,
                "cached_tokens": 0,
                "total_tokens": 0,
                "cost": 0.0,
            },
        )
        model_summary["requests"] += 1

        for key in ["prompt_tokens", "completion_tokens", "reasoning_tokens", "cached_tokens", "total_tokens"]:
            value = entry.get(key) or 0
            summary[key] += value
            model_summary[key] += value

        cost = entry.get("cost")
        if cost is not None:
            summary["cost"] += cost
            model_summary["cost"] += cost

    if reset:
        reset_llm_usage_log()

    return summary


def _usage_obj_to_dict(usage) -> Dict[str, Any]:
    if usage is None:
        return {}
    if isinstance(usage, dict):
        return usage
    if hasattr(usage, "model_dump"):
        return usage.model_dump()
    if hasattr(usage, "to_dict"):
        return usage.to_dict()

    result = {}
    for key in [
        "prompt_tokens",
        "completion_tokens",
        "total_tokens",
        "input_tokens",
        "output_tokens",
        "cost",
        "prompt_tokens_details",
        "completion_tokens_details",
        "cache_creation_input_tokens",
        "cache_read_input_tokens",
    ]:
        if hasattr(usage, key):
            result[key] = getattr(usage, key)
    return result


def _nested_get(mapping: Dict[str, Any], *path, default=None):
    current = mapping
    for key in path:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
    return current if current is not None else default


def _record_llm_usage(
    provider: str,
    model: str,
    usage=None,
    response_id: Optional[str] = None,
    response_model: Optional[str] = None,
) -> None:
    usage_dict = _usage_obj_to_dict(usage)

    prompt_tokens = usage_dict.get("prompt_tokens")
    if prompt_tokens is None:
        prompt_tokens = usage_dict.get("input_tokens")

    completion_tokens = usage_dict.get("completion_tokens")
    if completion_tokens is None:
        completion_tokens = usage_dict.get("output_tokens")

    total_tokens = usage_dict.get("total_tokens")
    if total_tokens is None and (prompt_tokens is not None or completion_tokens is not None):
        total_tokens = (prompt_tokens or 0) + (completion_tokens or 0)

    reasoning_tokens = _nested_get(
        usage_dict,
        "completion_tokens_details",
        "reasoning_tokens",
        default=usage_dict.get("reasoning_tokens"),
    )

    cached_tokens = _nested_get(
        usage_dict,
        "prompt_tokens_details",
        "cached_tokens",
        default=usage_dict.get("cache_read_input_tokens"),
    )

    entry = {
        "timestamp": time.time(),
        "provider": provider,
        "model": response_model or model,
        "requested_model": model,
        "response_id": response_id,
        "prompt_tokens": prompt_tokens or 0,
        "completion_tokens": completion_tokens or 0,
        "reasoning_tokens": reasoning_tokens or 0,
        "cached_tokens": cached_tokens or 0,
        "total_tokens": total_tokens or 0,
        "cost": usage_dict.get("cost"),
        "usage": usage_dict,
    }

    _LLM_USAGE_LOG.append(entry)

# Import model settings for automatic model name resolution
try:
    from .model_settings import resolve_model_name
except ImportError:
    # Fallback function if model_settings is not available
    def resolve_model_name(model_name: str, provider: str = None):
        return model_name, provider or "openrouter"

# Import proxy configuration
try:
    from .proxy_config import resolve_proxy_url
except ImportError:
    def resolve_proxy_url(provider):
        return None


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
            f"Or try a current model like 'gpt-5.6-terra' or 'claude-sonnet-5'."
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

def _extract_anthropic_text(response) -> str:
    """Return the first text block, skipping any leading thinking blocks."""
    content = getattr(response, "content", None) or []
    for block in content:
        text = getattr(block, "text", None)
        if text is None and isinstance(block, dict):
            text = block.get("text")
        if text is not None:
            return text
    return str(content) if content else "No content returned from Anthropic API"


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
        model: Specific model from the provider to use (e.g., "gpt-5.6-terra" for OpenAI)
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
            - OpenAI: Uses Responses API with reasoning parameter (GPT-5.6/GPT-6 series)
            - Anthropic: Uses messages API with effort parameter when requested
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
    
    # Use default model from model_settings.json if not specified
    if not model:
        try:
            from .model_settings import get_model_settings
            providers = get_model_settings().settings.get("providers", {})
            provider_settings = providers.get(provider, {})
            model = provider_settings.get("balanced") or provider_settings.get("recommended")
        except Exception:
            pass
        if not model:
            raise ValueError(f"No model specified and no default available for provider: {provider}")
    
    # Get API key from environment if not provided
    if not api_key:
        # Custom OpenAI-compatible endpoint (provider is a base URL):
        # read CUSTOMIZED_API_KEY set via CASSIA.set_api_key(key, provider=url)
        if provider.startswith("http"):
            api_key = os.environ.get("CUSTOMIZED_API_KEY")

    if not api_key and not provider.startswith("http"):
        env_var_names = {
            "openai": "OPENAI_API_KEY",
            "anthropic": "ANTHROPIC_API_KEY",
            "openrouter": "OPENROUTER_API_KEY",
            "google": "GOOGLE_API_KEY",
            "together": "TOGETHER_API_KEY",
        }
        env_var = env_var_names.get(provider)
        if env_var:
            api_key = os.environ.get(env_var)

        # If no API key in environment, try free API fallback for supported providers
        if not api_key:
            try:
                from .free_api import get_free_api_key, FREE_API_PROVIDERS, _thread_local

                if provider in FREE_API_PROVIDERS:
                    # Get cluster count from thread-local context (set by runCASSIA_batch)
                    num_clusters = getattr(_thread_local, 'num_clusters', 1)
                    api_key, error = get_free_api_key(provider, num_clusters)

                    if error:
                        raise ValueError(
                            f"No API key configured and free API access unavailable:\n"
                            f"{error}"
                        )
                else:
                    env_var_display = env_var if env_var else f"{provider.upper()}_API_KEY"
                    raise ValueError(
                        f"API key not provided and {env_var_display} not found in environment.\n"
                        f"Free API access is only available for: {', '.join(sorted(FREE_API_PROVIDERS))}.\n"
                        f"Set your API key with: CASSIA.set_api_key('{provider}', 'your-key')"
                    )
            except ImportError:
                # free_api module not available, raise original error
                env_var_display = env_var if env_var else f"{provider.upper()}_API_KEY"
                raise ValueError(f"API key not provided and {env_var_display} not found in environment")
    
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

        # Use proxy base URL if active (e.g., Cloudflare Worker for China users)
        proxy_url = resolve_proxy_url("openai")
        client_kwargs = {"api_key": api_key}
        if proxy_url:
            client_kwargs["base_url"] = proxy_url
        client = openai.OpenAI(**client_kwargs)

        # Handle message history from additional_params (for conversation history)
        params_copy = additional_params.copy() if additional_params else {}
        api_messages = messages.copy()

        if 'messages' in params_copy:
            # Use the full conversation history from additional_params
            api_messages = params_copy.pop('messages')
            # Only add system prompt if not already in history
            if system_prompt and not any(msg.get('role') == 'system' for msg in api_messages):
                api_messages.insert(0, {"role": "system", "content": system_prompt})

        model_lower = model.lower() if model else ""
        if reasoning is None and ("gpt-5.6" in model_lower or "gpt5.6" in model_lower):
            reasoning = {"effort": "medium"}
        elif reasoning is None and ("gpt-6" in model_lower or "gpt6" in model_lower):
            reasoning = {"effort": "low"}

        # Use Responses API when reasoning is specified.
        if reasoning:
            try:
                # Convert messages to input format for Responses API
                # Note: Responses API uses "input" not "messages", and "developer" not "system"
                input_messages = []
                for msg in api_messages:
                    role = "developer" if msg["role"] == "system" else msg["role"]
                    input_messages.append({"role": role, "content": msg["content"]})

                response_params = {
                    "model": model,
                    "input": input_messages,
                    "reasoning": reasoning,
                    "max_output_tokens": max_tokens,
                    **params_copy,
                }
                response = client.responses.create(**response_params)
                _record_llm_usage(
                    provider=provider,
                    model=model,
                    usage=getattr(response, "usage", None),
                    response_id=getattr(response, "id", None),
                    response_model=getattr(response, "model", None),
                )
                return response.output_text
            except Exception as e:
                _handle_api_error(e, provider, model)
                raise

        # Standard Chat Completions API (no reasoning)
        # GPT-5/GPT-6 models and o-series require max_completion_tokens.
        uses_max_completion_tokens = any(m in model_lower for m in ["gpt-5", "gpt5", "gpt-6", "gpt6", "o1", "o3", "o4"])

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
            _record_llm_usage(
                provider=provider,
                model=model,
                usage=getattr(response, "usage", None),
                response_id=getattr(response, "id", None),
                response_model=getattr(response, "model", None),
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
            _record_llm_usage(
                provider="custom",
                model=model,
                usage=getattr(response, "usage", None),
                response_id=getattr(response, "id", None),
                response_model=getattr(response, "model", None),
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

        # Use proxy base URL if active (e.g., Cloudflare Worker for China users)
        proxy_url = resolve_proxy_url("anthropic")
        client_kwargs = {"api_key": api_key}
        if proxy_url:
            client_kwargs["base_url"] = proxy_url
        client = anthropic.Anthropic(**client_kwargs)

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
            "messages": api_messages
        }

        model_lower = model.lower() if model else ""
        if "claude-sonnet-5" not in model_lower and "claude-opus-5" not in model_lower:
            message_params["temperature"] = temperature

        # Add system prompt if provided (Anthropic uses separate system parameter)
        if system_prompt:
            message_params["system"] = system_prompt

        # Add any remaining additional parameters (skip model to prevent override)
        for key, value in params_copy.items():
            if key != "model":
                message_params[key] = value

        # Claude 5 exposes output_config on the standard Messages API; older
        # models keep the beta path for backwards compatibility.
        if reasoning and reasoning.get("effort"):
            try:
                if "claude-sonnet-5" in model_lower or "claude-opus-5" in model_lower:
                    message_params["output_config"] = {"effort": reasoning["effort"]}
                    response = client.messages.create(**message_params)
                else:
                    response = client.beta.messages.create(
                        betas=["effort-2025-11-24"],
                        output_config={"effort": reasoning["effort"]},
                        **message_params
                    )
                _record_llm_usage(
                    provider=provider,
                    model=model,
                    usage=getattr(response, "usage", None),
                    response_id=getattr(response, "id", None),
                    response_model=getattr(response, "model", None),
                )

                return _extract_anthropic_text(response)
            except Exception as e:
                _handle_api_error(e, provider, model)
                raise

        # Standard API call (no effort/reasoning)
        try:
            response = client.messages.create(**message_params)
            _record_llm_usage(
                provider=provider,
                model=model,
                usage=getattr(response, "usage", None),
                response_id=getattr(response, "id", None),
                response_model=getattr(response, "model", None),
            )

            return _extract_anthropic_text(response)
        except Exception as e:
            _handle_api_error(e, provider, model)
            raise
    
    # OpenRouter API call
    elif provider == "openrouter":
        # Use proxy URL if active (e.g., Cloudflare Worker for China users)
        proxy_url = resolve_proxy_url("openrouter")
        url = f"{proxy_url}/chat/completions" if proxy_url else "https://openrouter.ai/api/v1/chat/completions"

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

        # GPT-5/GPT-6 and reasoning models require max_completion_tokens.
        model_lower = model.lower() if model else ""
        uses_max_completion_tokens = any(m in model_lower for m in ["gpt-5", "gpt5", "gpt-6", "gpt6", "o1", "o3", "o4"])

        # Auto-default reasoning to "medium" for GPT-5 series models if not explicitly set
        # Covers: gpt-5, gpt5, gpt-5.6-terra, openai/gpt-5.6-terra, etc.
        if reasoning is None and ("gpt-5" in model_lower or "gpt5" in model_lower):
            reasoning = {"effort": "medium"}
        elif reasoning is None and ("gpt-6" in model_lower or "gpt6" in model_lower):
            reasoning = {"effort": "low"}

        # Kimi K2.6 can spend the full OpenRouter output budget in reasoning
        # and return message.content=None unless reasoning is explicitly
        # disabled. CASSIA expects plain text for downstream parsers.
        if reasoning is None and "kimi-k2.6" in model_lower:
            reasoning = {"effort": "none", "exclude": True}

        data = {
            **params_copy,
            "model": model,
            "messages": api_messages,
        }

        omits_sampling = any(m in model_lower for m in ["gpt-6-astra", "claude-sonnet-5", "claude-opus-5"])
        if not omits_sampling:
            data["temperature"] = temperature

        # Use appropriate token parameter based on model
        if uses_max_completion_tokens:
            data["max_completion_tokens"] = max_tokens
        else:
            data["max_tokens"] = max_tokens

        # Add reasoning configuration only for models that support it.
        # OpenRouter's Anthropic Claude models do not accept the reasoning parameter.
        if reasoning:
            supports_reasoning = any(m in model_lower for m in ["gpt-5", "gpt5", "gpt-6", "gpt6", "o1", "o3", "o4", "kimi-k2.6"])
            if supports_reasoning:
                data["reasoning"] = reasoning

        try:
            response = requests.post(url, headers=headers, data=json.dumps(data), timeout=180)
            response.raise_for_status()
            response_json = response.json()
            _record_llm_usage(
                provider=provider,
                model=model,
                usage=response_json.get("usage"),
                response_id=response_json.get("id"),
                response_model=response_json.get("model"),
            )
            message = response_json["choices"][0]["message"]
            content = message.get("content")
            if content is not None:
                return content
            reasoning_text = message.get("reasoning")
            if reasoning_text:
                return reasoning_text
            return ""
        except Exception as e:
            _handle_api_error(e, provider, model)
            raise

    else:
        raise ValueError(f"Unsupported provider: {provider}")
