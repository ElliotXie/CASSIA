"""
API Key Validation Module

This module provides functionality to validate API keys for various LLM providers
(OpenAI, Anthropic, OpenRouter) by making minimal test calls. It includes an
intelligent caching system to avoid redundant validations and minimize costs.

Key features:
- Hash-based caching: Only validates each unique key once per session
- Thread-safe: Safe for concurrent use in batch operations
- Minimal cost: Uses cheapest models for testing (~$0.000001 per validation)
- Clear error messages: Actionable guidance for fixing issues
"""

import os
import hashlib
import threading
import time
from typing import Optional, Dict, Union, Tuple

from .logging_config import get_logger
from .exceptions import APIValidationError
from .llm_utils import call_llm

# Module logger
logger = get_logger(__name__)

# Validation models (cheapest options for each provider)
VALIDATION_MODELS = {
    "openai": "gpt-5-nano",
    "anthropic": "claude-3-haiku-20240307",
    "openrouter": "google/gemini-2.5-flash-lite"
}

# Test prompt configuration
VALIDATION_PROMPT = "Reply: OK"
VALIDATION_MAX_TOKENS = 10
VALIDATION_TEMPERATURE = 0.0

# Module-level cache for validated API keys
# Format: {provider:key_hash -> {validated: bool, timestamp: float, model_used: str}}
_api_key_validation_cache: Dict[str, Dict] = {}

# Thread lock for cache access
_cache_lock = threading.Lock()


def _hash_api_key(api_key: str) -> str:
    """
    Create a SHA-256 hash of an API key for cache lookup.

    Args:
        api_key: The API key to hash

    Returns:
        First 16 characters of the SHA-256 hash
    """
    if not api_key:
        return ""
    return hashlib.sha256(api_key.encode('utf-8')).hexdigest()[:16]


def _get_validation_model(provider: str) -> str:
    """
    Get the cheapest model for validation testing.

    Args:
        provider: Provider name ('openai', 'anthropic', 'openrouter')

    Returns:
        Model name string for validation

    Raises:
        ValueError: If provider is not supported
    """
    model = VALIDATION_MODELS.get(provider.lower())
    if not model:
        raise ValueError(
            f"Unsupported provider for validation: {provider}. "
            f"Supported providers: {', '.join(VALIDATION_MODELS.keys())}"
        )
    return model


def _get_cached_validation(provider: str, key_hash: str) -> Optional[Dict]:
    """
    Get cached validation result with thread safety.

    Args:
        provider: Provider name
        key_hash: Hash of the API key

    Returns:
        Cached validation dict or None if not found
    """
    with _cache_lock:
        cache_key = f"{provider}:{key_hash}"
        return _api_key_validation_cache.get(cache_key)


def _set_cached_validation(provider: str, key_hash: str, result: Dict) -> None:
    """
    Set cached validation result with thread safety.

    Args:
        provider: Provider name
        key_hash: Hash of the API key
        result: Validation result dict to cache
    """
    with _cache_lock:
        cache_key = f"{provider}:{key_hash}"
        _api_key_validation_cache[cache_key] = result


def _validate_single_provider(
    provider: str,
    api_key: Optional[str],
    force_revalidate: bool,
    verbose: bool
) -> Tuple[bool, Optional[str]]:
    """
    Validate a single provider's API key.

    Makes a minimal test API call to verify the key works. Uses intelligent
    caching to avoid redundant validations.

    Args:
        provider: Provider name ('openai', 'anthropic', 'openrouter')
        api_key: API key to validate (if None, reads from environment)
        force_revalidate: Skip cache and re-test even if previously validated
        verbose: Print validation status

    Returns:
        Tuple of (is_valid: bool, error_message: Optional[str])
    """
    # Get API key from environment if not provided
    if not api_key:
        env_var_names = {
            "openai": "OPENAI_API_KEY",
            "anthropic": "ANTHROPIC_API_KEY",
            "openrouter": "OPENROUTER_API_KEY",
        }
        env_var = env_var_names.get(provider.lower())
        if not env_var:
            return False, f"Unknown provider: {provider}"

        api_key = os.environ.get(env_var)
        if not api_key:
            error_msg = (
                f"API key not found. Environment variable {env_var} is not set. "
                f"Set it with: CASSIA.set_api_key('{provider}', 'your-key')"
            )
            return False, error_msg

    # Handle empty string
    api_key = api_key.strip()
    if not api_key:
        return False, "API key is empty"

    # Check cache (unless force_revalidate)
    key_hash = _hash_api_key(api_key)
    if not force_revalidate:
        cached = _get_cached_validation(provider, key_hash)
        if cached and cached.get("validated"):
            if verbose:
                logger.info(f"{provider.title()}: Already validated ✓")
            return True, None

    # Make test API call
    try:
        if verbose:
            logger.info(f"{provider.title()} ({_get_validation_model(provider)}): Testing...")

        model = _get_validation_model(provider)

        # gpt-5-nano only supports temperature=1.0
        test_temperature = 1.0 if "gpt-5-nano" in model.lower() else VALIDATION_TEMPERATURE

        # Make minimal test call
        response = call_llm(
            prompt=VALIDATION_PROMPT,
            provider=provider,
            model=model,
            api_key=api_key,
            temperature=test_temperature,
            max_tokens=VALIDATION_MAX_TOKENS
        )

        # Success - cache the result
        cache_result = {
            "validated": True,
            "timestamp": time.time(),
            "model_used": model
        }
        _set_cached_validation(provider, key_hash, cache_result)

        if verbose:
            logger.info(f"{provider.title()}: OK ✓")

        return True, None

    except Exception as e:
        error_str = str(e).lower()

        # Categorize errors with actionable messages
        if "401" in str(e) or "unauthorized" in error_str or "invalid api key" in error_str or "invalid_api_key" in error_str:
            error_msg = (
                f"Invalid API key. Please check your API key is valid. "
                f"Set it with: CASSIA.set_api_key('{provider}', 'your-key')"
            )
        elif "429" in str(e) or "rate limit" in error_str or "rate_limit" in error_str:
            error_msg = "Rate limit exceeded. Please wait and try again later."
        elif "timeout" in error_str or "timed out" in error_str:
            error_msg = "Request timed out. Please check your network connection and try again."
        elif "quota" in error_str or "insufficient" in error_str or "balance" in error_str:
            error_msg = "Insufficient credits or quota exceeded. Please check your account balance."
        elif "404" in str(e) or "not found" in error_str:
            error_msg = f"Model not found. The validation model '{_get_validation_model(provider)}' may not be available."
        else:
            # Generic error with details
            error_msg = f"Validation failed: {str(e)[:150]}"

        # Don't cache failures - always retry
        if verbose:
            logger.error(f"{provider.title()}: Failed - {error_msg}")

        return False, error_msg


def validate_api_keys(
    provider: Optional[str] = None,
    api_key: Optional[str] = None,
    force_revalidate: bool = False,
    verbose: bool = True
) -> Union[bool, Dict[str, bool]]:
    """
    Validate API keys by making a minimal test call to the provider(s).

    Uses intelligent caching to avoid redundant API calls. Each unique key is
    only tested once per session unless the key changes or force_revalidate=True.

    The validation makes a minimal test call using the cheapest model for each
    provider (costing ~$0.000001 per validation):
    - OpenAI: gpt-5-nano
    - Anthropic: claude-3-haiku-20240307
    - OpenRouter: openai/gpt-5-nano

    Args:
        provider: Specific provider to validate ('openai', 'anthropic',
                 'openrouter'). If None, validates all providers with keys set
                 in environment variables.
        api_key: Specific API key to validate. If None, reads from environment
                 variables. If provided, doesn't update environment variables.
        force_revalidate: Force revalidation even if key was previously
                          validated successfully. Default: False.
        verbose: Print validation status and results to logger. Default: True.

    Returns:
        - If provider specified: bool (True if valid, False if invalid)
        - If provider=None: dict mapping provider names to bool values

    Raises:
        ValueError: If provider is specified but not supported

    Examples:
        >>> import CASSIA
        >>>
        >>> # Validate all configured providers
        >>> CASSIA.validate_api_keys()
        Validating API keys...
        OpenAI (gpt-5-nano): Testing... OK ✓
        Anthropic (claude-3-haiku-20240307): Testing... OK ✓
        {'openai': True, 'anthropic': True}
        >>>
        >>> # Validate specific provider
        >>> CASSIA.validate_api_keys("openai")
        OpenAI: Already validated ✓
        True
        >>>
        >>> # Validate a key without setting it in environment
        >>> CASSIA.validate_api_keys("openai", api_key="sk-test...")
        OpenAI (gpt-5-nano): Testing... Failed
        False
        >>>
        >>> # Force revalidation (skip cache)
        >>> CASSIA.validate_api_keys("openai", force_revalidate=True)
        Validating API keys...
        OpenAI (gpt-5-nano): Testing... OK ✓
        True

    Notes:
        - First validation per key makes an API call (~$0.000001 cost)
        - Subsequent validations use cache (instant, no cost)
        - Cache automatically detects if key changes (via hash comparison)
        - Thread-safe for concurrent use in batch operations
        - Failed validations are NOT cached (always retry)
    """
    # If specific provider requested
    if provider is not None:
        provider = provider.lower()

        if verbose:
            logger.info("Validating API key...")

        is_valid, error_msg = _validate_single_provider(
            provider, api_key, force_revalidate, verbose
        )

        if not is_valid and error_msg and verbose:
            logger.error(f"Validation failed: {error_msg}")

        return is_valid

    # Otherwise, validate all providers with keys set
    if verbose:
        logger.info("Validating API keys...")

    env_var_names = {
        "openai": "OPENAI_API_KEY",
        "anthropic": "ANTHROPIC_API_KEY",
        "openrouter": "OPENROUTER_API_KEY",
    }

    results = {}
    providers_to_check = []

    # Find which providers have keys set
    for prov, env_var in env_var_names.items():
        if api_key or os.environ.get(env_var):
            providers_to_check.append(prov)

    if not providers_to_check:
        if verbose:
            logger.warning(
                "No API keys configured. Set API keys with: "
                "CASSIA.set_api_key(provider, 'your-key')"
            )
        return {}

    # Validate each provider
    for prov in providers_to_check:
        is_valid, error_msg = _validate_single_provider(
            prov, api_key, force_revalidate, verbose
        )
        results[prov] = is_valid

        if not is_valid and error_msg and verbose:
            logger.error(f"{prov.title()} validation failed: {error_msg}")

    return results


def clear_validation_cache(provider: Optional[str] = None) -> None:
    """
    Clear the API key validation cache.

    This can be useful if:
    - You've rotated your API keys and want to force revalidation
    - You suspect the cache is stale
    - You're debugging validation issues

    Args:
        provider: Specific provider to clear from cache. If None, clears
                 entire cache for all providers.

    Examples:
        >>> import CASSIA
        >>>
        >>> # Clear cache for specific provider
        >>> CASSIA.clear_validation_cache("openai")
        >>>
        >>> # Clear entire cache
        >>> CASSIA.clear_validation_cache()

    Notes:
        - This doesn't affect environment variables, only the validation cache
        - Next validation call will make an actual API test call
    """
    with _cache_lock:
        if provider is None:
            _api_key_validation_cache.clear()
            logger.info("Cleared entire validation cache")
        else:
            provider = provider.lower()
            # Remove all cache entries for this provider
            keys_to_remove = [k for k in _api_key_validation_cache if k.startswith(f"{provider}:")]
            for key in keys_to_remove:
                del _api_key_validation_cache[key]
            logger.info(f"Cleared validation cache for {provider}")
