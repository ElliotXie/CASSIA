"""
Simplified model settings for CASSIA with fuzzy model name matching.

Usage:
    from model_settings import resolve_model_name

    # Use tier shortcuts
    model, provider = resolve_model_name("best", "openai")      # -> ("gpt-5.1", "openai")
    model, provider = resolve_model_name("balanced", "anthropic") # -> ("claude-sonnet-4-5", "anthropic")
    model, provider = resolve_model_name("fast", "openrouter")    # -> ("google/gemini-2.5-flash", "openrouter")

    # Use fuzzy aliases (prints: "Note: Resolved 'gpt' to 'gpt-5.1' for openai")
    model, provider = resolve_model_name("gpt", "openai")         # -> ("gpt-5.1", "openai")
    model, provider = resolve_model_name("claude", "anthropic")   # -> ("claude-sonnet-4-5", "anthropic")
    model, provider = resolve_model_name("gemini", "openrouter")  # -> ("google/gemini-2.5-flash", "openrouter")

    # Or use exact model names (no resolution note printed)
    model, provider = resolve_model_name("gpt-4o", "openai")      # -> ("gpt-4o", "openai")

Resolution Priority:
    1. Tier shortcuts (best, balanced, fast, recommended)
    2. Provider-specific aliases (gpt, claude, gemini, etc.)
    3. Global aliases (sonnet, opus, haiku)
    4. Pass-through (exact model name)
"""

import json
from typing import Dict, Tuple, Optional
from pathlib import Path


# Valid tier shortcuts
VALID_TIERS = {"best", "balanced", "fast", "recommended"}

# Valid providers
VALID_PROVIDERS = {"openai", "anthropic", "openrouter"}


class ModelSettings:
    """Simple model settings manager."""

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize ModelSettings.

        Args:
            config_path: Path to model_settings.json. If None, uses default location.
        """
        if config_path is None:
            current_dir = Path(__file__).parent
            possible_paths = [
                current_dir / "data" / "model_settings.json",
                current_dir / "model_settings.json",
            ]
            for path in possible_paths:
                if path.exists():
                    config_path = path
                    break
            if config_path is None:
                config_path = current_dir / "data" / "model_settings.json"

        self.config_path = Path(config_path)
        self.settings = self._load_settings()

    def _load_settings(self) -> Dict:
        """Load settings from JSON file."""
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return self._get_fallback_settings()

    def _get_fallback_settings(self) -> Dict:
        """Fallback settings if JSON file is not available."""
        return {
            "providers": {
                "openai": {
                    "best": "gpt-5.1",
                    "balanced": "gpt-4o",
                    "fast": "gpt-5-mini",
                    "recommended": "gpt-5.1"
                },
                "anthropic": {
                    "best": "claude-opus-4-5",
                    "balanced": "claude-sonnet-4-5",
                    "fast": "claude-haiku-4-5",
                    "recommended": "claude-sonnet-4-5"
                },
                "openrouter": {
                    "best": "anthropic/claude-sonnet-4.5",
                    "balanced": "openai/gpt-5.1",
                    "fast": "google/gemini-2.5-flash",
                    "recommended": "anthropic/claude-sonnet-4.5"
                }
            },
            "aliases": self._get_fallback_aliases()
        }

    def _get_fallback_aliases(self) -> Dict:
        """Fallback aliases if not in JSON config."""
        return {
            "provider_specific": {
                "openai": {
                    "gpt": "gpt-5.1",
                    "gpt4": "gpt-4o",
                    "gpt-4": "gpt-4o",
                    "4o": "gpt-4o",
                    "gpt4o": "gpt-4o",
                    "gpt5": "gpt-5.1",
                    "gpt-5": "gpt-5.1",
                    "mini": "gpt-5-mini",
                    "gpt-mini": "gpt-5-mini"
                },
                "anthropic": {
                    "claude": "claude-sonnet-4-5",
                    "sonnet": "claude-sonnet-4-5",
                    "opus": "claude-opus-4-5",
                    "haiku": "claude-haiku-4-5"
                },
                "openrouter": {
                    "gpt": "openai/gpt-5.1",
                    "claude": "anthropic/claude-sonnet-4.5",
                    "sonnet": "anthropic/claude-sonnet-4.5",
                    "opus": "anthropic/claude-opus-4.5",
                    "haiku": "anthropic/claude-haiku-4.5",
                    "gemini": "google/gemini-2.5-flash",
                    "flash": "google/gemini-2.5-flash",
                    "deepseek": "deepseek/deepseek-chat"
                }
            },
            "global": {
                "sonnet": "claude-sonnet-4-5",
                "opus": "claude-opus-4-5",
                "haiku": "claude-haiku-4-5"
            }
        }

    def _resolve_alias(self, model_name: str, provider: str) -> Tuple[Optional[str], bool]:
        """
        Resolve an alias to a model name.

        Args:
            model_name: The alias to resolve
            provider: The provider name

        Returns:
            Tuple of (resolved_name, was_resolved)
            - resolved_name: The resolved model name, or None if not found
            - was_resolved: True if an alias was matched
        """
        model_lower = model_name.lower().strip()
        aliases = self.settings.get("aliases", self._get_fallback_aliases())

        # Try provider-specific alias first
        provider_aliases = aliases.get("provider_specific", {}).get(provider, {})
        if model_lower in provider_aliases:
            return provider_aliases[model_lower], True

        # Try global aliases
        global_aliases = aliases.get("global", {})
        if model_lower in global_aliases:
            return global_aliases[model_lower], True

        return None, False

    def resolve_model_name(self, model_name: str, provider: str, verbose: bool = True) -> Tuple[str, str]:
        """
        Resolve model name to actual model string.

        Resolution priority:
            1. Tier shortcuts (best, balanced, fast, recommended)
            2. Provider-specific aliases (gpt, claude, gemini, etc.)
            3. Global aliases (sonnet, opus, haiku)
            4. Pass-through (exact model name)

        Args:
            model_name: Model name, tier shortcut, or alias
            provider: Provider name ("openai", "anthropic", "openrouter")
            verbose: Print resolution messages when alias is used (default: True)

        Returns:
            Tuple of (resolved_model_name, provider)

        Examples:
            >>> resolve_model_name("best", "openai")
            ("gpt-5.1", "openai")
            >>> resolve_model_name("gpt", "openai")  # prints: Note: Resolved 'gpt' to 'gpt-5.1' for openai
            ("gpt-5.1", "openai")
            >>> resolve_model_name("gpt-4o", "openai")  # exact name, no note
            ("gpt-4o", "openai")
        """
        if not model_name:
            raise ValueError("Model name cannot be empty")

        if not provider:
            raise ValueError("Provider must be specified (openai, anthropic, openrouter, or an HTTP base URL)")

        provider_clean = provider.strip()
        if not provider_clean:
            raise ValueError("Provider must be specified (openai, anthropic, openrouter, or an HTTP base URL)")

        provider_lower = provider_clean.lower()
        model_name_lower = model_name.lower().strip()

        is_custom_provider = provider_lower.startswith("http://") or provider_lower.startswith("https://")

        if not is_custom_provider and provider_lower not in VALID_PROVIDERS:
            raise ValueError(f"Unknown provider: {provider_clean}. Must be one of: {VALID_PROVIDERS} or an HTTP/HTTPS base URL")

        # For standard providers, support tier shortcuts and aliases
        if not is_custom_provider:
            if model_name_lower in VALID_TIERS:
                provider_settings = self.settings.get("providers", {}).get(provider_lower, {})
                resolved = provider_settings.get(model_name_lower)
                if resolved:
                    return resolved, provider_lower
                else:
                    raise ValueError(f"Tier '{model_name}' not found for provider '{provider_clean}'")

            resolved, was_alias = self._resolve_alias(model_name, provider_lower)
            if was_alias and resolved:
                if verbose:
                    print(f"Note: Resolved '{model_name}' to '{resolved}' for {provider_lower}")
                return resolved, provider_lower

            return model_name, provider_lower

        # Custom HTTP providers: pass through model name to preserve user-specified values
        resolved, was_alias = self._resolve_alias(model_name, provider_lower)
        if was_alias and resolved:
            if verbose:
                print(f"Note: Resolved '{model_name}' to '{resolved}' for {provider_clean}")
            return resolved, provider_clean

        return model_name, provider_clean

    def get_available_tiers(self) -> list:
        """Get list of available tier shortcuts."""
        return list(VALID_TIERS)

    def get_available_providers(self) -> list:
        """Get list of available providers."""
        return list(VALID_PROVIDERS)

    def get_model_for_tier(self, tier: str, provider: str) -> str:
        """
        Get the model name for a specific tier and provider.

        Args:
            tier: One of "best", "balanced", "fast", "recommended"
            provider: One of "openai", "anthropic", "openrouter"

        Returns:
            Model name string
        """
        model, _ = self.resolve_model_name(tier, provider)
        return model

    def get_available_aliases(self, provider: Optional[str] = None) -> Dict:
        """
        Get available aliases for fuzzy model matching.

        Args:
            provider: Optional provider to filter aliases

        Returns:
            Dictionary of aliases
        """
        aliases = self.settings.get("aliases", self._get_fallback_aliases())

        if provider:
            provider = provider.lower()
            return {
                "provider_specific": aliases.get("provider_specific", {}).get(provider, {}),
                "global": aliases.get("global", {})
            }
        return aliases

    def print_available_models(self):
        """Print all available models in a readable format."""
        print("\n=== Available Models ===\n")
        print("Tiers: best, balanced, fast, recommended\n")

        providers = self.settings.get("providers", {})
        for provider, tiers in providers.items():
            print(f"{provider.upper()}:")
            for tier, model in tiers.items():
                print(f"  {tier:12} -> {model}")
            print()

    def print_available_aliases(self):
        """Print all available aliases in a readable format."""
        print("\n=== Available Model Aliases ===\n")
        aliases = self.settings.get("aliases", self._get_fallback_aliases())

        print("Provider-Specific Aliases:")
        for provider, provider_aliases in aliases.get("provider_specific", {}).items():
            print(f"\n  {provider.upper()}:")
            for alias, model in sorted(provider_aliases.items()):
                print(f"    {alias:15} -> {model}")

        print("\n\nGlobal Aliases (work with any provider):")
        for alias, model in sorted(aliases.get("global", {}).items()):
            print(f"    {alias:15} -> {model}")
        print()


# Global instance
_model_settings = None


def get_model_settings() -> ModelSettings:
    """Get the global ModelSettings instance."""
    global _model_settings
    if _model_settings is None:
        _model_settings = ModelSettings()
    return _model_settings


def resolve_model_name(model_name: str, provider: str, verbose: bool = True) -> Tuple[str, str]:
    """
    Resolve model name to actual model string.

    Supports:
    - Tier shortcuts: "best", "balanced", "fast", "recommended"
    - Aliases: "gpt", "claude", "sonnet", "opus", "haiku", "gemini", "flash", "deepseek"
    - Exact model names: passed through unchanged

    Args:
        model_name: Model name, tier shortcut, or alias
        provider: Provider name ("openai", "anthropic", "openrouter")
        verbose: Print resolution messages when alias is used (default: True)

    Returns:
        Tuple of (resolved_model_name, provider)

    Examples:
        >>> resolve_model_name("best", "openai")
        ('gpt-5.1', 'openai')
        >>> resolve_model_name("gpt", "openai")  # prints: Note: Resolved 'gpt' to 'gpt-5.1' for openai
        ('gpt-5.1', 'openai')
        >>> resolve_model_name("claude", "anthropic")  # prints: Note: Resolved 'claude' to 'claude-sonnet-4-5' for anthropic
        ('claude-sonnet-4-5', 'anthropic')
        >>> resolve_model_name("gemini", "openrouter")  # prints: Note: Resolved 'gemini' to 'google/gemini-2.5-flash' for openrouter
        ('google/gemini-2.5-flash', 'openrouter')
        >>> resolve_model_name("gpt-4o", "openai")  # no message, exact match
        ('gpt-4o', 'openai')
    """
    return get_model_settings().resolve_model_name(model_name, provider, verbose)


def get_recommended_model(provider: str) -> Tuple[str, str]:
    """
    Get the recommended model for a provider.

    Args:
        provider: Provider name ("openai", "anthropic", "openrouter")

    Returns:
        Tuple of (model_name, provider)
    """
    return get_model_settings().resolve_model_name("recommended", provider, verbose=False)


def get_available_aliases(provider: Optional[str] = None) -> Dict:
    """
    Get available aliases for fuzzy model matching.

    Args:
        provider: Optional provider to filter aliases

    Returns:
        Dictionary of aliases
    """
    return get_model_settings().get_available_aliases(provider)


def print_available_models():
    """Print all available models."""
    get_model_settings().print_available_models()


def print_available_aliases():
    """Print all available aliases."""
    get_model_settings().print_available_aliases()
