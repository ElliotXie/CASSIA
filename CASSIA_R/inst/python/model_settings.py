"""
Simplified model settings for CASSIA.

Usage:
    from model_settings import resolve_model_name

    # Use tier shortcuts
    model, provider = resolve_model_name("best", "openai")      # -> ("gpt-5.1", "openai")
    model, provider = resolve_model_name("balanced", "anthropic") # -> ("claude-sonnet-4-5", "anthropic")
    model, provider = resolve_model_name("fast", "openrouter")    # -> ("google/gemini-2.5-flash", "openrouter")

    # Or use exact model names
    model, provider = resolve_model_name("gpt-4o", "openai")      # -> ("gpt-4o", "openai")
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
            }
        }

    def resolve_model_name(self, model_name: str, provider: str) -> Tuple[str, str]:
        """
        Resolve model name to actual model string.

        Args:
            model_name: Either a tier ("best", "balanced", "fast", "recommended")
                       or an exact model name
            provider: Provider name ("openai", "anthropic", "openrouter")

        Returns:
            Tuple of (resolved_model_name, provider)

        Examples:
            >>> resolve_model_name("best", "openai")
            ("gpt-5.1", "openai")
            >>> resolve_model_name("gpt-4o", "openai")
            ("gpt-4o", "openai")
        """
        if not model_name:
            raise ValueError("Model name cannot be empty")

        if not provider:
            raise ValueError("Provider must be specified (openai, anthropic, or openrouter)")

        provider = provider.lower()
        model_name_lower = model_name.lower()

        if provider not in VALID_PROVIDERS:
            raise ValueError(f"Unknown provider: {provider}. Must be one of: {VALID_PROVIDERS}")

        # Check if it's a tier shortcut
        if model_name_lower in VALID_TIERS:
            provider_settings = self.settings.get("providers", {}).get(provider, {})
            resolved = provider_settings.get(model_name_lower)
            if resolved:
                return resolved, provider
            else:
                raise ValueError(f"Tier '{model_name}' not found for provider '{provider}'")

        # Otherwise, return the model name as-is
        return model_name, provider

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


# Global instance
_model_settings = None


def get_model_settings() -> ModelSettings:
    """Get the global ModelSettings instance."""
    global _model_settings
    if _model_settings is None:
        _model_settings = ModelSettings()
    return _model_settings


def resolve_model_name(model_name: str, provider: str) -> Tuple[str, str]:
    """
    Resolve model name to actual model string.

    Args:
        model_name: Either a tier ("best", "balanced", "fast", "recommended")
                   or an exact model name
        provider: Provider name ("openai", "anthropic", "openrouter")

    Returns:
        Tuple of (resolved_model_name, provider)

    Examples:
        >>> resolve_model_name("best", "openai")
        ("gpt-5.1", "openai")
        >>> resolve_model_name("balanced", "anthropic")
        ("claude-sonnet-4-5", "anthropic")
        >>> resolve_model_name("fast", "openrouter")
        ("google/gemini-2.5-flash", "openrouter")
        >>> resolve_model_name("gpt-4o", "openai")  # exact name
        ("gpt-4o", "openai")
    """
    return get_model_settings().resolve_model_name(model_name, provider)


def get_recommended_model(provider: str) -> Tuple[str, str]:
    """
    Get the recommended model for a provider.

    Args:
        provider: Provider name ("openai", "anthropic", "openrouter")

    Returns:
        Tuple of (model_name, provider)
    """
    return get_model_settings().resolve_model_name("recommended", provider)


def print_available_models():
    """Print all available models."""
    get_model_settings().print_available_models()
