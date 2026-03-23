"""
Simplified model settings for CASSIA with fuzzy model name matching.

All model configuration is stored in data/model_settings.json (the single source of truth).
No hardcoded model names exist in this file — everything is read from the JSON config.

Usage:
    from model_settings import resolve_model_name

    # Use tier shortcuts
    model, provider = resolve_model_name("best", "openai")
    model, provider = resolve_model_name("balanced", "anthropic")
    model, provider = resolve_model_name("fast", "openrouter")

    # Use fuzzy aliases (prints: "Note: Resolved 'gpt' to '...' for openai")
    model, provider = resolve_model_name("gpt", "openai")
    model, provider = resolve_model_name("claude", "anthropic")
    model, provider = resolve_model_name("gemini", "openrouter")

    # Or use exact model names (no resolution note printed)
    model, provider = resolve_model_name("gpt-4o", "openai")

Resolution Priority:
    1. Tier shortcuts (best, balanced, fast, recommended)
    2. Provider-specific aliases (gpt, claude, gemini, etc.)
    3. Global aliases (sonnet, opus, haiku)
    4. Pass-through (exact model name)
"""

import json
from typing import Dict, List, Tuple, Optional
from pathlib import Path


# Valid tier shortcuts
VALID_TIERS = {"best", "balanced", "fast", "recommended"}

# Valid providers
VALID_PROVIDERS = {"openai", "anthropic", "openrouter"}


class ModelSettings:
    """Simple model settings manager. Reads all config from model_settings.json."""

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize ModelSettings.

        Args:
            config_path: Path to model_settings.json. If None, uses default location.
        """
        if config_path is None:
            current_dir = Path(__file__).parent
            package_dir = current_dir.parent  # CASSIA/ package root
            possible_paths = [
                package_dir / "data" / "model_settings.json",
                current_dir / "data" / "model_settings.json",
                current_dir / "model_settings.json",
            ]
            for path in possible_paths:
                if path.exists():
                    config_path = path
                    break
            if config_path is None:
                config_path = package_dir / "data" / "model_settings.json"

        self.config_path = Path(config_path)
        self.settings = self._load_settings()

    def _load_settings(self) -> Dict:
        """Load settings from JSON file. Raises error if file is missing or invalid."""
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Model settings file not found: {self.config_path}\n"
                f"This file is required for CASSIA to work. "
                f"Please reinstall the package or restore data/model_settings.json."
            )
        except json.JSONDecodeError as e:
            raise ValueError(
                f"Invalid JSON in model settings file: {self.config_path}\n"
                f"Error: {e}"
            )

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
        aliases = self.settings.get("aliases", {})

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
        aliases = self.settings.get("aliases", {})

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
        aliases = self.settings.get("aliases", {})

        print("Provider-Specific Aliases:")
        for provider, provider_aliases in aliases.get("provider_specific", {}).items():
            print(f"\n  {provider.upper()}:")
            for alias, model in sorted(provider_aliases.items()):
                print(f"    {alias:15} -> {model}")

        print("\n\nGlobal Aliases (work with any provider):")
        for alias, model in sorted(aliases.get("global", {}).items()):
            print(f"    {alias:15} -> {model}")
        print()

    def get_pipeline_defaults(self, provider: str) -> Dict[str, str]:
        """
        Get default models for each pipeline stage for a given provider.

        Args:
            provider: Provider name ("openai", "anthropic", "openrouter")

        Returns:
            Dict mapping stage names to model names
        """
        provider_lower = provider.lower().strip()
        defaults = self.settings.get("pipeline_defaults", {})
        provider_defaults = defaults.get(provider_lower, {})

        if not provider_defaults:
            # Build from provider tiers as fallback
            providers = self.settings.get("providers", {}).get(provider_lower, {})
            return {
                "annotation": providers.get("recommended", providers.get("best", "")),
                "score": providers.get("recommended", providers.get("best", "")),
                "merge": providers.get("fast", ""),
                "annotationboost": providers.get("recommended", providers.get("best", ""))
            }

        return provider_defaults

    def get_agent_default(self, agent_name: str, provider: str) -> Dict[str, any]:
        """
        Get default model and temperature for a specific agent.

        Args:
            agent_name: One of 'annotation', 'scoring', 'merging', 'subclustering',
                       'subclustering_n', 'annotation_boost', 'uncertainty'
            provider: Provider name ("openai", "anthropic", "openrouter")

        Returns:
            Dict with 'model' and 'temperature' keys
        """
        provider_lower = provider.lower().strip()
        agent_lower = agent_name.lower().strip()

        # Try agent_defaults from JSON
        agent_defaults = self.settings.get("agent_defaults", {})
        provider_agents = agent_defaults.get(provider_lower, {})

        if agent_lower in provider_agents:
            return provider_agents[agent_lower]

        # Build from provider tiers as fallback
        providers = self.settings.get("providers", {}).get(provider_lower, {})
        default_model = providers.get("recommended", providers.get("best", ""))
        return {"model": default_model, "temperature": 0}

    def get_symphony_presets(self) -> Dict[str, List[str]]:
        """
        Get symphony compare model presets.

        Returns:
            Dict with 'premium' and 'budget' keys, each containing a list of model IDs
        """
        return self.settings.get("symphony_presets", {"premium": [], "budget": []})

    def get_model_personas(self) -> Dict[str, str]:
        """
        Get model persona mappings (model ID -> scientist name).

        Returns:
            Dict mapping model IDs to persona names
        """
        return self.settings.get("model_personas", {})

    def get_validation_models(self) -> Dict[str, str]:
        """
        Get validation models for API key testing (cheapest models per provider).

        Returns:
            Dict mapping provider names to validation model IDs
        """
        return self.settings.get("validation_models", {})


# Global instance with thread-safe initialization
_model_settings = None
_model_settings_lock = __import__('threading').Lock()


def get_model_settings() -> ModelSettings:
    """Get the global ModelSettings instance (thread-safe)."""
    global _model_settings
    if _model_settings is None:
        with _model_settings_lock:
            if _model_settings is None:  # Double-check after acquiring lock
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


def get_pipeline_defaults(provider: str) -> Dict[str, str]:
    """
    Get default models for each pipeline stage for a given provider.

    Args:
        provider: Provider name ("openai", "anthropic", "openrouter")

    Returns:
        Dict mapping stage names to model names
    """
    return get_model_settings().get_pipeline_defaults(provider)


def get_agent_default(agent_name: str, provider: str) -> Dict[str, any]:
    """
    Get default model and temperature for a specific agent.

    Args:
        agent_name: One of 'annotation', 'scoring', 'merging', 'subclustering',
                   'subclustering_n', 'annotation_boost', 'uncertainty'
        provider: Provider name ("openai", "anthropic", "openrouter")

    Returns:
        Dict with 'model' and 'temperature' keys
    """
    return get_model_settings().get_agent_default(agent_name, provider)
