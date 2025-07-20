import json
import os
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path

class ModelSettings:
    """
    A class to manage model settings and provide unified model name resolution.
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the ModelSettings class.
        
        Args:
            config_path: Path to the model settings JSON file. If None, uses default location.
        """
        if config_path is None:
            # Try to find the model_settings.json file in multiple locations
            current_dir = Path(__file__).parent
            
            # Try locations in order of preference
            possible_paths = [
                current_dir / "data" / "model_settings.json",  # Package data directory (preferred)
                current_dir / "model_settings.json",  # Same directory as model_settings.py
                current_dir.parent / "model_settings.json",  # Parent directory
                current_dir.parent.parent / "model_settings.json",  # Project root
                current_dir.parent.parent.parent / "model_settings.json",  # Project root parent
            ]
            
            config_path = None
            for path in possible_paths:
                if path.exists():
                    config_path = path
                    break
            
            if config_path is None:
                # If not found anywhere, default to package data directory
                config_path = current_dir / "data" / "model_settings.json"
        
        self.config_path = Path(config_path)
        self.settings = self._load_settings()
    
    def _load_settings(self) -> Dict[str, Any]:
        """Load model settings from JSON file."""
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: Model settings file not found at {self.config_path}")
            print("Using built-in fallback configuration.")
            return self._get_fallback_settings()
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in model settings file: {e}")
            print("Using built-in fallback configuration.")
            return self._get_fallback_settings()
    
    def _get_fallback_settings(self) -> Dict[str, Any]:
        """Get fallback model settings when JSON file is not available."""
        return {
            "providers": {
                "openai": {
                    "name": "OpenAI",
                    "default_model": "gpt-4o",
                    "recommended_model": "gpt-4o",
                    "models": {
                        "gpt-4o": {
                            "description": "OpenAI's most capable model",
                            "cost_tier": "high",
                            "context_window": 128000,
                            "use_cases": ["annotation", "scoring", "annotation_boost", "merging"]
                        },
                        "gpt-4o-mini": {
                            "description": "OpenAI's fast and efficient model",
                            "cost_tier": "medium",
                            "context_window": 128000,
                            "use_cases": ["annotation", "scoring", "merging"]
                        },
                        "gpt-3.5-turbo": {
                            "description": "OpenAI's fast, cost-effective model",
                            "cost_tier": "low",
                            "context_window": 16000,
                            "use_cases": ["annotation", "scoring"]
                        }
                    }
                },
                "anthropic": {
                    "name": "Anthropic",
                    "default_model": "claude-3-5-sonnet-latest",
                    "recommended_model": "claude-3-5-sonnet-latest",
                    "models": {
                        "claude-3-5-sonnet-latest": {
                            "description": "Anthropic's most capable model",
                            "cost_tier": "high",
                            "context_window": 200000,
                            "use_cases": ["annotation", "scoring", "annotation_boost", "merging"]
                        },
                        "claude-3-5-haiku-latest": {
                            "description": "Anthropic's fast and efficient model",
                            "cost_tier": "low",
                            "context_window": 200000,
                            "use_cases": ["annotation", "scoring", "merging"]
                        }
                    }
                },
                "openrouter": {
                    "name": "OpenRouter",
                    "default_model": "google/gemini-2.5-flash",
                    "recommended_model": "google/gemini-2.5-flash",
                    "models": {
                        "google/gemini-2.5-flash": {
                            "description": "Google's fast and efficient model",
                            "cost_tier": "low",
                            "context_window": 1000000,
                            "use_cases": ["annotation", "scoring", "annotation_boost", "merging"]
                        },
                        "deepseek/deepseek-chat-v3-0324": {
                            "description": "DeepSeek's very cost-effective model",
                            "cost_tier": "very_low",
                            "context_window": 64000,
                            "use_cases": ["annotation", "scoring"]
                        },
                        "anthropic/claude-3-5-sonnet": {
                            "description": "Anthropic Claude via OpenRouter",
                            "cost_tier": "high",
                            "context_window": 200000,
                            "use_cases": ["annotation", "scoring", "annotation_boost", "merging"]
                        },
                        "openai/gpt-4o": {
                            "description": "OpenAI GPT-4o via OpenRouter",
                            "cost_tier": "high",
                            "context_window": 128000,
                            "use_cases": ["annotation", "scoring", "annotation_boost", "merging"]
                        }
                    }
                }
            },
            "provider_shortcuts": {
                "openai": {
                    "gpt4": "gpt-4o",
                    "gpt": "gpt-4o",
                    "best": "gpt-4o",
                    "premium": "gpt-4o",
                    "cheap": "gpt-3.5-turbo",
                    "fast": "gpt-4o-mini"
                },
                "anthropic": {
                    "claude": "claude-3-5-sonnet-latest",
                    "sonnet": "claude-3-5-sonnet-latest",
                    "best": "claude-3-5-sonnet-latest",
                    "premium": "claude-3-5-sonnet-latest",
                    "cheap": "claude-3-5-haiku-latest",
                    "fast": "claude-3-5-haiku-latest"
                },
                "openrouter": {
                    "gpt4": "openai/gpt-4o",
                    "gpt": "openai/gpt-4o",
                    "claude": "anthropic/claude-3-5-sonnet",
                    "sonnet": "anthropic/claude-3-5-sonnet",
                    "gemini": "google/gemini-2.5-flash",
                    "flash": "google/gemini-2.5-flash",
                    "deepseek": "deepseek/deepseek-chat-v3-0324",
                    "best": "google/gemini-2.5-flash",
                    "fast": "google/gemini-2.5-flash",
                    "premium": "anthropic/claude-3-5-sonnet",
                    "cheap": "deepseek/deepseek-chat-v3-0324"
                }
            }
        }
    
    def resolve_model_name(self, model_name: str, provider: str = None) -> Tuple[str, str]:
        """
        Resolve a user-provided model name to the actual provider model name.
        
        Args:
            model_name: User-provided model name (can be alias or actual name)
            provider: Provider name (REQUIRED - users must specify which API they want to use)
        
        Returns:
            Tuple of (actual_model_name, provider_name)
        """
        if not model_name:
            raise ValueError("Model name cannot be empty")
        
        if not provider:
            raise ValueError("Provider must be specified (openai, anthropic, or openrouter)")
        
        provider = provider.lower()
        
        # Check provider shortcuts first
        provider_shortcuts = self.settings.get("provider_shortcuts", {}).get(provider, {})
        if model_name in provider_shortcuts:
            resolved_name = provider_shortcuts[model_name]
            return resolved_name, provider
        
        # If provider is specified, search within that provider
        provider_info = self.settings.get("providers", {}).get(provider, {})
        if not provider_info:
            raise ValueError(f"Unknown provider: {provider}")
        
        models = provider_info.get("models", {})
        
        # Check if it's an actual model name
        if model_name in models:
            return model_name, provider
        
        # Check if it's an alias
        for actual_name, model_info in models.items():
            aliases = model_info.get("aliases", [])
            if model_name in aliases:
                return actual_name, provider
        
        # If still not found, check if it looks like a provider/model format
        if "/" in model_name:
            return model_name, provider
        
        # If nothing found, return as-is with the specified provider
        return model_name, provider
    
    def _infer_provider_from_model(self, model_name: str) -> Optional[str]:
        """Infer provider from model name patterns."""
        if model_name.startswith("gpt-") or model_name.startswith("openai/"):
            return "openai"
        elif model_name.startswith("claude-") or model_name.startswith("anthropic/"):
            return "anthropic"
        elif "/" in model_name:
            return "openrouter"
        else:
            return None
    
    def get_recommended_model(self, provider: str = None, use_case: str = None) -> Tuple[str, str]:
        """
        Get the recommended model for a provider or use case.
        
        Args:
            provider: Provider name (REQUIRED - users must specify which API they want to use)
            use_case: Use case (annotation, scoring, annotation_boost, merging)
        
        Returns:
            Tuple of (model_name, provider_name)
        """
        if not provider:
            raise ValueError("Provider must be specified (openai, anthropic, or openrouter)")
        
        provider = provider.lower()
        provider_info = self.settings.get("providers", {}).get(provider, {})
        
        if not provider_info:
            raise ValueError(f"Unknown provider: {provider}")
        
        # If use_case is specified, get the provider-specific best model for that use case
        if use_case:
            # Check if the provider has a specific shortcut for this use case
            shortcuts = self.settings.get("provider_shortcuts", {}).get(provider, {})
            if "best" in shortcuts:
                return shortcuts["best"], provider
        
        # Return provider's recommended model
        recommended_model = provider_info.get("recommended_model")
        if recommended_model:
            return recommended_model, provider
        
        # Return provider's default model
        default_model = provider_info.get("default_model")
        if default_model:
            return default_model, provider
        
        # Fallback to first model in the provider
        models = provider_info.get("models", {})
        if models:
            first_model = list(models.keys())[0]
            return first_model, provider
        
        raise ValueError(f"No models found for provider: {provider}")
    
    def get_model_info(self, model_name: str, provider: str = None) -> Optional[Dict[str, Any]]:
        """
        Get detailed information about a model.
        
        Args:
            model_name: Model name
            provider: Provider name (if None, will try to infer)
        
        Returns:
            Dictionary with model information or None if not found
        """
        actual_name, provider_name = self.resolve_model_name(model_name, provider)
        
        provider_info = self.settings.get("providers", {}).get(provider_name, {})
        models = provider_info.get("models", {})
        
        return models.get(actual_name)
    
    def list_models(self, provider: str = None, cost_tier: str = None, use_case: str = None) -> List[Dict[str, Any]]:
        """
        List available models with optional filtering.
        
        Args:
            provider: Filter by provider
            cost_tier: Filter by cost tier (very_low, low, medium, high)
            use_case: Filter by use case
        
        Returns:
            List of model dictionaries with name, provider, and info
        """
        results = []
        
        providers_to_check = [provider] if provider else self.settings.get("providers", {}).keys()
        
        for provider_name in providers_to_check:
            provider_info = self.settings.get("providers", {}).get(provider_name, {})
            models = provider_info.get("models", {})
            
            for model_name, model_info in models.items():
                # Apply filters
                if cost_tier and model_info.get("cost_tier") != cost_tier:
                    continue
                
                if use_case and use_case not in model_info.get("use_cases", []):
                    continue
                
                results.append({
                    "name": model_name,
                    "provider": provider_name,
                    "info": model_info
                })
        
        return results
    
    def get_model_aliases(self, model_name: str, provider: str = None) -> List[str]:
        """
        Get all aliases for a model.
        
        Args:
            model_name: Model name
            provider: Provider name (if None, will try to infer)
        
        Returns:
            List of aliases
        """
        model_info = self.get_model_info(model_name, provider)
        if model_info:
            return model_info.get("aliases", [])
        return []
    
    def get_use_case_recommendations(self, use_case: str) -> Dict[str, Any]:
        """
        Get recommendations for a specific use case.
        
        Args:
            use_case: Use case name
        
        Returns:
            Dictionary with recommendations
        """
        return self.settings.get("use_case_recommendations", {}).get(use_case, {})
    
    def get_cost_tier_models(self, cost_tier: str) -> List[str]:
        """
        Get models in a specific cost tier.
        
        Args:
            cost_tier: Cost tier name
        
        Returns:
            List of model names
        """
        cost_info = self.settings.get("cost_tiers", {}).get(cost_tier, {})
        return cost_info.get("models", [])

# Global instance
_model_settings = None

def get_model_settings() -> ModelSettings:
    """Get the global ModelSettings instance."""
    global _model_settings
    if _model_settings is None:
        _model_settings = ModelSettings()
    return _model_settings

def resolve_model_name(model_name: str, provider: str = None) -> Tuple[str, str]:
    """
    Resolve a user-provided model name to the actual provider model name.
    
    Args:
        model_name: User-provided model name (can be alias or actual name)
        provider: Provider name (REQUIRED - must specify openai, anthropic, or openrouter)
    
    Returns:
        Tuple of (actual_model_name, provider_name)
    """
    return get_model_settings().resolve_model_name(model_name, provider)

def get_recommended_model(provider: str = None, use_case: str = None) -> Tuple[str, str]:
    """
    Get the recommended model for a provider or use case.
    
    Args:
        provider: Provider name (REQUIRED - must specify openai, anthropic, or openrouter)
        use_case: Use case (annotation, scoring, annotation_boost, merging)
    
    Returns:
        Tuple of (model_name, provider_name)
    """
    return get_model_settings().get_recommended_model(provider, use_case)

def get_model_info(model_name: str, provider: str = None) -> Optional[Dict[str, Any]]:
    """
    Get detailed information about a model.
    
    Args:
        model_name: Model name
        provider: Provider name (if None, will try to infer)
    
    Returns:
        Dictionary with model information or None if not found
    """
    return get_model_settings().get_model_info(model_name, provider)

def list_models(provider: str = None, cost_tier: str = None, use_case: str = None) -> List[Dict[str, Any]]:
    """
    List available models with optional filtering.
    
    Args:
        provider: Filter by provider
        cost_tier: Filter by cost tier (very_low, low, medium, high)
        use_case: Filter by use case
    
    Returns:
        List of model dictionaries with name, provider, and info
    """
    return get_model_settings().list_models(provider, cost_tier, use_case)

def get_use_case_recommendations(use_case: str) -> Dict[str, Any]:
    """
    Get recommendations for a specific use case.
    
    Args:
        use_case: Use case name
    
    Returns:
        Dictionary with recommendations
    """
    return get_model_settings().get_use_case_recommendations(use_case)

def print_model_recommendations(use_case: str = None):
    """
    Print model recommendations in a user-friendly format.
    
    Args:
        use_case: Specific use case to show recommendations for
    """
    settings = get_model_settings()
    
    if use_case:
        recommendations = settings.get_use_case_recommendations(use_case)
        if recommendations:
            print(f"\n=== Recommendations for {use_case} ===")
            print(f"Best: {recommendations.get('best', 'N/A')}")
            print(f"Description: {recommendations.get('description', 'N/A')}")
            alternatives = recommendations.get('alternatives', [])
            if alternatives:
                print(f"Alternatives: {', '.join(alternatives)}")
        else:
            print(f"No recommendations found for use case: {use_case}")
    else:
        print("\n=== Model Recommendations by Use Case ===")
        use_cases = settings.settings.get("use_case_recommendations", {})
        for case, info in use_cases.items():
            print(f"\n{case.title()}:")
            print(f"  Best: {info.get('best', 'N/A')}")
            print(f"  Description: {info.get('description', 'N/A')}")
            alternatives = info.get('alternatives', [])
            if alternatives:
                print(f"  Alternatives: {', '.join(alternatives[:3])}")
    
    print("\n=== Cost Tiers ===")
    cost_tiers = settings.settings.get("cost_tiers", {})
    for tier, info in cost_tiers.items():
        print(f"{tier.title()}: {info.get('description', 'N/A')}")
        models = info.get('models', [])
        if models:
            print(f"  Models: {', '.join(models[:3])}")
    
    print("\n=== Unified Names ===")
    print("You can use these simple names instead of full model names:")
    unified = settings.settings.get("unified_names", {})
    for alias, actual in unified.items():
        print(f"  {alias} -> {actual}")