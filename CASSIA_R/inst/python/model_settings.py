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
            # Try to find the model_settings.json file in the data directory
            current_dir = Path(__file__).parent
            
            # First try data directory in the same location (for CASSIA_python structure)
            config_path = current_dir / "data" / "model_settings.json"
            
            # If not found, try the CASSIA_python/CASSIA/data structure
            if not config_path.exists():
                cassia_python_dir = current_dir.parent.parent.parent / "CASSIA_python" / "CASSIA" / "data"
                config_path = cassia_python_dir / "model_settings.json"
            
            # If still not found, try relative to current directory
            if not config_path.exists():
                project_root = current_dir.parent.parent
                config_path = project_root / "model_settings.json"
        
        self.config_path = Path(config_path)
        self.settings = self._load_settings()
    
    def _load_settings(self) -> Dict[str, Any]:
        """Load model settings from JSON file."""
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(f"Model settings file not found at {self.config_path}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in model settings file: {e}")
    
    def resolve_model_name(self, model_name: str, provider: str = None) -> Tuple[str, str]:
        """
        Resolve a user-provided model name to the actual provider model name.
        
        Args:
            model_name: User-provided model name (can be alias or actual name)
            provider: Provider name (if None, will try to infer from model_name)
        
        Returns:
            Tuple of (actual_model_name, provider_name)
        """
        if not model_name:
            raise ValueError("Model name cannot be empty")
        
        # First, check if it's a unified name
        unified_names = self.settings.get("unified_names", {})
        if model_name in unified_names:
            resolved_name = unified_names[model_name]
            # Determine provider from the resolved name
            provider_name = self._infer_provider_from_model(resolved_name)
            return resolved_name, provider_name
        
        # Check provider shortcuts for all providers
        provider_shortcuts = self.settings.get("provider_shortcuts", {})
        for provider_name, shortcuts in provider_shortcuts.items():
            if model_name in shortcuts:
                resolved_name = shortcuts[model_name]
                return resolved_name, provider_name
        
        # If provider is specified, search within that provider
        if provider:
            provider_info = self.settings.get("providers", {}).get(provider, {})
            models = provider_info.get("models", {})
            
            # Check if it's an actual model name
            if model_name in models:
                return model_name, provider
            
            # Check if it's an alias
            for actual_name, model_info in models.items():
                aliases = model_info.get("aliases", [])
                if model_name in aliases:
                    return actual_name, provider
        
        # If no provider specified or not found, search all providers
        for provider_name, provider_info in self.settings.get("providers", {}).items():
            models = provider_info.get("models", {})
            
            # Check if it's an actual model name
            if model_name in models:
                return model_name, provider_name
            
            # Check if it's an alias
            for actual_name, model_info in models.items():
                aliases = model_info.get("aliases", [])
                if model_name in aliases:
                    return actual_name, provider_name
        
        # If still not found, check if it looks like a provider/model format
        if "/" in model_name:
            return model_name, self._infer_provider_from_model(model_name)
        
        # If nothing found, return as-is with inferred provider
        inferred_provider = self._infer_provider_from_model(model_name) or provider or "openrouter"
        return model_name, inferred_provider
    
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
            provider: Provider name (if None, returns overall best model)
            use_case: Use case (annotation, scoring, annotation_boost, merging)
        
        Returns:
            Tuple of (model_name, provider_name)
        """
        if use_case:
            use_case_info = self.settings.get("use_case_recommendations", {}).get(use_case, {})
            if use_case_info:
                best_model = use_case_info.get("best")
                if best_model:
                    provider_name = self._infer_provider_from_model(best_model)
                    return best_model, provider_name
        
        if provider:
            provider_info = self.settings.get("providers", {}).get(provider, {})
            recommended_model = provider_info.get("recommended_model")
            if recommended_model:
                return recommended_model, provider
        
        # Return overall default - check unified_names first, then provider_shortcuts, then fallback
        unified_names = self.settings.get("unified_names", {})
        if "default" in unified_names:
            default_model = unified_names["default"]
        else:
            # Try to get default from openrouter provider shortcuts as it's the main provider
            provider_shortcuts = self.settings.get("provider_shortcuts", {})
            openrouter_shortcuts = provider_shortcuts.get("openrouter", {})
            default_model = openrouter_shortcuts.get("default", "google/gemini-2.5-flash")
        
        provider_name = self._infer_provider_from_model(default_model)
        return default_model, provider_name
    
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
        provider: Provider name (if None, will try to infer from model_name)
    
    Returns:
        Tuple of (actual_model_name, provider_name)
    """
    return get_model_settings().resolve_model_name(model_name, provider)

def get_recommended_model(provider: str = None, use_case: str = None) -> Tuple[str, str]:
    """
    Get the recommended model for a provider or use case.
    
    Args:
        provider: Provider name (if None, returns overall best model)
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
    if unified:
        for alias, actual in unified.items():
            print(f"  {alias} -> {actual}")
    else:
        # Show provider shortcuts instead
        provider_shortcuts = settings.settings.get("provider_shortcuts", {})
        for provider_name, shortcuts in provider_shortcuts.items():
            print(f"\n  {provider_name.title()} shortcuts:")
            for alias, actual in shortcuts.items():
                print(f"    {alias} -> {actual}")