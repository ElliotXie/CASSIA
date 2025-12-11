"""
Tests for provider-specific default models in CASSIA pipeline.
"""

import pytest
import sys
import os

# Add CASSIA to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'CASSIA_python'))

from CASSIA.core.model_settings import get_pipeline_defaults, get_model_settings


class TestPipelineDefaults:
    """Tests for get_pipeline_defaults function."""

    def test_openai_defaults(self):
        """Test OpenAI provider defaults."""
        defaults = get_pipeline_defaults("openai")

        assert defaults["annotation"] == "gpt-5.1"
        assert defaults["score"] == "gpt-5-mini"
        assert defaults["merge"] == "gpt-5-mini"
        assert defaults["annotationboost"] == "gpt-5.1"

    def test_anthropic_defaults(self):
        """Test Anthropic provider defaults."""
        defaults = get_pipeline_defaults("anthropic")

        assert defaults["annotation"] == "claude-sonnet-4-5"
        assert defaults["score"] == "claude-sonnet-4-5"
        assert defaults["merge"] == "claude-haiku-4-5"
        assert defaults["annotationboost"] == "claude-sonnet-4-5"

    def test_openrouter_defaults(self):
        """Test OpenRouter provider defaults."""
        defaults = get_pipeline_defaults("openrouter")

        assert defaults["annotation"] == "openai/gpt-5.1"
        assert defaults["score"] == "anthropic/claude-sonnet-4.5"
        assert defaults["merge"] == "google/gemini-2.5-flash"
        assert defaults["annotationboost"] == "anthropic/claude-sonnet-4.5"

    def test_case_insensitive(self):
        """Test that provider names are case-insensitive."""
        defaults_lower = get_pipeline_defaults("openai")
        defaults_upper = get_pipeline_defaults("OPENAI")
        defaults_mixed = get_pipeline_defaults("OpenAI")

        assert defaults_lower == defaults_upper
        assert defaults_lower == defaults_mixed

    def test_unknown_provider_fallback(self):
        """Test that unknown provider falls back to openrouter defaults."""
        defaults = get_pipeline_defaults("unknown_provider")
        openrouter_defaults = get_pipeline_defaults("openrouter")

        assert defaults == openrouter_defaults

    def test_all_keys_present(self):
        """Test that all expected keys are present for each provider."""
        expected_keys = {"annotation", "score", "merge", "annotationboost"}

        for provider in ["openai", "anthropic", "openrouter"]:
            defaults = get_pipeline_defaults(provider)
            assert set(defaults.keys()) >= expected_keys, f"Missing keys for {provider}"

    def test_model_settings_instance_method(self):
        """Test that ModelSettings instance method works."""
        settings = get_model_settings()

        defaults = settings.get_pipeline_defaults("openai")
        assert defaults["annotation"] == "gpt-5.1"

    def test_defaults_not_empty(self):
        """Test that no default values are empty strings."""
        for provider in ["openai", "anthropic", "openrouter"]:
            defaults = get_pipeline_defaults(provider)
            for key, value in defaults.items():
                assert value, f"Empty value for {provider}.{key}"


class TestPipelineDefaultsIntegration:
    """Integration tests for pipeline defaults with model resolution."""

    def test_defaults_are_valid_models(self):
        """Test that default models are valid and can be resolved."""
        from CASSIA.core.model_settings import resolve_model_name

        for provider in ["openai", "anthropic"]:
            defaults = get_pipeline_defaults(provider)
            for stage, model in defaults.items():
                # Should not raise an error
                resolved_model, resolved_provider = resolve_model_name(model, provider, verbose=False)
                assert resolved_model, f"Failed to resolve {model} for {provider}"

    def test_openrouter_defaults_have_namespace(self):
        """Test that OpenRouter defaults include provider namespace."""
        defaults = get_pipeline_defaults("openrouter")

        for key, value in defaults.items():
            assert "/" in value, f"OpenRouter default {key}={value} missing provider namespace"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
