"""
Unit tests for custom provider functionality in CASSIA.
"""

import unittest
import os
from CASSIA import (
    set_custom_provider,
    list_custom_providers,
    get_custom_base_url,
    get_custom_api_key,
    set_custom_base_url,
    set_custom_api_key,
    set_api_key
)

class TestCustomProviders(unittest.TestCase):
    """Test custom provider functionality."""
    
    def setUp(self):
        """Set up test environment before each test."""
        # Clear any environment variables that might affect tests
        for var in ['OPENAI_API_KEY', 'ANTHROPIC_API_KEY', 'OPENROUTER_API_KEY']:
            if var in os.environ:
                del os.environ[var]
        
        # Sample provider data
        self.providers = {
            "provider1": {
                "api_key": "test-api-key-1",
                "base_url": "https://api.provider1.example.com"
            },
            "provider2": {
                "api_key": "test-api-key-2",
                "base_url": "https://api.provider2.example.com"
            }
        }
    
    def test_set_custom_provider(self):
        """Test setting a custom provider with the convenience function."""
        for name, config in self.providers.items():
            result = set_custom_provider(
                api_key=config["api_key"],
                provider_name=name,
                base_url=config["base_url"]
            )
            self.assertIsInstance(result, str)
            self.assertIn(name, result)
            self.assertIn(config["base_url"], result)
    
    def test_list_custom_providers(self):
        """Test listing all custom providers."""
        # First set up some providers
        for name, config in self.providers.items():
            set_custom_provider(
                api_key=config["api_key"],
                provider_name=name,
                base_url=config["base_url"]
            )
        
        # Get the list of providers
        providers_list = list_custom_providers()
        
        # Verify all our providers are in the list
        self.assertIsInstance(providers_list, dict)
        for name, config in self.providers.items():
            self.assertIn(name, providers_list)
            self.assertEqual(providers_list[name], config["base_url"])
    
    def test_get_custom_base_url(self):
        """Test getting a custom base URL."""
        # First set up a provider
        provider_name = "test_provider"
        base_url = "https://api.testprovider.example.com"
        set_custom_base_url(provider_name, base_url)
        
        # Verify we can get the base URL
        retrieved_url = get_custom_base_url(provider_name)
        self.assertEqual(retrieved_url, base_url)
        
        # Test with non-existent provider
        self.assertIsNone(get_custom_base_url("non_existent_provider"))
    
    def test_get_custom_api_key(self):
        """Test getting a custom API key."""
        # First set up a provider
        provider_name = "test_provider"
        api_key = "test-api-key-secret"
        set_custom_api_key(provider_name, api_key)
        
        # Verify we can get the API key
        retrieved_key = get_custom_api_key(provider_name)
        self.assertEqual(retrieved_key, api_key)
        
        # Test with non-existent provider
        self.assertIsNone(get_custom_api_key("non_existent_provider"))
    
    def test_set_api_key(self):
        """Test setting API keys with different providers."""
        # Test with OpenAI
        openai_key = "sk-openai-test-key"
        set_api_key(openai_key, provider="openai")
        self.assertEqual(os.environ.get("OPENAI_API_KEY"), openai_key)
        
        # Test with Anthropic
        anthropic_key = "sk-ant-test-key"
        set_api_key(anthropic_key, provider="anthropic")
        self.assertEqual(os.environ.get("ANTHROPIC_API_KEY"), anthropic_key)
        
        # Test with OpenRouter
        openrouter_key = "sk-or-test-key"
        set_api_key(openrouter_key, provider="openrouter")
        self.assertEqual(os.environ.get("OPENROUTER_API_KEY"), openrouter_key)
        
        # Test with custom provider
        custom_key = "sk-custom-test-key"
        custom_url = "https://api.custom.example.com"
        set_api_key(custom_key, provider="custom", base_url=custom_url)
        self.assertEqual(get_custom_api_key("custom"), custom_key)
        self.assertEqual(get_custom_base_url("custom"), custom_url)
        
        # Test error case: custom provider without base_url
        with self.assertRaises(ValueError):
            set_api_key("test-key", provider="error_provider")
    
    def test_case_insensitivity(self):
        """Test that provider names are case-insensitive."""
        api_key = "test-api-key"
        base_url = "https://api.test.example.com"
        
        # Set with lowercase
        set_custom_provider(api_key, "test_case", base_url)
        
        # Retrieve with different case
        self.assertEqual(get_custom_api_key("TEST_CASE"), api_key)
        self.assertEqual(get_custom_base_url("Test_Case"), base_url)
    
    def test_multiple_providers(self):
        """Test that multiple providers don't interfere with each other."""
        providers = {
            "provider_a": {"key": "key-a", "url": "https://api.a.example.com"},
            "provider_b": {"key": "key-b", "url": "https://api.b.example.com"},
            "provider_c": {"key": "key-c", "url": "https://api.c.example.com"}
        }
        
        # Set up all providers
        for name, config in providers.items():
            set_custom_provider(config["key"], name, config["url"])
        
        # Verify each provider independently
        for name, config in providers.items():
            self.assertEqual(get_custom_api_key(name), config["key"])
            self.assertEqual(get_custom_base_url(name), config["url"])


if __name__ == "__main__":
    unittest.main() 