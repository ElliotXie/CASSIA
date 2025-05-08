import pytest
import os
from ..tools_function import (
    set_openai_api_key,
    set_anthropic_api_key,
    set_openrouter_api_key,
    set_api_key
)

def test_set_openai_api_key():
    test_key = "test-openai-key-123"
    set_openai_api_key(test_key)
    assert os.environ["OPENAI_API_KEY"] == test_key

def test_set_anthropic_api_key():
    test_key = "test-anthropic-key-456"
    set_anthropic_api_key(test_key)
    assert os.environ["ANTHROPIC_API_KEY"] == test_key

def test_set_openrouter_api_key():
    test_key = "test-openrouter-key-789"
    set_openrouter_api_key(test_key)
    assert os.environ["OPENROUTER_API_KEY"] == test_key

def test_set_api_key():
    # Test setting OpenAI key
    openai_key = "openai-via-general-function"
    set_api_key(openai_key, provider="openai")
    assert os.environ["OPENAI_API_KEY"] == openai_key
    
    # Test setting Anthropic key
    anthropic_key = "anthropic-via-general-function"
    set_api_key(anthropic_key, provider="anthropic")
    assert os.environ["ANTHROPIC_API_KEY"] == anthropic_key
    
    # Test setting OpenRouter key
    openrouter_key = "openrouter-via-general-function"
    set_api_key(openrouter_key, provider="openrouter")
    assert os.environ["OPENROUTER_API_KEY"] == openrouter_key
    
    # Test invalid provider
    with pytest.raises(ValueError, match="Provider must be either 'openai' or 'anthropic' or 'openrouter'"):
        set_api_key("some-key", provider="invalid-provider") 