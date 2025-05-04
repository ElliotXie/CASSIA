"""
Simple test to verify CASSIA imports and custom provider functionality.
"""

import os
import sys
from CASSIA import (
    set_custom_provider,
    list_custom_providers, 
    runCASSIA
)

def test_cassia_imports():
    """Test basic CASSIA imports"""
    print("CASSIA import test: PASSED")
    return True

def test_custom_provider():
    """Test custom provider functionality with DeepSeek API"""
    # DeepSeek API credentials
    deepseek_api_key = "sk-e967ebf66bc24069a5dfa642792bc491"
    deepseek_model = "deepseek-chat"
    deepseek_url = "https://api.deepseek.com"
    
    # Set up the custom provider
    set_custom_provider(deepseek_api_key, "deepseek", deepseek_url)
    
    # Verify the provider was set up correctly
    providers = list_custom_providers()
    if "deepseek" in providers and providers["deepseek"] == deepseek_url:
        print(f"Custom provider 'deepseek' successfully configured with URL: {providers['deepseek']}")
    else:
        print("Failed to configure custom provider")
        return False
    
    # Simple test markers that should be identifiable by any LLM
    test_markers = ["CD3D", "CD4", "IL7R", "CD3E", "CD3G"]
    
    print("Running simple CASSIA test with DeepSeek custom API...")
    try:
        # Run a simple analysis
        result, _ = runCASSIA(
            model=deepseek_model,
            marker_list=test_markers,
            tissue="blood",
            species="human",
            provider="deepseek"
        )
        
        # Check if we got a valid result
        if isinstance(result, dict) and "main_cell_type" in result:
            print(f"Successfully received cell type annotation: {result['main_cell_type']}")
            print(f"Sub cell types: {result.get('sub_cell_types', [])}")
            return True
        else:
            print(f"Unexpected result format: {result}")
            return False
    except Exception as e:
        print(f"Error running CASSIA with custom provider: {str(e)}")
        return False

if __name__ == "__main__":
    test_cassia_imports()
    test_custom_provider() 