#!/usr/bin/env python3
"""
Test script for CASSIA batch analysis with model settings
Demonstrates using simple model names like 'gemini' and 'claude'
"""

import sys
import os
import pandas as pd
import numpy as np

# Add parent directory to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

# Import CASSIA and model settings
import CASSIA
from model_settings import resolve_model_name

def create_sample_data():
    """Create sample marker data for testing"""
    print("📊 Creating sample marker data...")
    
    # Sample marker data with realistic gene expression values
    sample_markers = pd.DataFrame({
        'Gene': ['CD3D', 'CD4', 'CD8A', 'CD19', 'CD14', 'FCGR3A', 'IL7R', 'KLRB1', 'GNLY', 'NKG7'],
        'Cluster_0': [0.1, 0.2, 0.05, 0.01, 0.8, 0.1, 0.3, 0.02, 0.01, 0.02],
        'Cluster_1': [0.9, 0.8, 0.1, 0.02, 0.1, 0.05, 0.7, 0.03, 0.02, 0.03],
        'Cluster_2': [0.8, 0.2, 0.9, 0.01, 0.05, 0.02, 0.4, 0.01, 0.02, 0.01],
        'Cluster_3': [0.2, 0.1, 0.05, 0.9, 0.02, 0.01, 0.1, 0.05, 0.01, 0.02],
        'Cluster_4': [0.1, 0.05, 0.02, 0.05, 0.9, 0.7, 0.2, 0.1, 0.05, 0.1],
        'Cluster_5': [0.05, 0.02, 0.01, 0.02, 0.3, 0.8, 0.1, 0.6, 0.9, 0.8]
    })
    
    print(f"✅ Created sample data: {sample_markers.shape[0]} genes, {sample_markers.shape[1]-1} clusters")
    return sample_markers

def test_model_resolution():
    """Test model name resolution"""
    print("\n🔧 Testing model name resolution...")
    
    test_cases = [
        ("gemini", "openrouter"),
        ("claude", "anthropic"),
        ("cheap", "openrouter"),
        ("fast", "anthropic"),
        ("best", "openai")
    ]
    
    for model, provider in test_cases:
        try:
            resolved = resolve_model_name(model, provider)
            print(f"✅ {model:8} + {provider:10} → {resolved[0]}")
        except Exception as e:
            print(f"❌ {model:8} + {provider:10} → Error: {e}")

def run_batch_analysis(model_name, provider, output_name):
    """Run CASSIA batch analysis with specified model settings"""
    print(f"\n🚀 Running batch analysis with {model_name} ({provider})...")
    
    # Create sample data
    unprocessed_markers = create_sample_data()
    
    # Show model resolution
    try:
        resolved = resolve_model_name(model_name, provider)
        print(f"🔧 Model resolution: {model_name} + {provider} → {resolved[0]}")
    except Exception as e:
        print(f"❌ Model resolution failed: {e}")
        return
    
    # Prepare batch analysis parameters
    print(f"\n📋 Batch analysis configuration:")
    print(f"   Output name: {output_name}")
    print(f"   Model: {model_name}")
    print(f"   Provider: {provider}")
    print(f"   Tissue: large intestine")
    print(f"   Species: human")
    print(f"   Max workers: 6")
    print(f"   N genes: 50")
    
    # Check if API key is set
    api_key_vars = {
        'openai': 'OPENAI_API_KEY',
        'anthropic': 'ANTHROPIC_API_KEY',
        'openrouter': 'OPENROUTER_API_KEY'
    }
    
    api_key_var = api_key_vars.get(provider)
    if api_key_var and not os.environ.get(api_key_var):
        print(f"⚠️  Warning: {api_key_var} not set in environment")
        print(f"   To run the analysis, set: export {api_key_var}='your_api_key'")
        print("   Skipping actual API call for now...")
        return
    
    # Run the batch analysis
    try:
        print("🚀 Starting CASSIA batch analysis...")
        
        result = CASSIA.runCASSIA_batch(
            marker=unprocessed_markers,
            output_name=output_name,
            model=model_name,
            tissue="large intestine",
            species="human",
            max_workers=6,
            n_genes=50,
            additional_info=None,
            provider=provider
        )
        
        print("✅ Batch analysis completed successfully!")
        print(f"📊 Results saved to: {output_name}")
        
    except Exception as e:
        print(f"❌ Error during batch analysis: {e}")
        if "API key" in str(e):
            print(f"   Make sure {api_key_var} is set correctly")

def main():
    """Main function to run tests"""
    print("🧪 CASSIA Batch Analysis with Model Settings Test")
    print("=" * 50)
    
    # Test model resolution
    test_model_resolution()
    
    # Test with different models
    test_configs = [
        ("gemini", "openrouter", "test_gemini_analysis"),
        ("cheap", "openrouter", "test_cheap_analysis"),
        ("claude", "anthropic", "test_claude_analysis"),
        ("fast", "anthropic", "test_fast_analysis")
    ]
    
    for model, provider, output in test_configs:
        run_batch_analysis(model, provider, output)
        print("-" * 50)

if __name__ == "__main__":
    main()