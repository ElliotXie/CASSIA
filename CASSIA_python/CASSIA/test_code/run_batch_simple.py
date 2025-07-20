#!/usr/bin/env python3
"""
Simple script to run CASSIA batch analysis with model settings
Usage: python run_batch_simple.py
"""

import sys
import os
import pandas as pd
import numpy as np

# Add parent directory to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

# Direct imports from local files, not from the installed package
try:
    from .tools_function import *
except ImportError:
    from tools_function import *

try:
    from .main_function_code import *
except ImportError:
    from main_function_code import *

try:
    from .model_settings import *
except ImportError:
    from model_settings import *

try:
    from .llm_utils import call_llm
except ImportError:
    from llm_utils import call_llm

print("✅ Successfully imported CASSIA modules locally")

# Set API key for OpenRouter
os.environ['OPENROUTER_API_KEY'] = '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'
print("🔑 OpenRouter API key set")

def load_real_data():
    """Load the real CASSIA dataset"""
    # Load the unprocessed dataset from CASSIA data directory
    data_path = os.path.join(parent_dir, "data", "unprocessed.csv")
    print(f"📁 Loading dataset from: {data_path}")
    
    try:
        df = pd.read_csv(data_path)
        print(f"✅ Loaded dataset: {df.shape[0]} rows, {df.shape[1]} columns")
        
        # Show basic info about the dataset
        if 'cluster' in df.columns:
            clusters = df['cluster'].unique()
            print(f"📊 Clusters found: {list(clusters)}")
        if 'gene' in df.columns:
            genes = df['gene'].nunique()
            print(f"🧬 Unique genes: {genes}")
            
        return df
        
    except FileNotFoundError:
        print(f"❌ Dataset not found at {data_path}")
        print("📝 Make sure you're running from the correct directory")
        return None
    except Exception as e:
        print(f"❌ Error loading dataset: {e}")
        return None

def main():
    print("🚀 Running CASSIA batch analysis with model settings...")
    
    # Load real data
    unprocessed_markers = load_real_data()
    if unprocessed_markers is None:
        print("❌ Failed to load dataset. Exiting.")
        return
    
    output_name = "model_settings_test_real_data"
    
    print(f"📊 Real dataset: {unprocessed_markers.shape[0]} rows, {unprocessed_markers.shape[1]} columns")
    
    # Test with gemini (OpenRouter) - using resolved model name
    print("\n🔧 Testing with gemini (OpenRouter)...")
    try:
        resolved_model = resolve_model_name("gemini", "openrouter")
        print(f"🔧 Resolved: gemini → {resolved_model[0]}")
        
        result = runCASSIA_batch(
            marker=unprocessed_markers,
            output_name=f"{output_name}_gemini",
            model=resolved_model[0],  # Use resolved model name
            tissue="large intestine",
            species="human",
            max_workers=6,
            n_genes=50,
            additional_info=None,
            provider="openrouter"
        )
        print("✅ Gemini analysis completed!")
    except Exception as e:
        print(f"❌ Gemini analysis failed: {e}")
        if "API key" in str(e):
            print("💡 Set OPENROUTER_API_KEY environment variable")
    
    # Test with cheap (OpenRouter - DeepSeek) - using resolved model name
    print("\n🔧 Testing with cheap (OpenRouter - DeepSeek)...")
    try:
        resolved_model = resolve_model_name("cheap", "openrouter")
        print(f"🔧 Resolved: cheap → {resolved_model[0]}")
        
        result = runCASSIA_batch(
            marker=unprocessed_markers,
            output_name=f"{output_name}_cheap",
            model=resolved_model[0],  # Use resolved model name
            tissue="large intestine",
            species="human",
            max_workers=6,
            n_genes=50,
            additional_info=None,
            provider="openrouter"
        )
        print("✅ Cheap analysis completed!")
    except Exception as e:
        print(f"❌ Cheap analysis failed: {e}")
        if "API key" in str(e):
            print("💡 Set OPENROUTER_API_KEY environment variable")

if __name__ == "__main__":
    main()