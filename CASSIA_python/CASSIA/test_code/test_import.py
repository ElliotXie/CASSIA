#!/usr/bin/env python3
"""Test the import from the test_code directory"""

import sys
import os

# Add parent directory to path
parent_dir = os.path.dirname(os.getcwd())
sys.path.insert(0, parent_dir)

print(f"Current directory: {os.getcwd()}")
print(f"Parent directory: {parent_dir}")
print(f"Files in parent directory: {os.listdir(parent_dir)}")

# Test the import
try:
    from model_settings import resolve_model_name, get_model_settings
    print("✅ Import successful")
    
    # Test basic functionality
    settings = get_model_settings()
    print(f"📁 Configuration path: {settings.config_path}")
    
    # Test resolution
    resolved = resolve_model_name('gpt4', 'openai')
    print(f"🧪 Test: gpt4 + openai → {resolved[0]}")
    
    print("🎉 Notebook import is working!")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()