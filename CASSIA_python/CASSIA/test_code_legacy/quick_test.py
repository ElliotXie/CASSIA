#!/usr/bin/env python3
"""
Quick test script to verify model settings functionality before running the notebook.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.getcwd()), 'CASSIA'))

def quick_test():
    """Run a quick test of model settings functionality"""
    
    print("🧪 Quick Model Settings Test")
    print("=" * 40)
    
    # Test 1: Import
    try:
        from model_settings import (
            resolve_model_name,
            get_recommended_model,
            get_model_settings
        )
        print("✅ Import successful")
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    
    # Test 2: Configuration location
    try:
        settings = get_model_settings()
        print(f"✅ Configuration loaded from: {settings.config_path}")
        
        if "data/model_settings.json" in str(settings.config_path):
            print("✅ Using package data directory")
        else:
            print("⚠️  Using alternative location")
    except Exception as e:
        print(f"❌ Configuration error: {e}")
        return False
    
    # Test 3: Basic resolution
    try:
        resolved = resolve_model_name("gpt4", "openai")
        print(f"✅ Resolution test: gpt4 + openai → {resolved[0]}")
    except Exception as e:
        print(f"❌ Resolution failed: {e}")
        return False
    
    # Test 4: Provider requirement
    try:
        resolve_model_name("gpt4")
        print("❌ Provider requirement test failed")
        return False
    except Exception as e:
        print("✅ Provider requirement enforced")
    
    # Test 5: Recommendations
    try:
        recommended = get_recommended_model("openai")
        print(f"✅ Recommendation test: openai → {recommended[0]}")
    except Exception as e:
        print(f"❌ Recommendation failed: {e}")
        return False
    
    print("\n🎉 All tests passed! Notebook is ready to run.")
    return True

if __name__ == "__main__":
    success = quick_test()
    if not success:
        print("\n❌ Some tests failed. Check the errors above.")
        sys.exit(1)
    else:
        print("\n📚 Run the notebook: test_model_settings.ipynb")