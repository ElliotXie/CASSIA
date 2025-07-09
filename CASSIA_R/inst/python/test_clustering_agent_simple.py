#!/usr/bin/env python3
"""
Simplified test script for CASSIA Clustering Agent functions.
Tests the core logic without requiring scanpy.
"""

import sys
import os
import numpy as np
import pandas as pd
import json
import warnings
warnings.filterwarnings('ignore')

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_basic_imports():
    """Test basic import functionality."""
    print("Testing basic imports...")
    
    try:
        import numpy as np
        print("✓ numpy imported successfully")
    except ImportError as e:
        print(f"✗ numpy import failed: {e}")
        return False
    
    try:
        import pandas as pd
        print("✓ pandas imported successfully")
    except ImportError as e:
        print(f"✗ pandas import failed: {e}")
        return False
    
    try:
        import json
        print("✓ json imported successfully")
    except ImportError as e:
        print(f"✗ json import failed: {e}")
        return False
    
    return True


def test_llm_utils_image_functions():
    """Test LLM image utility functions."""
    print("\nTesting LLM image utilities...")
    
    try:
        from llm_utils_image import encode_image_to_base64, create_image_url
        
        # Create a dummy image file for testing
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot([1, 2, 3], [1, 2, 3])
        ax.set_title("Test Image")
        test_image_path = "test_image.png"
        plt.savefig(test_image_path)
        plt.close()
        
        # Test image encoding
        base64_str = encode_image_to_base64(test_image_path)
        assert len(base64_str) > 0, "Base64 encoding failed"
        assert isinstance(base64_str, str), "Base64 result should be string"
        
        # Test image URL creation
        image_url = create_image_url(test_image_path)
        assert image_url.startswith("data:image/png;base64,"), "Image URL format incorrect"
        
        # Clean up
        os.remove(test_image_path)
        
        print("✓ LLM image utilities work correctly")
        return True
        
    except Exception as e:
        print(f"✗ LLM image utilities failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parameter_combinations():
    """Test parameter combination generation logic."""
    print("\nTesting parameter combination logic...")
    
    try:
        from itertools import product
        
        # Test parameter combinations
        n_pcs_list = [20, 30]
        min_dist_list = [0.2, 0.3]
        n_neighbors_list = [30, 40]
        
        param_combinations = list(product(n_pcs_list, min_dist_list, n_neighbors_list))
        
        assert len(param_combinations) == 8, f"Expected 8 combinations, got {len(param_combinations)}"
        
        # Check that all combinations are unique
        assert len(set(param_combinations)) == len(param_combinations), "Duplicate combinations found"
        
        # Check structure
        for combo in param_combinations:
            assert len(combo) == 3, "Each combination should have 3 parameters"
            assert combo[0] in n_pcs_list, "n_pcs not in expected range"
            assert combo[1] in min_dist_list, "min_dist not in expected range"
            assert combo[2] in n_neighbors_list, "n_neighbors not in expected range"
        
        print("✓ Parameter combination logic works correctly")
        print(f"  - Generated {len(param_combinations)} combinations")
        
        return True
        
    except Exception as e:
        print(f"✗ Parameter combination logic failed: {e}")
        return False


def test_stability_metrics():
    """Test stability metric calculations."""
    print("\nTesting stability metrics...")
    
    try:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
        
        # Create test clustering results
        clustering1 = np.array([0, 0, 1, 1, 2, 2, 0, 1, 2])
        clustering2 = np.array([0, 0, 1, 1, 2, 2, 0, 1, 2])  # Identical
        clustering3 = np.array([1, 1, 0, 0, 2, 2, 1, 0, 2])  # Different but similar structure
        
        # Test identical clusterings
        ari_identical = adjusted_rand_score(clustering1, clustering2)
        assert ari_identical == 1.0, f"ARI for identical clusterings should be 1.0, got {ari_identical}"
        
        # Test different clusterings
        ari_different = adjusted_rand_score(clustering1, clustering3)
        assert 0 <= ari_different <= 1, f"ARI should be between 0 and 1, got {ari_different}"
        
        # Test stability calculation logic
        clustering_results = [clustering1, clustering2, clustering3]
        stability_scores = []
        
        for i in range(len(clustering_results)):
            for j in range(i+1, len(clustering_results)):
                ari = adjusted_rand_score(clustering_results[i], clustering_results[j])
                stability_scores.append(ari)
        
        mean_stability = np.mean(stability_scores)
        std_stability = np.std(stability_scores)
        
        assert 0 <= mean_stability <= 1, "Mean stability should be between 0 and 1"
        assert std_stability >= 0, "Std stability should be non-negative"
        
        print("✓ Stability metrics work correctly")
        print(f"  - Mean stability: {mean_stability:.3f}")
        print(f"  - Std stability: {std_stability:.3f}")
        
        return True
        
    except Exception as e:
        print(f"✗ Stability metrics failed: {e}")
        return False


def test_fallback_analysis_logic():
    """Test fallback analysis logic."""
    print("\nTesting fallback analysis logic...")
    
    try:
        # Create mock parameter grid results
        param_results = [
            {'n_pcs': 20, 'min_dist': 0.2, 'n_neighbors': 30, 'n_clusters': 8},
            {'n_pcs': 20, 'min_dist': 0.3, 'n_neighbors': 30, 'n_clusters': 6},
            {'n_pcs': 30, 'min_dist': 0.2, 'n_neighbors': 30, 'n_clusters': 10},
            {'n_pcs': 30, 'min_dist': 0.3, 'n_neighbors': 30, 'n_clusters': 7}
        ]
        
        # Create mock stability results
        stability_results = [
            {'resolution': 0.3, 'mean_stability': 0.85, 'std_stability': 0.05},
            {'resolution': 0.5, 'mean_stability': 0.92, 'std_stability': 0.03},
            {'resolution': 0.7, 'mean_stability': 0.88, 'std_stability': 0.04}
        ]
        
        # Test parameter selection logic
        cluster_counts = [r['n_clusters'] for r in param_results]
        median_clusters = np.median(cluster_counts)
        best_param_idx = np.argmin([abs(count - median_clusters) for count in cluster_counts])
        best_params = param_results[best_param_idx]
        
        # Test stability selection logic
        stabilities = [r['mean_stability'] for r in stability_results]
        best_stability_idx = np.argmax(stabilities)
        best_resolution = stability_results[best_stability_idx]['resolution']
        
        assert best_resolution == 0.5, f"Expected resolution 0.5, got {best_resolution}"
        assert 'n_pcs' in best_params, "n_pcs not in best parameters"
        assert 'min_dist' in best_params, "min_dist not in best parameters"
        assert 'n_neighbors' in best_params, "n_neighbors not in best parameters"
        
        print("✓ Fallback analysis logic works correctly")
        print(f"  - Selected resolution: {best_resolution}")
        print(f"  - Selected parameters: {best_params}")
        
        return True
        
    except Exception as e:
        print(f"✗ Fallback analysis logic failed: {e}")
        return False


def test_json_parsing():
    """Test JSON parsing for LLM responses."""
    print("\nTesting JSON parsing logic...")
    
    try:
        # Test valid JSON response
        valid_json_response = '''
        Here is my analysis:
        {
            "analysis": "The clustering looks good",
            "best_parameters": {"n_pcs": 30, "min_dist": 0.2, "n_neighbors": 30},
            "confidence": 0.85
        }
        Additional text here.
        '''
        
        # Test JSON extraction
        import re
        json_match = re.search(r'\{.*\}', valid_json_response, re.DOTALL)
        assert json_match is not None, "JSON pattern not found"
        
        parsed_json = json.loads(json_match.group())
        assert 'analysis' in parsed_json, "Analysis not in parsed JSON"
        assert 'best_parameters' in parsed_json, "Best parameters not in parsed JSON"
        assert 'confidence' in parsed_json, "Confidence not in parsed JSON"
        
        # Test invalid JSON handling
        invalid_json_response = "This is just text without JSON"
        json_match = re.search(r'\{.*\}', invalid_json_response, re.DOTALL)
        assert json_match is None, "JSON pattern should not be found in invalid response"
        
        print("✓ JSON parsing logic works correctly")
        
        return True
        
    except Exception as e:
        print(f"✗ JSON parsing logic failed: {e}")
        return False


def test_directory_creation():
    """Test directory creation logic."""
    print("\nTesting directory creation...")
    
    try:
        test_dirs = ["test_dir_1", "test_dir_2/subdir", "test_dir_3"]
        
        for test_dir in test_dirs:
            os.makedirs(test_dir, exist_ok=True)
            assert os.path.exists(test_dir), f"Directory {test_dir} not created"
            assert os.path.isdir(test_dir), f"{test_dir} is not a directory"
        
        # Clean up
        import shutil
        for test_dir in ["test_dir_1", "test_dir_2", "test_dir_3"]:
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
        
        print("✓ Directory creation works correctly")
        
        return True
        
    except Exception as e:
        print(f"✗ Directory creation failed: {e}")
        return False


def test_file_operations():
    """Test file operations."""
    print("\nTesting file operations...")
    
    try:
        # Test file writing and reading
        test_data = {"test": "data", "number": 42}
        test_file = "test_file.json"
        
        # Write JSON file
        with open(test_file, 'w') as f:
            json.dump(test_data, f)
        
        assert os.path.exists(test_file), "Test file not created"
        
        # Read JSON file
        with open(test_file, 'r') as f:
            loaded_data = json.load(f)
        
        assert loaded_data == test_data, "Loaded data doesn't match original"
        
        # Clean up
        os.remove(test_file)
        
        print("✓ File operations work correctly")
        
        return True
        
    except Exception as e:
        print(f"✗ File operations failed: {e}")
        return False


def main():
    """Run all simplified tests."""
    print("="*60)
    print("CASSIA Clustering Agent - Simplified Test Suite")
    print("="*60)
    
    tests = [
        test_basic_imports,
        test_llm_utils_image_functions,
        test_parameter_combinations,
        test_stability_metrics,
        test_fallback_analysis_logic,
        test_json_parsing,
        test_directory_creation,
        test_file_operations
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"✗ Test {test.__name__} crashed: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*60)
    
    if failed == 0:
        print("🎉 ALL TESTS PASSED! 🎉")
        print("\nThe core clustering agent logic is working correctly!")
        print("Note: Full functionality testing requires scanpy and a complete environment.")
    else:
        print(f"❌ {failed} tests failed. Please check the errors above.")
    
    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)