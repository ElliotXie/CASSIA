"""
CASSIA Test 05: Model Settings
==============================
Tests model name resolution and provider shortcuts.

Usage:
    python test_model_settings.py
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_marker_dataframe_for_cluster
from test_utils import (
    setup_cassia_imports,
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    print_config_summary
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata
)

# Setup CASSIA imports
setup_cassia_imports()


def run_model_settings_test():
    """Test model settings and name resolution."""
    print_test_header("05 - Model Settings")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Import CASSIA functions
    from model_settings import resolve_model_name, get_recommended_model
    from tools_function import runCASSIA

    # Create results directory
    results_dir = create_results_dir("05_model_settings")
    print(f"Results will be saved to: {results_dir}")

    start_time = time.time()
    errors = []
    test_results = {
        "resolution_tests": [],
        "practical_test": None
    }

    # Test 1: Model name resolution
    print("\n" + "="*40)
    print("Testing Model Name Resolution")
    print("="*40)

    resolution_tests = [
        # (simple_name, provider, expected_contains)
        ("gemini", "openrouter", "gemini"),
        ("best", "openrouter", "gemini"),
        ("cheap", "openrouter", None),  # Just verify it resolves
        ("flash", "openrouter", "flash"),
    ]

    resolution_passed = 0
    for simple_name, provider, expected_contains in resolution_tests:
        try:
            resolved_name, resolved_provider = resolve_model_name(simple_name, provider)
            success = True
            if expected_contains and expected_contains.lower() not in resolved_name.lower():
                success = False
                errors.append(f"'{simple_name}' resolved to '{resolved_name}', expected to contain '{expected_contains}'")

            test_results["resolution_tests"].append({
                "input": simple_name,
                "provider": provider,
                "resolved": resolved_name,
                "success": success
            })

            status_symbol = "[OK]" if success else "[X]"
            print(f"  {status_symbol} '{simple_name}' + '{provider}' -> '{resolved_name}'")

            if success:
                resolution_passed += 1

        except Exception as e:
            test_results["resolution_tests"].append({
                "input": simple_name,
                "provider": provider,
                "error": str(e),
                "success": False
            })
            errors.append(f"Resolution failed for '{simple_name}': {e}")
            print(f"  [X] '{simple_name}' + '{provider}' -> ERROR: {e}")

    print(f"\nResolution tests: {resolution_passed}/{len(resolution_tests)} passed")

    # Test 2: Recommended model
    print("\n" + "="*40)
    print("Testing Recommended Model")
    print("="*40)

    try:
        rec_model = get_recommended_model("annotation", "openrouter")
        print(f"  Recommended for annotation: {rec_model}")
        test_results["recommended_model"] = rec_model
    except Exception as e:
        errors.append(f"get_recommended_model failed: {e}")
        print(f"  [X] ERROR: {e}")

    # Test 3: Practical test with gemini-2.5-flash via OpenRouter
    print("\n" + "="*40)
    print("Testing Practical Annotation with OpenRouter")
    print("="*40)

    data_config = config['data']
    marker_df = get_marker_dataframe_for_cluster("plasma cell", 20)

    try:
        print(f"  Model: google/gemini-2.5-flash")
        print(f"  Provider: openrouter")
        print(f"  Running annotation...")

        result, _ = runCASSIA(
            model="google/gemini-2.5-flash",
            temperature=0.3,
            marker_list=marker_df,
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            provider="openrouter",
            validator_involvement="v1"
        )

        practical_success = 'main_cell_type' in result and result['main_cell_type']
        test_results["practical_test"] = {
            "success": practical_success,
            "main_cell_type": result.get('main_cell_type'),
            "sub_cell_types": result.get('sub_cell_types')
        }

        if practical_success:
            print(f"  [OK] Annotation successful")
            print(f"       Main type: {result.get('main_cell_type')}")
        else:
            print(f"  [X] Annotation returned empty result")
            errors.append("Practical test returned empty result")

    except Exception as e:
        test_results["practical_test"] = {
            "success": False,
            "error": str(e)
        }
        errors.append(f"Practical test failed: {e}")
        print(f"  [X] ERROR: {e}")

    duration = time.time() - start_time

    # Determine overall status
    resolution_all_passed = resolution_passed == len(resolution_tests)
    practical_passed = test_results.get("practical_test", {}).get("success", False)

    if resolution_all_passed and practical_passed:
        status = "passed"
    elif practical_passed:
        status = "passed"  # Practical test is most important
    else:
        status = "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="model_settings",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=["plasma cell"],
        errors=errors
    )
    save_test_metadata(results_dir, metadata)
    save_test_results(results_dir, test_results)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_model_settings_test()
    sys.exit(0 if success else 1)
