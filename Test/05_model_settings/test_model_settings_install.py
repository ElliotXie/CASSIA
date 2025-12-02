"""
CASSIA Test 05: Model Settings (PIP INSTALL MODE)
==================================================
Tests model name resolution, provider shortcuts, and fuzzy alias matching
using pip-installed CASSIA.

Usage:
    python test_model_settings_install.py
"""

import sys
import time
import io
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_marker_dataframe_for_cluster
from test_utils import (
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    print_config_summary,
    verify_cassia_pip_install
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA


def run_model_settings_test(results_dir):
    """Test model settings and name resolution using pip-installed CASSIA."""
    print_test_header("05 - Model Settings (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Results directory passed in from main()
    print(f"Results will be saved to: {results_dir['base']}")

    start_time = time.time()
    errors = []
    test_results = {
        "resolution_tests": [],
        "practical_test": None
    }

    # Test 1: Tier shortcuts
    print("\n" + "="*40)
    print("Testing Tier Shortcuts")
    print("="*40)

    tier_tests = [
        ("best", "openai", "gpt"),
        ("balanced", "openai", "gpt-4o"),
        ("fast", "openai", "gpt"),
        ("best", "anthropic", "opus"),
        ("balanced", "anthropic", "sonnet"),
        ("fast", "anthropic", "haiku"),
        ("best", "openrouter", "claude"),
        ("fast", "openrouter", "gemini"),
    ]

    tier_passed = 0
    for tier, provider, expected_contains in tier_tests:
        try:
            resolved_name, resolved_provider = CASSIA.resolve_model_name(tier, provider, verbose=False)
            success = expected_contains.lower() in resolved_name.lower()
            if not success:
                errors.append(f"Tier '{tier}' resolved to '{resolved_name}', expected to contain '{expected_contains}'")

            test_results["resolution_tests"].append({
                "type": "tier",
                "input": tier,
                "provider": provider,
                "resolved": resolved_name,
                "success": success
            })

            status_symbol = "[OK]" if success else "[X]"
            print(f"  {status_symbol} '{tier}' + '{provider}' -> '{resolved_name}'")

            if success:
                tier_passed += 1

        except Exception as e:
            test_results["resolution_tests"].append({
                "type": "tier",
                "input": tier,
                "provider": provider,
                "error": str(e),
                "success": False
            })
            errors.append(f"Tier resolution failed for '{tier}': {e}")
            print(f"  [X] '{tier}' + '{provider}' -> ERROR: {e}")

    print(f"\nTier tests: {tier_passed}/{len(tier_tests)} passed")

    # Test 2: Practical test with model
    print("\n" + "="*40)
    print("Testing Practical Annotation with OpenRouter")
    print("="*40)

    data_config = config['data']
    marker_df = get_marker_dataframe_for_cluster("plasma cell", 20)

    try:
        print(f"  Model: google/gemini-2.5-flash")
        print(f"  Provider: openrouter")
        print(f"  Running annotation...")

        result, _, _ = CASSIA.runCASSIA(
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
    practical_passed = test_results.get("practical_test", {}).get("success", False)
    status = "passed" if practical_passed else "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="model_settings_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=["plasma cell"],
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)
    save_test_results(results_dir['base'], test_results)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("05_model_settings")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_model_settings_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
