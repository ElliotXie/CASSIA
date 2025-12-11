"""
CASSIA Test 14: Customized Provider Testing (PIP INSTALL MODE)
===============================================================
Tests that call_llm correctly connects to customized/custom API providers.
Uses pip-installed CASSIA package.

This test demonstrates how to use a custom OpenAI-compatible API endpoint
with CASSIA, using DeepSeek as an example.

Usage:
    python test_customized_provider_install.py
"""

import os
import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

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
    create_test_metadata
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA

# Try to load local API keys file if it exists
try:
    import set_api_keys  # Load API keys from local file
    print("Loaded API keys from set_api_keys.py")
except ImportError:
    pass  # No local keys file, will use setup_api_keys() instead


def run_customized_provider_test():
    """Test call_llm function with customized/custom API providers."""
    print_test_header("14 - Customized Provider Testing (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Location: {pip_info['location']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    if not pip_info['is_pip_install']:
        print("\nWARNING: CASSIA may not be from pip install!")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Create results directory
    results_dirs = create_results_dir("14_customized_provider")
    results_dir = results_dirs['base']
    print(f"Results will be saved to: {results_dir}")

    # Test cases for customized providers: (provider_url, model, env_var_name, display_name)
    # DeepSeek is used as an example of a custom OpenAI-compatible API
    test_cases = [
        ("https://api.deepseek.com", "deepseek-chat", "CUSTOMIZED_API_KEY", "DeepSeek Custom API"),
    ]

    # Simple test prompt
    test_prompt = "Reply with exactly one word: 'success'"

    start_time = time.time()
    errors = []
    test_results = {
        "provider_tests": []
    }

    passed_count = 0
    skipped_count = 0
    failed_count = 0
    total_tests = len(test_cases)

    for provider, model, env_var, display_name in test_cases:
        print("\n" + "=" * 50)
        print(f"Provider: {display_name}")
        print(f"  Endpoint: {provider}")
        print(f"  Model: {model}")
        print("=" * 50)

        test_result = {
            "provider": provider,
            "display_name": display_name,
            "model": model,
            "env_var": env_var,
            "status": "pending",
            "response_length": 0,
            "error": None
        }

        # Check if API key is set
        api_key = os.environ.get(env_var)
        if not api_key:
            print(f"  Status: [SKIP] {env_var} not set in environment")
            test_result["status"] = "skipped"
            skipped_count += 1
            test_results["provider_tests"].append(test_result)
            continue

        print(f"  Status: API key found")
        print(f"  Sending test prompt...")

        try:
            response = CASSIA.call_llm(
                prompt=test_prompt,
                provider=provider,
                model=model,
                temperature=0.3,
                max_tokens=50
            )

            if response and len(response) > 0:
                test_result["status"] = "passed"
                test_result["response_length"] = len(response)
                test_result["response_preview"] = response[:100] if len(response) > 100 else response
                print(f"  [OK] Response received ({len(response)} chars)")
                print(f"  Response: {response[:80]}{'...' if len(response) > 80 else ''}")
                passed_count += 1
            else:
                test_result["status"] = "failed"
                test_result["error"] = "Empty response received"
                print(f"  [X] FAILED: Empty response received")
                errors.append(f"{display_name}: Empty response")
                failed_count += 1

        except Exception as e:
            test_result["status"] = "failed"
            test_result["error"] = str(e)
            error_msg = str(e)
            # Truncate long error messages for display
            if len(error_msg) > 200:
                error_msg = error_msg[:200] + "..."
            print(f"  [X] ERROR: {error_msg}")
            errors.append(f"{display_name}: {str(e)[:100]}")
            failed_count += 1

        test_results["provider_tests"].append(test_result)

    duration = time.time() - start_time

    # Summary
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"  Passed:  {passed_count}/{total_tests}")
    print(f"  Failed:  {failed_count}/{total_tests}")
    print(f"  Skipped: {skipped_count}/{total_tests}")
    print(f"  Duration: {duration:.2f}s")

    # Determine overall status
    # Test passes if all non-skipped tests passed
    tests_run = total_tests - skipped_count
    if tests_run == 0:
        status = "skipped"
        print("\n  Note: All tests skipped (no API keys configured)")
    elif passed_count == tests_run:
        status = "passed"
    else:
        status = "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="customized_provider_test_install",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[],  # Not applicable for this test
        errors=errors
    )
    metadata["test_summary"] = {
        "passed": passed_count,
        "failed": failed_count,
        "skipped": skipped_count,
        "total": total_tests
    }
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir, metadata)
    save_test_results(results_dir, test_results)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_customized_provider_test()
    sys.exit(0 if success else 1)
