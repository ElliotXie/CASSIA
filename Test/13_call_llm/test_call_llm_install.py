"""
CASSIA Test 13: call_llm Provider Testing (PIP INSTALL MODE)
=============================================================
Tests that call_llm correctly connects to different LLM providers
using pip-installed CASSIA.

Usage:
    python test_call_llm_install.py
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
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA

# Try to load local API keys file if it exists
try:
    import set_api_keys
    print("Loaded API keys from set_api_keys.py")
except ImportError:
    pass


def run_call_llm_test(results_dir):
    """Test call_llm function with different providers using pip-installed CASSIA."""
    print_test_header("13 - call_llm Provider Testing (PIP INSTALL MODE)")

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

    # Test cases
    test_cases = [
        ("openai", "gpt-4o", "OPENAI_API_KEY", "OpenAI"),
        ("anthropic", "claude-sonnet-4-5", "ANTHROPIC_API_KEY", "Anthropic"),
        ("openrouter", "google/gemini-2.5-flash", "OPENROUTER_API_KEY", "OpenRouter"),
    ]

    test_prompt = "Reply with exactly one word: 'success'"

    start_time = time.time()
    errors = []
    test_results = {"provider_tests": []}

    passed_count = 0
    skipped_count = 0
    failed_count = 0
    total_tests = len(test_cases)

    for provider, model, env_var, display_name in test_cases:
        print("\n" + "=" * 50)
        print(f"Provider: {display_name}")
        print(f"  Model: {model}")
        print("=" * 50)

        test_result = {
            "provider": provider,
            "display_name": display_name,
            "model": model,
            "env_var": env_var,
            "status": "pending",
        }

        api_key = os.environ.get(env_var)
        if not api_key:
            print(f"  Status: [SKIP] {env_var} not set")
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
                print(f"  [OK] Response received ({len(response)} chars)")
                print(f"  Response: {response[:80]}{'...' if len(response) > 80 else ''}")
                passed_count += 1
            else:
                test_result["status"] = "failed"
                test_result["error"] = "Empty response"
                print(f"  [X] FAILED: Empty response")
                errors.append(f"{display_name}: Empty response")
                failed_count += 1

        except Exception as e:
            test_result["status"] = "failed"
            test_result["error"] = str(e)[:200]
            print(f"  [X] ERROR: {str(e)[:200]}")
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

    tests_run = total_tests - skipped_count
    if tests_run == 0:
        status = "skipped"
    elif passed_count == tests_run:
        status = "passed"
    else:
        status = "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="call_llm_provider_test_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[],
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    metadata["test_summary"] = {
        "passed": passed_count,
        "failed": failed_count,
        "skipped": skipped_count,
        "total": total_tests
    }
    save_test_metadata(results_dir['base'], metadata)
    save_test_results(results_dir['base'], test_results)

    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    # Create results directory for logging
    results_dir = create_results_dir("13_call_llm")

    # Setup logging to capture all console output
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_call_llm_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
