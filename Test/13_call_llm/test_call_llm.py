"""
CASSIA Test 13: call_llm Provider Testing
==========================================
Tests that call_llm correctly connects to different LLM providers.

This test verifies API connectivity for:
- OpenAI (gpt-4o)
- Anthropic (claude-sonnet-4-5)
- OpenRouter (google/gemini-2.5-flash)
- Custom/DeepSeek (deepseek-chat via https://api.deepseek.com)

Usage:
    python test_call_llm.py
"""

import os
import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

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

# Try to load local API keys file if it exists
try:
    import set_api_keys  # Load API keys from local file
    print("Loaded API keys from set_api_keys.py")
except ImportError:
    pass  # No local keys file, will use setup_api_keys() instead


def run_call_llm_test():
    """Test call_llm function with different providers."""
    print_test_header("13 - call_llm Provider Testing")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Import CASSIA call_llm function
    from llm_utils import call_llm

    # Create results directory
    results_dir = create_results_dir("13_call_llm")
    print(f"Results will be saved to: {results_dir}")

    # Test cases: (provider, model, env_var_name, display_name)
    test_cases = [
        ("openai", "gpt-4o", "OPENAI_API_KEY", "OpenAI"),
        ("anthropic", "claude-sonnet-4-5", "ANTHROPIC_API_KEY", "Anthropic"),
        ("openrouter", "google/gemini-2.5-flash", "OPENROUTER_API_KEY", "OpenRouter"),
        ("https://api.deepseek.com", "deepseek-chat", "CUSTOMIZED_API_KEY", "Custom (DeepSeek)"),
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
            response = call_llm(
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
        test_name="call_llm_provider_test",
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
    save_test_metadata(results_dir, metadata)
    save_test_results(results_dir, test_results)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_call_llm_test()
    sys.exit(0 if success else 1)
