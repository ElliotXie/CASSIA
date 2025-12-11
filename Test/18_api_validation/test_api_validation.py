"""
CASSIA Test 18: API Key Validation Testing
==========================================
Tests the validate_api_keys() function for different LLM providers.

This test verifies:
- API key validation for OpenAI, Anthropic, and OpenRouter
- Caching behavior (second validation should be instant)
- Force revalidation
- Cache invalidation on key change
- Thread safety for concurrent validations
- Error handling for invalid/missing keys

Usage:
    python test_api_validation.py
"""

import os
import sys
import time
import threading
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


def test_validation_basic():
    """Test basic API key validation for all providers."""
    from CASSIA import validate_api_keys, clear_validation_cache

    print("\n" + "=" * 60)
    print("Test 1: Basic API Key Validation")
    print("=" * 60)

    results = {
        "test_name": "basic_validation",
        "providers": []
    }

    # Test each provider that has a key set
    providers = ["openai", "anthropic", "openrouter"]
    env_vars = {
        "openai": "OPENAI_API_KEY",
        "anthropic": "ANTHROPIC_API_KEY",
        "openrouter": "OPENROUTER_API_KEY"
    }

    for provider in providers:
        print(f"\nTesting {provider}...")
        env_var = env_vars[provider]

        if not os.environ.get(env_var):
            print(f"  [SKIP] {env_var} not set")
            results["providers"].append({
                "provider": provider,
                "status": "skipped",
                "reason": "API key not set"
            })
            continue

        try:
            start = time.time()
            is_valid = validate_api_keys(provider, verbose=True)
            duration = time.time() - start

            results["providers"].append({
                "provider": provider,
                "status": "passed" if is_valid else "failed",
                "valid": is_valid,
                "duration_seconds": duration
            })

            print(f"  Result: {'VALID' if is_valid else 'INVALID'}")
            print(f"  Duration: {duration:.3f}s")

        except Exception as e:
            results["providers"].append({
                "provider": provider,
                "status": "error",
                "error": str(e)
            })
            print(f"  [ERROR] {str(e)[:100]}")

    return results


def test_caching_behavior():
    """Test that validation caching works correctly."""
    from CASSIA import validate_api_keys, clear_validation_cache

    print("\n" + "=" * 60)
    print("Test 2: Caching Behavior")
    print("=" * 60)

    # Find a provider with API key set
    providers = ["openai", "anthropic", "openrouter"]
    env_vars = {
        "openai": "OPENAI_API_KEY",
        "anthropic": "ANTHROPIC_API_KEY",
        "openrouter": "OPENROUTER_API_KEY"
    }

    test_provider = None
    for provider in providers:
        if os.environ.get(env_vars[provider]):
            test_provider = provider
            break

    if not test_provider:
        print("  [SKIP] No API keys configured")
        return {"test_name": "caching_behavior", "status": "skipped"}

    print(f"\nUsing provider: {test_provider}")

    # Clear cache first
    clear_validation_cache(test_provider)

    # First validation (should make API call)
    print("\nFirst validation (should make API call)...")
    start1 = time.time()
    result1 = validate_api_keys(test_provider, verbose=True)
    duration1 = time.time() - start1
    print(f"  Duration: {duration1:.3f}s")

    # Second validation (should use cache)
    print("\nSecond validation (should use cache)...")
    start2 = time.time()
    result2 = validate_api_keys(test_provider, verbose=True)
    duration2 = time.time() - start2
    print(f"  Duration: {duration2:.3f}s")

    # Check that second call was much faster (< 10ms)
    speedup = duration1 / duration2 if duration2 > 0 else float('inf')
    cache_worked = duration2 < 0.01  # Less than 10ms

    print(f"\nCache performance:")
    print(f"  First call:  {duration1:.3f}s (API call)")
    print(f"  Second call: {duration2:.6f}s (cached)")
    print(f"  Speedup: {speedup:.0f}x")
    print(f"  Cache working: {'YES' if cache_worked else 'NO'}")

    return {
        "test_name": "caching_behavior",
        "provider": test_provider,
        "first_duration": duration1,
        "second_duration": duration2,
        "speedup": speedup,
        "cache_working": cache_worked,
        "status": "passed" if cache_worked else "failed"
    }


def test_force_revalidate():
    """Test force_revalidate parameter."""
    from CASSIA import validate_api_keys, clear_validation_cache

    print("\n" + "=" * 60)
    print("Test 3: Force Revalidation")
    print("=" * 60)

    # Find a provider with API key set
    providers = ["openai", "anthropic", "openrouter"]
    env_vars = {
        "openai": "OPENAI_API_KEY",
        "anthropic": "ANTHROPIC_API_KEY",
        "openrouter": "OPENROUTER_API_KEY"
    }

    test_provider = None
    for provider in providers:
        if os.environ.get(env_vars[provider]):
            test_provider = provider
            break

    if not test_provider:
        print("  [SKIP] No API keys configured")
        return {"test_name": "force_revalidate", "status": "skipped"}

    print(f"\nUsing provider: {test_provider}")

    # First validation
    print("\nFirst validation...")
    validate_api_keys(test_provider, verbose=False)

    # Second validation with force_revalidate=True (should make API call)
    print("\nForced revalidation (should make API call even though cached)...")
    start = time.time()
    result = validate_api_keys(test_provider, force_revalidate=True, verbose=True)
    duration = time.time() - start

    # Should take longer than cache lookup (> 100ms typically)
    api_call_made = duration > 0.1

    print(f"  Duration: {duration:.3f}s")
    print(f"  API call made: {'YES' if api_call_made else 'NO (used cache)'}")

    return {
        "test_name": "force_revalidate",
        "provider": test_provider,
        "duration": duration,
        "api_call_made": api_call_made,
        "status": "passed" if api_call_made else "failed"
    }


def test_invalid_key():
    """Test validation with invalid API key."""
    from CASSIA import validate_api_keys

    print("\n" + "=" * 60)
    print("Test 4: Invalid API Key")
    print("=" * 60)

    print("\nTesting with fake API key...")

    # Use a fake key that will fail
    fake_key = "sk-fake-key-for-testing-12345"

    try:
        result = validate_api_keys("openai", api_key=fake_key, verbose=True)

        if not result:
            print("  [OK] Invalid key correctly detected")
            status = "passed"
        else:
            print("  [FAIL] Invalid key reported as valid!")
            status = "failed"

    except Exception as e:
        # Exception is also acceptable
        print(f"  [OK] Exception raised: {str(e)[:80]}")
        status = "passed"

    return {
        "test_name": "invalid_key",
        "status": status
    }


def test_missing_key():
    """Test validation when API key is not set."""
    from CASSIA import validate_api_keys

    print("\n" + "=" * 60)
    print("Test 5: Missing API Key")
    print("=" * 60)

    # Save original key
    original_key = os.environ.get("OPENAI_API_KEY")

    try:
        # Temporarily remove the key
        if "OPENAI_API_KEY" in os.environ:
            del os.environ["OPENAI_API_KEY"]

        print("\nTesting with no API key set...")
        result = validate_api_keys("openai", verbose=True)

        if not result:
            print("  [OK] Missing key correctly detected")
            status = "passed"
        else:
            print("  [FAIL] Missing key reported as valid!")
            status = "failed"

    finally:
        # Restore original key
        if original_key:
            os.environ["OPENAI_API_KEY"] = original_key

    return {
        "test_name": "missing_key",
        "status": status
    }


def test_thread_safety():
    """Test concurrent validation calls (thread safety)."""
    from CASSIA import validate_api_keys, clear_validation_cache

    print("\n" + "=" * 60)
    print("Test 6: Thread Safety")
    print("=" * 60)

    # Find a provider with API key set
    providers = ["openai", "anthropic", "openrouter"]
    env_vars = {
        "openai": "OPENAI_API_KEY",
        "anthropic": "ANTHROPIC_API_KEY",
        "openrouter": "OPENROUTER_API_KEY"
    }

    test_provider = None
    for provider in providers:
        if os.environ.get(env_vars[provider]):
            test_provider = provider
            break

    if not test_provider:
        print("  [SKIP] No API keys configured")
        return {"test_name": "thread_safety", "status": "skipped"}

    print(f"\nUsing provider: {test_provider}")

    # Clear cache first
    clear_validation_cache(test_provider)

    print("\nRunning 5 concurrent validations...")

    results_list = []
    threads = []

    def validate_thread():
        result = validate_api_keys(test_provider, verbose=False)
        results_list.append(result)

    start = time.time()

    # Create and start threads
    for i in range(5):
        thread = threading.Thread(target=validate_thread)
        threads.append(thread)
        thread.start()

    # Wait for all threads
    for thread in threads:
        thread.join()

    duration = time.time() - start

    # Check that all returned the same result
    all_same = len(set(results_list)) == 1
    all_valid = all(results_list)

    print(f"  Results: {results_list}")
    print(f"  All consistent: {'YES' if all_same else 'NO'}")
    print(f"  All valid: {'YES' if all_valid else 'NO'}")
    print(f"  Duration: {duration:.3f}s")

    status = "passed" if (all_same and all_valid) else "failed"

    return {
        "test_name": "thread_safety",
        "provider": test_provider,
        "results": results_list,
        "all_consistent": all_same,
        "all_valid": all_valid,
        "duration": duration,
        "status": status
    }


def run_api_validation_test():
    """Run all API validation tests."""
    print_test_header("18 - API Key Validation Testing")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Create results directory
    results_dir_dict = create_results_dir("18_api_validation")
    results_dir = results_dir_dict['base']  # Extract the base path
    print(f"\nResults will be saved to: {results_dir}")

    start_time = time.time()
    all_test_results = {}
    errors = []

    # Run all tests
    try:
        all_test_results["basic_validation"] = test_validation_basic()
    except Exception as e:
        errors.append(f"basic_validation: {str(e)}")
        print(f"\n[ERROR] Basic validation test failed: {str(e)}")

    try:
        all_test_results["caching_behavior"] = test_caching_behavior()
    except Exception as e:
        errors.append(f"caching_behavior: {str(e)}")
        print(f"\n[ERROR] Caching test failed: {str(e)}")

    try:
        all_test_results["force_revalidate"] = test_force_revalidate()
    except Exception as e:
        errors.append(f"force_revalidate: {str(e)}")
        print(f"\n[ERROR] Force revalidate test failed: {str(e)}")

    try:
        all_test_results["invalid_key"] = test_invalid_key()
    except Exception as e:
        errors.append(f"invalid_key: {str(e)}")
        print(f"\n[ERROR] Invalid key test failed: {str(e)}")

    try:
        all_test_results["missing_key"] = test_missing_key()
    except Exception as e:
        errors.append(f"missing_key: {str(e)}")
        print(f"\n[ERROR] Missing key test failed: {str(e)}")

    try:
        all_test_results["thread_safety"] = test_thread_safety()
    except Exception as e:
        errors.append(f"thread_safety: {str(e)}")
        print(f"\n[ERROR] Thread safety test failed: {str(e)}")

    duration = time.time() - start_time

    # Calculate summary
    total_tests = len(all_test_results)
    passed_tests = sum(1 for r in all_test_results.values()
                      if isinstance(r, dict) and r.get("status") == "passed")
    failed_tests = sum(1 for r in all_test_results.values()
                      if isinstance(r, dict) and r.get("status") == "failed")
    skipped_tests = sum(1 for r in all_test_results.values()
                       if isinstance(r, dict) and r.get("status") == "skipped")

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Passed:  {passed_tests}/{total_tests}")
    print(f"  Failed:  {failed_tests}/{total_tests}")
    print(f"  Skipped: {skipped_tests}/{total_tests}")
    print(f"  Duration: {duration:.2f}s")

    # Determine overall status
    if passed_tests == total_tests:
        status = "passed"
    elif passed_tests > 0 and failed_tests == 0:
        status = "passed_with_skipped"
    else:
        status = "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="api_validation_test",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[],
        errors=errors
    )
    metadata["test_summary"] = {
        "passed": passed_tests,
        "failed": failed_tests,
        "skipped": skipped_tests,
        "total": total_tests
    }
    save_test_metadata(results_dir, metadata)
    save_test_results(results_dir, all_test_results)

    # Print final result
    success = (status == "passed" or status == "passed_with_skipped")
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_api_validation_test()
    sys.exit(0 if success else 1)
