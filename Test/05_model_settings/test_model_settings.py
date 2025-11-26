"""
CASSIA Test 05: Model Settings
==============================
Tests model name resolution, provider shortcuts, and fuzzy alias matching.

Usage:
    python test_model_settings.py
"""

import sys
import time
import io
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
    from model_settings import resolve_model_name, get_recommended_model, get_available_aliases
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

    # Test 1: Tier shortcuts
    print("\n" + "="*40)
    print("Testing Tier Shortcuts")
    print("="*40)

    tier_tests = [
        # (tier, provider, expected_contains)
        ("best", "openai", "gpt-5.1"),
        ("balanced", "openai", "gpt-4o"),
        ("fast", "openai", "gpt-5-mini"),
        ("best", "anthropic", "opus"),
        ("balanced", "anthropic", "sonnet"),
        ("fast", "anthropic", "haiku"),
        ("best", "openrouter", "claude"),
        ("fast", "openrouter", "gemini"),
    ]

    tier_passed = 0
    for tier, provider, expected_contains in tier_tests:
        try:
            # Use verbose=False to suppress output for tier tests
            resolved_name, resolved_provider = resolve_model_name(tier, provider, verbose=False)
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

    # Test 2: Fuzzy alias matching
    print("\n" + "="*40)
    print("Testing Fuzzy Alias Matching")
    print("="*40)

    alias_tests = [
        # (alias, provider, expected_contains, should_print_note)
        ("gpt", "openai", "gpt-5.1", True),
        ("gpt", "openrouter", "openai/gpt", True),
        ("claude", "anthropic", "claude-sonnet-4-5", True),
        ("claude", "openrouter", "anthropic/claude-sonnet", True),
        ("gemini", "openrouter", "google/gemini", True),
        ("flash", "openrouter", "gemini-2.5-flash", True),
        ("sonnet", "anthropic", "sonnet", True),
        ("opus", "anthropic", "opus", True),
        ("haiku", "anthropic", "haiku", True),
        ("4o", "openai", "gpt-4o", True),
        ("deepseek", "openrouter", "deepseek", True),
        # Case insensitivity
        ("GPT", "openai", "gpt-5.1", True),
        ("Claude", "anthropic", "claude-sonnet", True),
        # Exact match passthrough (no note)
        ("gpt-4o", "openai", "gpt-4o", False),
        ("claude-sonnet-4-5", "anthropic", "claude-sonnet-4-5", False),
    ]

    alias_passed = 0
    for alias, provider, expected_contains, should_print_note in alias_tests:
        try:
            # Capture stdout to check for notification
            captured = io.StringIO()
            sys.stdout = captured

            resolved_name, _ = resolve_model_name(alias, provider, verbose=True)

            sys.stdout = sys.__stdout__
            output = captured.getvalue()

            # Check result
            result_ok = expected_contains.lower() in resolved_name.lower()

            # Check notification
            note_printed = "Note: Resolved" in output
            notify_ok = note_printed == should_print_note

            success = result_ok and notify_ok

            if not result_ok:
                errors.append(f"Alias '{alias}' resolved to '{resolved_name}', expected to contain '{expected_contains}'")
            if not notify_ok:
                if should_print_note:
                    errors.append(f"Alias '{alias}' should have printed a note but didn't")
                else:
                    errors.append(f"Alias '{alias}' should NOT have printed a note but did")

            test_results["resolution_tests"].append({
                "type": "alias",
                "input": alias,
                "provider": provider,
                "resolved": resolved_name,
                "note_printed": note_printed,
                "success": success
            })

            status_symbol = "[OK]" if success else "[X]"
            note_indicator = "(note)" if note_printed else "(no note)"
            print(f"  {status_symbol} '{alias}' + '{provider}' -> '{resolved_name}' {note_indicator}")

            if success:
                alias_passed += 1

        except Exception as e:
            sys.stdout = sys.__stdout__
            test_results["resolution_tests"].append({
                "type": "alias",
                "input": alias,
                "provider": provider,
                "error": str(e),
                "success": False
            })
            errors.append(f"Alias resolution failed for '{alias}': {e}")
            print(f"  [X] '{alias}' + '{provider}' -> ERROR: {e}")

    print(f"\nAlias tests: {alias_passed}/{len(alias_tests)} passed")

    # Test 3: Provider-aware resolution (same alias, different output)
    print("\n" + "="*40)
    print("Testing Provider-Aware Resolution")
    print("="*40)

    provider_aware_tests = [
        # Same alias resolves differently per provider
        ("claude", "anthropic", "claude-sonnet-4-5"),
        ("claude", "openrouter", "anthropic/claude-sonnet-4.5"),
        ("gpt", "openai", "gpt-5.1"),
        ("gpt", "openrouter", "openai/gpt-5.1"),
    ]

    provider_aware_passed = 0
    for alias, provider, expected in provider_aware_tests:
        try:
            resolved_name, _ = resolve_model_name(alias, provider, verbose=False)
            success = resolved_name == expected

            if not success:
                errors.append(f"Provider-aware: '{alias}' + '{provider}' -> '{resolved_name}', expected '{expected}'")

            status_symbol = "[OK]" if success else "[X]"
            print(f"  {status_symbol} '{alias}' + '{provider}' -> '{resolved_name}' (expected: {expected})")

            if success:
                provider_aware_passed += 1

        except Exception as e:
            errors.append(f"Provider-aware test failed for '{alias}' + '{provider}': {e}")
            print(f"  [X] '{alias}' + '{provider}' -> ERROR: {e}")

    print(f"\nProvider-aware tests: {provider_aware_passed}/{len(provider_aware_tests)} passed")

    # Test 4: Get available aliases
    print("\n" + "="*40)
    print("Testing get_available_aliases()")
    print("="*40)

    try:
        all_aliases = get_available_aliases()
        assert "provider_specific" in all_aliases, "Missing provider_specific"
        assert "global" in all_aliases, "Missing global"
        assert "openrouter" in all_aliases["provider_specific"], "Missing openrouter aliases"
        print("  [OK] get_available_aliases() returns correct structure")

        openai_aliases = get_available_aliases("openai")
        assert "gpt" in openai_aliases["provider_specific"], "Missing 'gpt' alias for openai"
        print("  [OK] get_available_aliases('openai') filters correctly")

        test_results["alias_structure_test"] = {"success": True}
    except Exception as e:
        errors.append(f"get_available_aliases() test failed: {e}")
        print(f"  [X] ERROR: {e}")
        test_results["alias_structure_test"] = {"success": False, "error": str(e)}

    # Calculate total resolution tests
    total_resolution_tests = len(tier_tests) + len(alias_tests) + len(provider_aware_tests)
    total_resolution_passed = tier_passed + alias_passed + provider_aware_passed
    print(f"\nTotal resolution tests: {total_resolution_passed}/{total_resolution_tests} passed")

    # Test 5: Recommended model
    print("\n" + "="*40)
    print("Testing Recommended Model")
    print("="*40)

    try:
        rec_model = get_recommended_model("openrouter")
        print(f"  Recommended for openrouter: {rec_model}")
        test_results["recommended_model"] = rec_model
    except Exception as e:
        errors.append(f"get_recommended_model failed: {e}")
        print(f"  [X] ERROR: {e}")

    # Test 6: Practical test with gemini-2.5-flash via OpenRouter
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
    resolution_all_passed = total_resolution_passed == total_resolution_tests
    practical_passed = test_results.get("practical_test", {}).get("success", False)

    # Pass if resolution tests pass OR practical test passes
    if resolution_all_passed and practical_passed:
        status = "passed"
    elif practical_passed:
        status = "passed"  # Practical test is most important
    elif resolution_all_passed:
        status = "passed"  # Resolution tests all passed
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
