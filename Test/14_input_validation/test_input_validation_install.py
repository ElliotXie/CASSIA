"""
CASSIA Test 14: Input Validation (PIP INSTALL MODE)
====================================================
Tests input validation for runCASSIA and runCASSIA_batch functions
using pip-installed CASSIA.

Usage:
    python test_input_validation_install.py
"""

import sys
import warnings
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import (
    print_test_header,
    print_test_result,
    verify_cassia_pip_install
)
from result_manager import (
    create_results_dir,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA

# Import validation functions and exceptions from pip-installed package
from CASSIA.core.validation import (
    validate_marker_list,
    validate_temperature,
    validate_tissue,
    validate_species,
    validate_provider,
    validate_model,
    validate_positive_int,
    validate_runCASSIA_inputs
)
from CASSIA.core.exceptions import (
    CASSIAValidationError,
    MarkerValidationError,
    TemperatureValidationError,
    TissueSpeciesValidationError,
    ProviderValidationError,
    ModelValidationError,
    BatchParameterValidationError
)


def test_marker_list_validation():
    """Test marker_list validation."""
    print("\n" + "="*50)
    print("Testing Marker List Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: Valid list of strings
    try:
        result = validate_marker_list(["CD4", "CD8", "FOXP3"])
        assert result == ["CD4", "CD8", "FOXP3"]
        print("  [OK] Valid list of strings")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] Valid list of strings: {e}")
        tests_failed += 1

    # Test 2: None should raise error
    try:
        validate_marker_list(None)
        print("  [X] None raises error: Expected error")
        tests_failed += 1
    except MarkerValidationError:
        print("  [OK] None raises error")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] None raises error: Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 3: Empty list should raise error
    try:
        validate_marker_list([])
        print("  [X] Empty list raises error: Expected error")
        tests_failed += 1
    except MarkerValidationError:
        print("  [OK] Empty list raises error")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] Empty list raises error: Wrong exception: {type(e)}")
        tests_failed += 1

    return tests_passed, tests_failed


def test_temperature_validation():
    """Test temperature validation."""
    print("\n" + "="*50)
    print("Testing Temperature Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: Valid temperature
    try:
        result = validate_temperature(0)
        assert result == 0.0
        print("  [OK] Valid temperature (0)")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] Valid temperature (0): {e}")
        tests_failed += 1

    # Test 2: Negative temperature should raise error
    try:
        validate_temperature(-0.5)
        print("  [X] Negative raises error: Expected error")
        tests_failed += 1
    except TemperatureValidationError:
        print("  [OK] Negative raises error")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] Negative raises error: Wrong exception: {type(e)}")
        tests_failed += 1

    return tests_passed, tests_failed


def test_provider_validation():
    """Test provider validation."""
    print("\n" + "="*50)
    print("Testing Provider Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: Valid providers
    for provider in ["openai", "anthropic", "openrouter"]:
        try:
            result = validate_provider(provider)
            assert result == provider
            print(f"  [OK] Valid provider '{provider}'")
            tests_passed += 1
        except Exception as e:
            print(f"  [X] Valid provider '{provider}': {e}")
            tests_failed += 1

    # Test 2: Unknown provider raises error
    try:
        validate_provider("invalid_provider")
        print("  [X] Unknown provider raises error: Expected error")
        tests_failed += 1
    except ProviderValidationError:
        print("  [OK] Unknown provider raises error")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] Unknown provider raises error: Wrong exception: {type(e)}")
        tests_failed += 1

    return tests_passed, tests_failed


def test_model_validation():
    """Test model validation."""
    print("\n" + "="*50)
    print("Testing Model Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: Valid model names
    try:
        result = validate_model("gpt-4")
        assert result == "gpt-4"
        print("  [OK] Valid model 'gpt-4'")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] Valid model 'gpt-4': {e}")
        tests_failed += 1

    # Test 2: None raises error
    try:
        validate_model(None)
        print("  [X] None model raises error: Expected error")
        tests_failed += 1
    except ModelValidationError:
        print("  [OK] None model raises error")
        tests_passed += 1
    except Exception as e:
        print(f"  [X] None model raises error: Wrong exception: {type(e)}")
        tests_failed += 1

    return tests_passed, tests_failed


def run_all_tests():
    """Run all validation tests."""
    print_test_header("14 - Input Validation (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    total_passed = 0
    total_failed = 0

    # Run all test groups
    passed, failed = test_marker_list_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_temperature_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_provider_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_model_validation()
    total_passed += passed
    total_failed += failed

    # Print summary
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    print(f"Tests passed: {total_passed}")
    print(f"Tests failed: {total_failed}")
    print(f"Total tests:  {total_passed + total_failed}")

    if total_failed == 0:
        print("\n*** ALL TESTS PASSED ***")
        return 0
    else:
        print(f"\n*** {total_failed} TEST(S) FAILED ***")
        return 1


def main():
    """Main entry point with logging."""
    # Create results directory for logging
    results_dir = create_results_dir("14_input_validation")

    # Setup logging to capture all console output
    logging_context = setup_logging(results_dir['logs'])

    try:
        exit_code = run_all_tests()
    finally:
        cleanup_logging(logging_context)

    return exit_code


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
