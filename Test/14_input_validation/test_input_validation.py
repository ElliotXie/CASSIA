"""
CASSIA Test 14: Input Validation
================================
Tests input validation for runCASSIA and runCASSIA_batch functions.

This test suite verifies that:
1. Invalid inputs are rejected with clear error messages
2. Valid inputs are accepted and processed correctly
3. Edge cases are handled appropriately (warnings vs errors)

Usage:
    python test_input_validation.py
"""

import sys
import warnings
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import (
    setup_cassia_imports,
    print_test_header,
    print_test_result
)

# Setup CASSIA imports
setup_cassia_imports()

# Import validation functions and exceptions
from validation import (
    validate_marker_list,
    validate_temperature,
    validate_tissue,
    validate_species,
    validate_provider,
    validate_model,
    validate_marker_dataframe,
    validate_positive_int,
    validate_ranking_method,
    validate_runCASSIA_inputs
)
from exceptions import (
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
        assert result == ["CD4", "CD8", "FOXP3"], f"Expected normalized list, got {result}"
        print_test_result("Valid list of strings", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Valid list of strings", False, str(e))
        tests_failed += 1

    # Test 2: Valid comma-separated string
    try:
        result = validate_marker_list("CD4, CD8, FOXP3")
        assert result == ["CD4", "CD8", "FOXP3"], f"Expected ['CD4', 'CD8', 'FOXP3'], got {result}"
        print_test_result("Valid comma-separated string", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Valid comma-separated string", False, str(e))
        tests_failed += 1

    # Test 3: List with one comma-separated string
    try:
        result = validate_marker_list(["CD4, CD8, FOXP3"])
        assert result == ["CD4", "CD8", "FOXP3"], f"Expected ['CD4', 'CD8', 'FOXP3'], got {result}"
        print_test_result("List with comma-separated string", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("List with comma-separated string", False, str(e))
        tests_failed += 1

    # Test 4: None should raise error
    try:
        validate_marker_list(None)
        print_test_result("None raises error", False, "Expected MarkerValidationError")
        tests_failed += 1
    except MarkerValidationError as e:
        assert "cannot be None" in str(e), f"Wrong error message: {e}"
        print_test_result("None raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("None raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 5: Empty list should raise error
    try:
        validate_marker_list([])
        print_test_result("Empty list raises error", False, "Expected MarkerValidationError")
        tests_failed += 1
    except MarkerValidationError as e:
        assert "empty" in str(e).lower(), f"Wrong error message: {e}"
        print_test_result("Empty list raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Empty list raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 6: Non-string in list should raise error
    try:
        validate_marker_list(["CD4", 123, "CD8"])
        print_test_result("Non-string raises error", False, "Expected MarkerValidationError")
        tests_failed += 1
    except MarkerValidationError as e:
        assert "not a string" in str(e), f"Wrong error message: {e}"
        print_test_result("Non-string raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Non-string raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 7: Wrong type should raise error
    try:
        validate_marker_list(12345)
        print_test_result("Wrong type raises error", False, "Expected MarkerValidationError")
        tests_failed += 1
    except MarkerValidationError as e:
        assert "Invalid type" in str(e), f"Wrong error message: {e}"
        print_test_result("Wrong type raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Wrong type raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 8: Too many markers should raise error
    try:
        many_markers = [f"GENE{i}" for i in range(600)]
        validate_marker_list(many_markers)
        print_test_result("Too many markers raises error", False, "Expected MarkerValidationError")
        tests_failed += 1
    except MarkerValidationError as e:
        assert "Too many markers" in str(e), f"Wrong error message: {e}"
        print_test_result("Too many markers raises error (>500)", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Too many markers raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 9: Ensembl IDs should raise error
    try:
        validate_marker_list(["ENSG00000141510", "ENSG00000134644", "ENSG00000087086"])
        print_test_result("Ensembl IDs raise error", False, "Expected MarkerValidationError")
        tests_failed += 1
    except MarkerValidationError as e:
        assert "Ensembl" in str(e), f"Wrong error message: {e}"
        print_test_result("Ensembl IDs raise error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Ensembl IDs raise error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 10: Entrez IDs (>30% numeric) should raise error
    try:
        # 100% numeric - should error
        validate_marker_list(["7157", "675", "472"])
        print_test_result("Entrez IDs (>30% numeric) raise error", False, "Expected MarkerValidationError")
        tests_failed += 1
    except MarkerValidationError as e:
        assert "numeric" in str(e).lower() or "Entrez" in str(e), f"Wrong error message: {e}"
        print_test_result("Entrez IDs (>30% numeric) raise error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Entrez IDs (>30% numeric) raise error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 10b: Few numeric markers (<30%) should warn but continue
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # 1 numeric out of 10 = 10% - should warn but continue
            result = validate_marker_list(["CD4", "CD8", "FOXP3", "CD3", "7157", "CD19", "CD20", "CD56", "CD14", "HLA-DR"])
            assert len(result) == 10
            # Should have warning about numeric marker
            numeric_warnings = [x for x in w if "numeric" in str(x.message).lower()]
            assert len(numeric_warnings) >= 1, "Expected warning about numeric marker"
        print_test_result("Few numeric markers (<30%) warn but continue", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Few numeric markers (<30%) warn but continue", False, str(e))
        tests_failed += 1

    # Test 11: Few markers should warn but continue
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_marker_list(["CD4", "CD8", "CD3"])
            # Should succeed with warning
            assert len(result) == 3
            assert len(w) == 1
            assert "Recommended minimum" in str(w[0].message)
        print_test_result("Few markers warns but continues (<10)", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Few markers warns but continues", False, str(e))
        tests_failed += 1

    return tests_passed, tests_failed


def test_temperature_validation():
    """Test temperature validation."""
    print("\n" + "="*50)
    print("Testing Temperature Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: Valid temperature (0)
    try:
        result = validate_temperature(0)
        assert result == 0.0
        print_test_result("Valid temperature (0)", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Valid temperature (0)", False, str(e))
        tests_failed += 1

    # Test 2: Valid temperature (1.5)
    try:
        result = validate_temperature(1.5)
        assert result == 1.5
        print_test_result("Valid temperature (1.5)", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Valid temperature (1.5)", False, str(e))
        tests_failed += 1

    # Test 3: Integer coerced to float
    try:
        result = validate_temperature(1)
        assert result == 1.0
        print_test_result("Integer coerced to float", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Integer coerced to float", False, str(e))
        tests_failed += 1

    # Test 4: Negative temperature should raise error
    try:
        validate_temperature(-0.5)
        print_test_result("Negative raises error", False, "Expected TemperatureValidationError")
        tests_failed += 1
    except TemperatureValidationError as e:
        assert ">= 0" in str(e), f"Wrong error message: {e}"
        print_test_result("Negative raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Negative raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 5: None should raise error
    try:
        validate_temperature(None)
        print_test_result("None raises error", False, "Expected TemperatureValidationError")
        tests_failed += 1
    except TemperatureValidationError as e:
        assert "cannot be None" in str(e), f"Wrong error message: {e}"
        print_test_result("None raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("None raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 6: String should raise error
    try:
        validate_temperature("high")
        print_test_result("String raises error", False, "Expected TemperatureValidationError")
        tests_failed += 1
    except TemperatureValidationError as e:
        assert "must be a number" in str(e), f"Wrong error message: {e}"
        print_test_result("String raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("String raises error", False, f"Wrong exception type: {type(e)}")
        tests_failed += 1

    # Test 7: High temperature (5.0) should be accepted (no upper bound per user request)
    try:
        result = validate_temperature(5.0)
        assert result == 5.0
        print_test_result("High temperature accepted (no upper bound)", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("High temperature accepted", False, str(e))
        tests_failed += 1

    return tests_passed, tests_failed


def test_tissue_species_validation():
    """Test tissue and species validation."""
    print("\n" + "="*50)
    print("Testing Tissue/Species Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: Valid tissue
    try:
        result = validate_tissue("lung")
        assert result == "lung"
        print_test_result("Valid tissue", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Valid tissue", False, str(e))
        tests_failed += 1

    # Test 2: Special tissue values
    try:
        result = validate_tissue("none")
        assert result == "none"
        print_test_result("Special tissue 'none'", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Special tissue 'none'", False, str(e))
        tests_failed += 1

    # Test 3: Empty tissue warns but continues
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_tissue("")
            assert result == ""
            assert len(w) == 1
            assert "empty" in str(w[0].message).lower()
        print_test_result("Empty tissue warns but continues", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Empty tissue warns but continues", False, str(e))
        tests_failed += 1

    # Test 4: None tissue warns but continues
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_tissue(None)
            assert result == ""
            assert len(w) == 1
        print_test_result("None tissue warns but continues", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("None tissue warns but continues", False, str(e))
        tests_failed += 1

    # Test 5: Non-string tissue raises error
    try:
        validate_tissue(123)
        print_test_result("Non-string tissue raises error", False, "Expected error")
        tests_failed += 1
    except TissueSpeciesValidationError:
        print_test_result("Non-string tissue raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Non-string tissue raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 6: Valid species
    try:
        result = validate_species("human")
        assert result == "human"
        print_test_result("Valid species", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Valid species", False, str(e))
        tests_failed += 1

    # Test 7: Empty species warns but continues
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_species("")
            assert result == ""
            assert len(w) == 1
        print_test_result("Empty species warns but continues", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Empty species warns but continues", False, str(e))
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
            print_test_result(f"Valid provider '{provider}'", True)
            tests_passed += 1
        except Exception as e:
            print_test_result(f"Valid provider '{provider}'", False, str(e))
            tests_failed += 1

    # Test 2: Case insensitive
    try:
        result = validate_provider("OpenAI")
        assert result == "openai"
        print_test_result("Case insensitive provider", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Case insensitive provider", False, str(e))
        tests_failed += 1

    # Test 3: HTTP URL accepted
    try:
        url = "http://localhost:8000/v1"
        result = validate_provider(url)
        assert result == url
        print_test_result("HTTP URL accepted", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("HTTP URL accepted", False, str(e))
        tests_failed += 1

    # Test 4: HTTPS URL accepted
    try:
        url = "https://api.example.com/v1"
        result = validate_provider(url)
        assert result == url
        print_test_result("HTTPS URL accepted", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("HTTPS URL accepted", False, str(e))
        tests_failed += 1

    # Test 5: Unknown provider raises error
    try:
        validate_provider("gemini")
        print_test_result("Unknown provider raises error", False, "Expected error")
        tests_failed += 1
    except ProviderValidationError as e:
        assert "Unknown provider" in str(e)
        print_test_result("Unknown provider raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Unknown provider raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 6: None raises error
    try:
        validate_provider(None)
        print_test_result("None provider raises error", False, "Expected error")
        tests_failed += 1
    except ProviderValidationError:
        print_test_result("None provider raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("None provider raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 7: Empty string raises error
    try:
        validate_provider("")
        print_test_result("Empty provider raises error", False, "Expected error")
        tests_failed += 1
    except ProviderValidationError:
        print_test_result("Empty provider raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Empty provider raises error", False, f"Wrong exception: {type(e)}")
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
    for model in ["gpt-4", "claude-3", "gemini-pro"]:
        try:
            result = validate_model(model)
            assert result == model
            print_test_result(f"Valid model '{model}'", True)
            tests_passed += 1
        except Exception as e:
            print_test_result(f"Valid model '{model}'", False, str(e))
            tests_failed += 1

    # Test 2: None raises error
    try:
        validate_model(None)
        print_test_result("None model raises error", False, "Expected error")
        tests_failed += 1
    except ModelValidationError:
        print_test_result("None model raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("None model raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 3: Empty string raises error
    try:
        validate_model("")
        print_test_result("Empty model raises error", False, "Expected error")
        tests_failed += 1
    except ModelValidationError:
        print_test_result("Empty model raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Empty model raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 4: Non-string raises error
    try:
        validate_model(123)
        print_test_result("Non-string model raises error", False, "Expected error")
        tests_failed += 1
    except ModelValidationError:
        print_test_result("Non-string model raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Non-string model raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    return tests_passed, tests_failed


def test_batch_parameters():
    """Test batch-specific parameter validation."""
    print("\n" + "="*50)
    print("Testing Batch Parameter Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: Valid positive integers
    try:
        assert validate_positive_int(10, "n_genes") == 10
        assert validate_positive_int(5, "max_workers") == 5
        print_test_result("Valid positive integers", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Valid positive integers", False, str(e))
        tests_failed += 1

    # Test 2: Zero allowed with flag
    try:
        result = validate_positive_int(0, "max_retries", allow_zero=True)
        assert result == 0
        print_test_result("Zero allowed with allow_zero=True", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Zero allowed with allow_zero=True", False, str(e))
        tests_failed += 1

    # Test 3: Zero rejected without flag
    try:
        validate_positive_int(0, "n_genes", allow_zero=False)
        print_test_result("Zero rejected without allow_zero", False, "Expected error")
        tests_failed += 1
    except BatchParameterValidationError:
        print_test_result("Zero rejected without allow_zero", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Zero rejected without allow_zero", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 4: Negative raises error
    try:
        validate_positive_int(-1, "max_workers")
        print_test_result("Negative raises error", False, "Expected error")
        tests_failed += 1
    except BatchParameterValidationError:
        print_test_result("Negative raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Negative raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 5: Valid ranking methods
    for method in ["avg_log2FC", "p_val_adj", "pct_diff", "Score"]:
        try:
            result = validate_ranking_method(method)
            assert result == method
            print_test_result(f"Valid ranking method '{method}'", True)
            tests_passed += 1
        except Exception as e:
            print_test_result(f"Valid ranking method '{method}'", False, str(e))
            tests_failed += 1

    # Test 6: Invalid ranking method raises error
    try:
        validate_ranking_method("invalid_method")
        print_test_result("Invalid ranking method raises error", False, "Expected error")
        tests_failed += 1
    except BatchParameterValidationError:
        print_test_result("Invalid ranking method raises error", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Invalid ranking method raises error", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    return tests_passed, tests_failed


def test_integrated_validation():
    """Test the integrated validation function."""
    print("\n" + "="*50)
    print("Testing Integrated runCASSIA Validation")
    print("="*50)

    tests_passed = 0
    tests_failed = 0

    # Test 1: All valid inputs
    try:
        markers = ["CD4", "CD8", "FOXP3", "CD3", "CD19", "CD20", "CD56", "CD14", "CD68", "HLA-DR"]
        result = validate_runCASSIA_inputs(
            model="gpt-4",
            temperature=0,
            marker_list=markers,
            tissue="lung",
            species="human",
            provider="openai",
            additional_info="Test analysis"
        )
        assert 'marker_list' in result
        assert 'temperature' in result
        assert len(result['marker_list']) == 10
        print_test_result("All valid inputs accepted", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("All valid inputs accepted", False, str(e))
        tests_failed += 1

    # Test 2: Invalid marker_list fails early
    try:
        validate_runCASSIA_inputs(
            model="gpt-4",
            temperature=0,
            marker_list=None,  # Invalid
            tissue="lung",
            species="human",
            provider="openai"
        )
        print_test_result("Invalid marker_list caught", False, "Expected error")
        tests_failed += 1
    except MarkerValidationError:
        print_test_result("Invalid marker_list caught", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Invalid marker_list caught", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 3: Invalid temperature fails
    try:
        validate_runCASSIA_inputs(
            model="gpt-4",
            temperature=-1,  # Invalid
            marker_list=["CD4", "CD8", "FOXP3", "CD3", "CD19", "CD20", "CD56", "CD14", "CD68", "HLA-DR"],
            tissue="lung",
            species="human",
            provider="openai"
        )
        print_test_result("Invalid temperature caught", False, "Expected error")
        tests_failed += 1
    except TemperatureValidationError:
        print_test_result("Invalid temperature caught", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Invalid temperature caught", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    # Test 4: Invalid provider fails
    try:
        validate_runCASSIA_inputs(
            model="gpt-4",
            temperature=0,
            marker_list=["CD4", "CD8", "FOXP3", "CD3", "CD19", "CD20", "CD56", "CD14", "CD68", "HLA-DR"],
            tissue="lung",
            species="human",
            provider="invalid"  # Invalid
        )
        print_test_result("Invalid provider caught", False, "Expected error")
        tests_failed += 1
    except ProviderValidationError:
        print_test_result("Invalid provider caught", True)
        tests_passed += 1
    except Exception as e:
        print_test_result("Invalid provider caught", False, f"Wrong exception: {type(e)}")
        tests_failed += 1

    return tests_passed, tests_failed


def run_all_tests():
    """Run all validation tests."""
    print_test_header("14 - Input Validation")

    total_passed = 0
    total_failed = 0

    # Run all test groups
    passed, failed = test_marker_list_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_temperature_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_tissue_species_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_provider_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_model_validation()
    total_passed += passed
    total_failed += failed

    passed, failed = test_batch_parameters()
    total_passed += passed
    total_failed += failed

    passed, failed = test_integrated_validation()
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


if __name__ == "__main__":
    exit_code = run_all_tests()
    sys.exit(exit_code)
