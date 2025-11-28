"""
Centralized input validation for CASSIA.

This module provides validation functions for all CASSIA public API parameters.
All validation happens early (fail-fast) with clear, actionable error messages.
"""

import os
import re
import warnings
from typing import List, Union, Any, Optional

import pandas as pd

try:
    from .exceptions import (
        CASSIAValidationError,
        MarkerValidationError,
        TemperatureValidationError,
        ProviderValidationError,
        ModelValidationError,
        TissueSpeciesValidationError,
        BatchParameterValidationError
    )
except ImportError:
    from exceptions import (
        CASSIAValidationError,
        MarkerValidationError,
        TemperatureValidationError,
        ProviderValidationError,
        ModelValidationError,
        TissueSpeciesValidationError,
        BatchParameterValidationError
    )

# Valid values for enumerated parameters
VALID_PROVIDERS = {"openai", "anthropic", "openrouter"}
VALID_MODEL_TIERS = {"best", "balanced", "fast", "recommended"}
VALID_RANKING_METHODS = {"avg_log2FC", "p_val_adj", "pct_diff", "Score"}
VALID_VALIDATOR_INVOLVEMENTS = {"v0", "v1"}

# Gene ID patterns for detection
# Ensembl gene IDs: ENSG (human), ENSMUSG (mouse), etc.
ENSEMBL_PATTERN = re.compile(r'^ENS[A-Z]*G\d{11}', re.IGNORECASE)

# Entrez/NCBI Gene IDs: pure numeric
ENTREZ_PATTERN = re.compile(r'^\d+$')

# Recommended marker count
RECOMMENDED_MIN_MARKERS = 50
WARNING_MIN_MARKERS = 10
MAX_MARKERS = 500


# =============================================================================
# Marker List Validation
# =============================================================================

def validate_marker_list(
    marker_list: Any,
    parameter_name: str = "marker_list",
    min_markers: int = 1,
    max_markers: int = MAX_MARKERS,
    check_gene_format: bool = True
) -> List[str]:
    """
    Validate and normalize marker_list input.

    Accepts:
    - List of strings: ["CD4", "CD8", "FOXP3"]
    - Single comma-separated string: "CD4, CD8, FOXP3"
    - List with one comma-separated string: ["CD4, CD8, FOXP3"]

    Args:
        marker_list: Input marker list in any accepted format
        parameter_name: Name of parameter (for error messages)
        min_markers: Minimum number of markers required (default: 1)
        max_markers: Maximum number of markers allowed (default: 500)
        check_gene_format: Whether to check for Ensembl/Entrez IDs (default: True)

    Returns:
        List[str]: Normalized list of marker gene names

    Raises:
        MarkerValidationError: If validation fails
    """
    # Check for None
    if marker_list is None:
        raise MarkerValidationError(
            "Marker list cannot be None. Please provide a list of gene markers.",
            parameter=parameter_name,
            received_value=marker_list
        )

    # Handle different input formats
    markers = _normalize_marker_input(marker_list, parameter_name)

    # Check for empty result
    if not markers:
        raise MarkerValidationError(
            "Marker list is empty after processing. "
            "Please provide at least one gene marker.",
            parameter=parameter_name,
            received_value=marker_list
        )

    # Check count bounds
    if len(markers) < min_markers:
        raise MarkerValidationError(
            f"Too few markers provided ({len(markers)}). "
            f"At least {min_markers} marker(s) required.",
            parameter=parameter_name,
            received_value=f"{len(markers)} markers"
        )

    if len(markers) > max_markers:
        raise MarkerValidationError(
            f"Too many markers provided ({len(markers)}). "
            f"Maximum {max_markers} markers allowed. "
            "If you need to analyze more markers, consider using runCASSIA_batch() "
            "with a processed marker file.",
            parameter=parameter_name,
            received_value=f"{len(markers)} markers"
        )

    # Warn if too few markers for reliable annotation
    if len(markers) < WARNING_MIN_MARKERS:
        warnings.warn(
            f"Only {len(markers)} markers provided. "
            f"Recommended minimum is {RECOMMENDED_MIN_MARKERS} markers for reliable cell type annotation. "
            "Results may be less accurate with fewer markers.",
            UserWarning
        )

    # Check for invalid gene ID formats
    if check_gene_format:
        _validate_gene_symbols(markers, parameter_name)

    return markers


def _normalize_marker_input(marker_list: Any, parameter_name: str) -> List[str]:
    """Normalize marker input to a list of strings."""

    # Case 1: Already a list
    if isinstance(marker_list, list):
        # Check if it's a list with one comma-separated string
        if (len(marker_list) == 1 and
            isinstance(marker_list[0], str) and
            ',' in marker_list[0]):
            return _split_marker_string(marker_list[0])

        # Validate each element is a string
        normalized = []
        for i, marker in enumerate(marker_list):
            if not isinstance(marker, str):
                raise MarkerValidationError(
                    f"Marker at index {i} is not a string. "
                    f"All markers must be strings (gene names).",
                    parameter=parameter_name,
                    received_value=f"type at index {i}: {type(marker).__name__}"
                )
            cleaned = marker.strip()
            if cleaned:
                normalized.append(cleaned)
        return normalized

    # Case 2: Single string (comma-separated)
    if isinstance(marker_list, str):
        return _split_marker_string(marker_list)

    # Case 3: Invalid type
    raise MarkerValidationError(
        f"Invalid type for marker list. "
        f"Expected list of strings or comma-separated string, "
        f"got {type(marker_list).__name__}.",
        parameter=parameter_name,
        received_value=type(marker_list).__name__
    )


def _split_marker_string(marker_string: str) -> List[str]:
    """Split a marker string into individual markers."""
    # Try comma+space, then comma, then space
    markers = re.split(r',\s*', marker_string)
    if len(markers) == 1:
        markers = marker_string.split(',')
    if len(markers) == 1:
        markers = marker_string.split()

    return [m.strip() for m in markers if m.strip()]


def _validate_gene_symbols(markers: List[str], parameter_name: str) -> None:
    """Validate that markers are gene symbols, not Ensembl or Entrez IDs."""
    ensembl_ids = []
    entrez_ids = []

    for marker in markers:
        if ENSEMBL_PATTERN.match(marker):
            ensembl_ids.append(marker)
        elif ENTREZ_PATTERN.match(marker) and len(marker) > 2:
            # Only flag as Entrez if it's a longer numeric string
            # (short numbers might be valid gene names like "1" in some organisms)
            entrez_ids.append(marker)

    # Error if any Ensembl IDs detected (these are very distinctive)
    if ensembl_ids:
        examples = ensembl_ids[:3]
        raise MarkerValidationError(
            f"Detected {len(ensembl_ids)} Ensembl gene ID(s) in marker list. "
            f"CASSIA requires gene symbols (e.g., 'TP53', 'CD4'), not Ensembl IDs. "
            f"Examples found: {', '.join(examples)}. "
            "Please convert Ensembl IDs to gene symbols before running CASSIA.",
            parameter=parameter_name,
            received_value=f"{len(ensembl_ids)} Ensembl IDs"
        )

    # Only error if >30% of markers are purely numeric (likely Entrez IDs)
    # This avoids false positives for a few numeric-looking markers
    if entrez_ids:
        numeric_ratio = len(entrez_ids) / len(markers)
        if numeric_ratio > 0.3:
            examples = entrez_ids[:3]
            raise MarkerValidationError(
                f"Detected {len(entrez_ids)} ({numeric_ratio:.0%}) purely numeric marker(s), "
                f"which likely are Entrez/NCBI gene IDs. "
                f"CASSIA requires gene symbols (e.g., 'TP53', 'CD4'), not numeric Entrez IDs. "
                f"Examples found: {', '.join(examples)}. "
                "Please convert Entrez IDs to gene symbols before running CASSIA.",
                parameter=parameter_name,
                received_value=f"{len(entrez_ids)} numeric IDs ({numeric_ratio:.0%})"
            )
        elif len(entrez_ids) > 0:
            # Warn but continue if only a few numeric markers
            warnings.warn(
                f"Found {len(entrez_ids)} purely numeric marker(s): {', '.join(entrez_ids[:3])}. "
                "If these are Entrez IDs, please convert to gene symbols for best results.",
                UserWarning
            )


# =============================================================================
# Temperature Validation
# =============================================================================

def validate_temperature(
    temperature: Any,
    parameter_name: str = "temperature"
) -> float:
    """
    Validate temperature parameter.

    Args:
        temperature: Temperature value (must be >= 0)
        parameter_name: Name of parameter (for error messages)

    Returns:
        float: Validated temperature value

    Raises:
        TemperatureValidationError: If validation fails
    """
    # Check for None
    if temperature is None:
        raise TemperatureValidationError(
            "Temperature cannot be None. Please provide a numeric value >= 0.",
            parameter=parameter_name,
            received_value=temperature
        )

    # Convert to float
    try:
        temp_float = float(temperature)
    except (TypeError, ValueError):
        raise TemperatureValidationError(
            f"Temperature must be a number, got {type(temperature).__name__}.",
            parameter=parameter_name,
            received_value=temperature
        )

    # Check range (only lower bound)
    if temp_float < 0.0:
        raise TemperatureValidationError(
            f"Temperature must be >= 0, got {temp_float}.",
            parameter=parameter_name,
            received_value=temp_float
        )

    return temp_float


# =============================================================================
# Tissue and Species Validation
# =============================================================================

def validate_tissue(
    tissue: Any,
    parameter_name: str = "tissue"
) -> str:
    """
    Validate tissue parameter.

    Warns but continues if tissue is empty/None.

    Args:
        tissue: Tissue type string
        parameter_name: Name of parameter (for error messages)

    Returns:
        str: Validated tissue string (may be empty)

    Raises:
        TissueSpeciesValidationError: If tissue is not a string
    """
    # Handle None
    if tissue is None:
        warnings.warn(
            "Tissue is None. Specifying a tissue type (e.g., 'lung', 'brain') "
            "may improve annotation accuracy. Use 'none' for tissue-blind analysis.",
            UserWarning
        )
        return ""

    # Check type
    if not isinstance(tissue, str):
        raise TissueSpeciesValidationError(
            f"Tissue must be a string, got {type(tissue).__name__}.",
            parameter=parameter_name,
            received_value=type(tissue).__name__
        )

    tissue_clean = tissue.strip()

    # Warn if empty
    if not tissue_clean:
        warnings.warn(
            "Tissue is empty. Specifying a tissue type (e.g., 'lung', 'brain') "
            "may improve annotation accuracy. Use 'none' for tissue-blind analysis.",
            UserWarning
        )

    return tissue_clean


def validate_species(
    species: Any,
    parameter_name: str = "species"
) -> str:
    """
    Validate species parameter.

    Warns but continues if species is empty/None.

    Args:
        species: Species string (e.g., "human", "mouse")
        parameter_name: Name of parameter (for error messages)

    Returns:
        str: Validated species string (may be empty)

    Raises:
        TissueSpeciesValidationError: If species is not a string
    """
    # Handle None
    if species is None:
        warnings.warn(
            "Species is None. Specifying a species (e.g., 'human', 'mouse') "
            "may improve annotation accuracy.",
            UserWarning
        )
        return ""

    # Check type
    if not isinstance(species, str):
        raise TissueSpeciesValidationError(
            f"Species must be a string, got {type(species).__name__}.",
            parameter=parameter_name,
            received_value=type(species).__name__
        )

    species_clean = species.strip()

    # Warn if empty
    if not species_clean:
        warnings.warn(
            "Species is empty. Specifying a species (e.g., 'human', 'mouse') "
            "may improve annotation accuracy.",
            UserWarning
        )

    return species_clean


# =============================================================================
# Model and Provider Validation
# =============================================================================

def validate_provider(
    provider: Any,
    parameter_name: str = "provider"
) -> str:
    """
    Validate provider parameter.

    Args:
        provider: Provider string ('openai', 'anthropic', 'openrouter', or HTTP URL)
        parameter_name: Name of parameter (for error messages)

    Returns:
        str: Validated provider string

    Raises:
        ProviderValidationError: If validation fails
    """
    if provider is None:
        raise ProviderValidationError(
            "Provider cannot be None. "
            f"Use one of: {', '.join(sorted(VALID_PROVIDERS))} "
            "or an HTTP URL for custom endpoints.",
            parameter=parameter_name,
            received_value=provider
        )

    if not isinstance(provider, str):
        raise ProviderValidationError(
            f"Provider must be a string, got {type(provider).__name__}.",
            parameter=parameter_name,
            received_value=type(provider).__name__
        )

    provider_clean = provider.strip()

    if not provider_clean:
        raise ProviderValidationError(
            "Provider cannot be an empty string. "
            f"Use one of: {', '.join(sorted(VALID_PROVIDERS))} "
            "or an HTTP URL for custom endpoints.",
            parameter=parameter_name,
            received_value=repr(provider)
        )

    # Check if it's a valid provider name (case-insensitive)
    provider_lower = provider_clean.lower()
    if provider_lower in VALID_PROVIDERS:
        return provider_lower

    # Check if it's an HTTP URL (custom endpoint)
    if provider_clean.startswith("http://") or provider_clean.startswith("https://"):
        return provider_clean  # Return original case for URLs

    raise ProviderValidationError(
        f"Unknown provider '{provider}'. "
        f"Valid providers: {', '.join(sorted(VALID_PROVIDERS))} "
        "or an HTTP/HTTPS URL for custom endpoints.",
        parameter=parameter_name,
        received_value=provider
    )


def validate_model(
    model: Any,
    parameter_name: str = "model"
) -> str:
    """
    Validate model parameter format.

    Note: This only validates the format. The ModelSettings class
    handles alias resolution and provider-specific validation.

    Args:
        model: Model name or alias
        parameter_name: Name of parameter (for error messages)

    Returns:
        str: Validated model string

    Raises:
        ModelValidationError: If validation fails
    """
    if model is None:
        raise ModelValidationError(
            "Model cannot be None. "
            f"Use a model name, alias (e.g., 'gpt', 'claude'), "
            f"or tier (e.g., {', '.join(sorted(VALID_MODEL_TIERS))}).",
            parameter=parameter_name,
            received_value=model
        )

    if not isinstance(model, str):
        raise ModelValidationError(
            f"Model must be a string, got {type(model).__name__}.",
            parameter=parameter_name,
            received_value=type(model).__name__
        )

    model_clean = model.strip()

    if not model_clean:
        raise ModelValidationError(
            "Model cannot be an empty string. "
            f"Use a model name, alias (e.g., 'gpt', 'claude'), "
            f"or tier (e.g., {', '.join(sorted(VALID_MODEL_TIERS))}).",
            parameter=parameter_name,
            received_value=repr(model)
        )

    return model_clean


# =============================================================================
# Batch-Specific Parameter Validation
# =============================================================================

def validate_marker_dataframe(
    marker: Any,
    parameter_name: str = "marker"
) -> Union[pd.DataFrame, str]:
    """
    Validate marker input for batch functions.

    Args:
        marker: DataFrame or path to CSV file
        parameter_name: Name of parameter (for error messages)

    Returns:
        Union[pd.DataFrame, str]: Validated DataFrame or file path

    Raises:
        BatchParameterValidationError: If validation fails
    """
    if marker is None:
        raise BatchParameterValidationError(
            "Marker data cannot be None. "
            "Provide a pandas DataFrame or path to a CSV file.",
            parameter=parameter_name,
            received_value=marker
        )

    # Case 1: DataFrame
    if isinstance(marker, pd.DataFrame):
        if marker.empty:
            raise BatchParameterValidationError(
                "Marker DataFrame is empty. "
                "Please provide a DataFrame with marker data.",
                parameter=parameter_name,
                received_value="empty DataFrame"
            )
        return marker

    # Case 2: File path
    if isinstance(marker, str):
        if not marker.strip():
            raise BatchParameterValidationError(
                "Marker file path cannot be an empty string.",
                parameter=parameter_name,
                received_value=repr(marker)
            )

        if not os.path.exists(marker):
            raise BatchParameterValidationError(
                f"Marker file not found: '{marker}'. "
                "Please check the file path.",
                parameter=parameter_name,
                received_value=marker
            )

        return marker

    raise BatchParameterValidationError(
        f"Invalid marker type. Expected DataFrame or file path string, "
        f"got {type(marker).__name__}.",
        parameter=parameter_name,
        received_value=type(marker).__name__
    )


def validate_positive_int(
    value: Any,
    parameter_name: str,
    allow_zero: bool = False
) -> int:
    """
    Validate a positive integer parameter.

    Args:
        value: Value to validate
        parameter_name: Name of parameter (for error messages)
        allow_zero: Whether to allow zero (default: False)

    Returns:
        int: Validated integer

    Raises:
        BatchParameterValidationError: If validation fails
    """
    if value is None:
        raise BatchParameterValidationError(
            f"{parameter_name} cannot be None.",
            parameter=parameter_name,
            received_value=value
        )

    try:
        int_value = int(value)
    except (TypeError, ValueError):
        raise BatchParameterValidationError(
            f"{parameter_name} must be an integer, got {type(value).__name__}.",
            parameter=parameter_name,
            received_value=value
        )

    min_val = 0 if allow_zero else 1
    if int_value < min_val:
        raise BatchParameterValidationError(
            f"{parameter_name} must be {'non-negative' if allow_zero else 'positive'}, "
            f"got {int_value}.",
            parameter=parameter_name,
            received_value=int_value
        )

    return int_value


def validate_ranking_method(
    ranking_method: Any,
    parameter_name: str = "ranking_method"
) -> str:
    """
    Validate ranking method parameter.

    Args:
        ranking_method: Ranking method string
        parameter_name: Name of parameter (for error messages)

    Returns:
        str: Validated ranking method

    Raises:
        BatchParameterValidationError: If validation fails
    """
    if ranking_method is None:
        raise BatchParameterValidationError(
            f"Ranking method cannot be None. "
            f"Valid methods: {', '.join(sorted(VALID_RANKING_METHODS))}",
            parameter=parameter_name,
            received_value=ranking_method
        )

    if not isinstance(ranking_method, str):
        raise BatchParameterValidationError(
            f"Ranking method must be a string, got {type(ranking_method).__name__}.",
            parameter=parameter_name,
            received_value=type(ranking_method).__name__
        )

    if ranking_method not in VALID_RANKING_METHODS:
        raise BatchParameterValidationError(
            f"Invalid ranking method '{ranking_method}'. "
            f"Valid methods: {', '.join(sorted(VALID_RANKING_METHODS))}",
            parameter=parameter_name,
            received_value=ranking_method
        )

    return ranking_method


def validate_validator_involvement(
    validator_involvement: Any,
    parameter_name: str = "validator_involvement"
) -> str:
    """
    Validate validator_involvement parameter.

    Args:
        validator_involvement: Validator involvement level ('v0' or 'v1')
        parameter_name: Name of parameter (for error messages)

    Returns:
        str: Validated validator involvement level

    Raises:
        CASSIAValidationError: If validation fails
    """
    if validator_involvement is None:
        return "v1"  # Default value

    if not isinstance(validator_involvement, str):
        raise CASSIAValidationError(
            f"validator_involvement must be a string, got {type(validator_involvement).__name__}.",
            parameter=parameter_name,
            received_value=type(validator_involvement).__name__
        )

    if validator_involvement not in VALID_VALIDATOR_INVOLVEMENTS:
        raise CASSIAValidationError(
            f"Invalid validator_involvement '{validator_involvement}'. "
            f"Valid values: {', '.join(sorted(VALID_VALIDATOR_INVOLVEMENTS))}",
            parameter=parameter_name,
            received_value=validator_involvement
        )

    return validator_involvement


# =============================================================================
# High-Level Validation Functions
# =============================================================================

def validate_runCASSIA_inputs(
    model: Any,
    temperature: Any,
    marker_list: Any,
    tissue: Any,
    species: Any,
    provider: Any,
    additional_info: Any = None,
    validator_involvement: str = "v1"
) -> dict:
    """
    Validate all inputs for runCASSIA function.

    Args:
        All parameters from runCASSIA

    Returns:
        dict: Dictionary of validated parameters

    Raises:
        CASSIAValidationError: If any validation fails
    """
    validated = {}

    # Validate in order of likelihood to fail
    validated['marker_list'] = validate_marker_list(marker_list)
    validated['temperature'] = validate_temperature(temperature)
    validated['tissue'] = validate_tissue(tissue)
    validated['species'] = validate_species(species)
    validated['provider'] = validate_provider(provider)
    validated['model'] = validate_model(model)

    # Optional parameters
    if additional_info is not None and not isinstance(additional_info, str):
        raise CASSIAValidationError(
            f"additional_info must be a string or None, "
            f"got {type(additional_info).__name__}.",
            parameter="additional_info",
            received_value=type(additional_info).__name__
        )
    validated['additional_info'] = additional_info

    # Validator involvement
    validated['validator_involvement'] = validate_validator_involvement(validator_involvement)

    return validated


def validate_runCASSIA_batch_inputs(
    marker: Any,
    output_name: Any,
    n_genes: Any,
    model: Any,
    temperature: Any,
    tissue: Any,
    species: Any,
    max_workers: Any,
    provider: Any,
    max_retries: Any,
    ranking_method: Any,
    validator_involvement: str = "v1",
    **kwargs
) -> dict:
    """
    Validate all inputs for runCASSIA_batch function.

    Args:
        All parameters from runCASSIA_batch

    Returns:
        dict: Dictionary of validated parameters

    Raises:
        CASSIAValidationError: If any validation fails
    """
    validated = {}

    # Validate marker data first (most likely to fail)
    validated['marker'] = validate_marker_dataframe(marker)

    # Common parameters
    validated['temperature'] = validate_temperature(temperature)
    validated['tissue'] = validate_tissue(tissue)
    validated['species'] = validate_species(species)
    validated['provider'] = validate_provider(provider)
    validated['model'] = validate_model(model)

    # Batch-specific parameters
    validated['n_genes'] = validate_positive_int(n_genes, "n_genes")
    validated['max_workers'] = validate_positive_int(max_workers, "max_workers")
    validated['max_retries'] = validate_positive_int(max_retries, "max_retries", allow_zero=True)
    validated['ranking_method'] = validate_ranking_method(ranking_method)

    # Validator involvement
    validated['validator_involvement'] = validate_validator_involvement(validator_involvement)

    # Output name
    if output_name is not None and not isinstance(output_name, str):
        raise BatchParameterValidationError(
            f"output_name must be a string, got {type(output_name).__name__}.",
            parameter="output_name",
            received_value=type(output_name).__name__
        )
    validated['output_name'] = output_name

    return validated


def validate_runCASSIA_with_reference_inputs(
    model: Any,
    temperature: Any,
    marker_list: Any,
    tissue: Any,
    species: Any,
    provider: Any,
    additional_info: Any = None,
    validator_involvement: str = "v1",
    reference_threshold: Any = 40,
    **kwargs
) -> dict:
    """
    Validate all inputs for runCASSIA_with_reference function.

    Args:
        All parameters from runCASSIA_with_reference

    Returns:
        dict: Dictionary of validated parameters

    Raises:
        CASSIAValidationError: If any validation fails
    """
    # Start with base runCASSIA validation
    validated = validate_runCASSIA_inputs(
        model=model,
        temperature=temperature,
        marker_list=marker_list,
        tissue=tissue,
        species=species,
        provider=provider,
        additional_info=additional_info,
        validator_involvement=validator_involvement
    )

    # Additional reference-specific validation
    if reference_threshold is not None:
        try:
            threshold = int(reference_threshold)
            if threshold < 0 or threshold > 100:
                raise CASSIAValidationError(
                    f"reference_threshold must be between 0 and 100, got {threshold}.",
                    parameter="reference_threshold",
                    received_value=threshold
                )
            validated['reference_threshold'] = threshold
        except (TypeError, ValueError):
            raise CASSIAValidationError(
                f"reference_threshold must be an integer, got {type(reference_threshold).__name__}.",
                parameter="reference_threshold",
                received_value=reference_threshold
            )

    return validated
