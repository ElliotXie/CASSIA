"""
Custom exceptions for CASSIA input validation.

This module provides a hierarchy of exception classes for clear,
actionable error messages when validation fails.
"""


class CASSIAValidationError(ValueError):
    """Base exception for CASSIA validation errors."""

    def __init__(self, message: str, parameter: str = None, received_value=None):
        """
        Initialize validation error.

        Args:
            message: Clear description of the validation failure
            parameter: Name of the parameter that failed validation
            received_value: The invalid value that was received
        """
        self.parameter = parameter
        self.received_value = received_value
        self.message = message
        super().__init__(self._format_message())

    def _format_message(self) -> str:
        """Format the error message with parameter and value info."""
        parts = [f"CASSIA Validation Error: {self.message}"]
        if self.parameter:
            parts.append(f"  Parameter: {self.parameter}")
        if self.received_value is not None:
            # Truncate long values
            val_str = repr(self.received_value)
            if len(val_str) > 100:
                val_str = val_str[:100] + "..."
            parts.append(f"  Received: {val_str}")
        return "\n".join(parts)


class MarkerValidationError(CASSIAValidationError):
    """Raised when marker_list validation fails.

    This includes:
    - None or empty marker list
    - Non-string elements in marker list
    - Too few markers (< 10, warning only)
    - Too many markers (> 500)
    - Invalid gene ID formats (Ensembl, Entrez)
    """
    pass


class TemperatureValidationError(CASSIAValidationError):
    """Raised when temperature validation fails.

    Temperature must be >= 0.
    """
    pass


class ProviderValidationError(CASSIAValidationError):
    """Raised when provider validation fails.

    Valid providers: 'openai', 'anthropic', 'openrouter', or HTTP URL.
    """
    pass


class ModelValidationError(CASSIAValidationError):
    """Raised when model validation fails.

    Model must be a non-empty string.
    """
    pass


class TissueSpeciesValidationError(CASSIAValidationError):
    """Raised when tissue or species validation fails.

    Note: Empty tissue/species generates a warning, not an error.
    This exception is for type errors (e.g., non-string input).
    """
    pass


class BatchParameterValidationError(CASSIAValidationError):
    """Raised when batch-specific parameters are invalid.

    This includes:
    - Invalid marker DataFrame or file path
    - Invalid n_genes, max_workers, max_retries values
    - Invalid ranking_method
    """
    pass


class APIValidationError(CASSIAValidationError):
    """Raised when API key validation fails.

    This includes:
    - Missing API keys
    - Invalid API keys (authentication failures)
    - API quota or rate limit issues during validation
    - Network connectivity issues
    """
    pass
