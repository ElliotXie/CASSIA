# CASSIA Core Module
# Foundation modules for CASSIA

from .logging_config import get_logger
from .exceptions import CASSIAValidationError

# Alias for backward compatibility
CASSIAError = CASSIAValidationError
from .validation import (
    validate_runCASSIA_inputs,
    validate_runCASSIA_batch_inputs,
    validate_runCASSIA_with_reference_inputs
)
from .llm_utils import call_llm
from .model_settings import ModelSettings, resolve_model_name, get_recommended_model
from .progress_tracker import BatchProgressTracker
from .utils import safe_get, natural_sort_key, clean_conversation_history, write_csv
from .marker_utils import (
    split_markers,
    get_top_markers,
    loadmarker,
    list_available_markers,
    _validate_ranking_parameters,
    _prepare_ranking_column,
    _get_sort_direction
)

__all__ = [
    # Logging
    'get_logger',
    # Exceptions
    'CASSIAError',
    'CASSIAValidationError',
    # Validation
    'validate_runCASSIA_inputs',
    'validate_runCASSIA_batch_inputs',
    'validate_runCASSIA_with_reference_inputs',
    # LLM utilities
    'call_llm',
    # Model settings
    'ModelSettings',
    'resolve_model_name',
    'get_recommended_model',
    # Progress tracking
    'BatchProgressTracker',
    # Utilities
    'safe_get',
    'natural_sort_key',
    'clean_conversation_history',
    'write_csv',
    # Marker utilities
    'split_markers',
    'get_top_markers',
    'loadmarker',
    'list_available_markers',
]
