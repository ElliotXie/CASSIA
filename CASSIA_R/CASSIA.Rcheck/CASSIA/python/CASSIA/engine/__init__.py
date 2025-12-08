# CASSIA Engine Module
# Core annotation engine

from .tools_function import (
    runCASSIA,
    runCASSIA_batch,
    runCASSIA_with_reference,
    runCASSIA_batch_with_reference,
    set_api_key,
    set_openai_api_key,
    set_anthropic_api_key,
    set_openrouter_api_key,
)

# Re-export main_function_code components
from .main_function_code import (
    run_cell_type_analysis,
    run_cell_type_analysis_claude,
    run_cell_type_analysis_openrouter,
    run_cell_type_analysis_custom,
    extract_json_from_reply,
)

__all__ = [
    # Main functions
    'runCASSIA',
    'runCASSIA_batch',
    'runCASSIA_with_reference',
    'runCASSIA_batch_with_reference',
    # API key setters
    'set_api_key',
    'set_openai_api_key',
    'set_anthropic_api_key',
    'set_openrouter_api_key',
    # Analysis functions
    'run_cell_type_analysis',
    'run_cell_type_analysis_claude',
    'run_cell_type_analysis_openrouter',
    'run_cell_type_analysis_custom',
    'extract_json_from_reply',
]
