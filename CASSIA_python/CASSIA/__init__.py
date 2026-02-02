# CASSIA - Cell Annotation with Semantic Similarity for Intelligent Analysis
# Root module with backward-compatible exports from reorganized submodules

__version__ = "0.3.18"

# =============================================================================
# BACKWARD COMPATIBILITY LAYER
# =============================================================================
# This file re-exports all public functions from the new organized submodules
# to maintain backward compatibility with existing code that uses:
#   from CASSIA import runCASSIA_batch
#   from CASSIA import set_api_key
# etc.

# Import logging configuration (errors visible by default)
from .core.logging_config import get_logger, set_log_level, warn_user

# Track module availability for diagnostics
_module_availability = {}
_import_errors = {}

# -----------------------------------------------------------------------------
# CORE ANALYSIS FUNCTIONS (most commonly used)
# -----------------------------------------------------------------------------
from .engine.tools_function import (
    runCASSIA,
    runCASSIA_batch,
    runCASSIA_with_reference,
    runCASSIA_batch_with_reference,
    set_api_key,
    set_openai_api_key,
    set_anthropic_api_key,
    set_openrouter_api_key,
)

# Main function code exports
from .engine.main_function_code import *

# -----------------------------------------------------------------------------
# MARKER UTILITIES
# -----------------------------------------------------------------------------
from .core.marker_utils import loadmarker, list_available_markers, split_markers, get_top_markers

# -----------------------------------------------------------------------------
# LLM UTILITIES
# -----------------------------------------------------------------------------
from .core.llm_utils import call_llm

# -----------------------------------------------------------------------------
# API KEY VALIDATION
# -----------------------------------------------------------------------------
from .core.api_validation import validate_api_keys, clear_validation_cache

# -----------------------------------------------------------------------------
# FREE API ACCESS (for users without API keys)
# -----------------------------------------------------------------------------
try:
    from .core.free_api import (
        get_free_api_key,
        is_free_api_available,
        get_free_api_usage,
        get_free_api_info,
        clear_free_api_cache,
        is_using_free_api,
        MAX_CLUSTERS_PER_JOB,
        MAX_JOBS_PER_DAY,
        FREE_API_PROVIDERS
    )
except ImportError:
    pass  # requests package may not be installed

# -----------------------------------------------------------------------------
# MODEL SETTINGS (with fuzzy alias support)
# -----------------------------------------------------------------------------
from .core.model_settings import (
    ModelSettings,
    resolve_model_name,
    get_recommended_model,
    get_available_aliases,
    print_available_models,
    print_available_aliases
)

# -----------------------------------------------------------------------------
# VALIDATION AND EXCEPTIONS
# -----------------------------------------------------------------------------
from .core.exceptions import (
    CASSIAValidationError,
    MarkerValidationError,
    TemperatureValidationError,
    ProviderValidationError,
    ModelValidationError,
    TissueSpeciesValidationError,
    BatchParameterValidationError,
    APIValidationError
)

# Alias for backward compatibility
CASSIAError = CASSIAValidationError

from .core.validation import (
    validate_marker_list,
    validate_temperature,
    validate_tissue,
    validate_species,
    validate_provider,
    validate_model,
    validate_marker_dataframe,
    validate_runCASSIA_inputs,
    validate_runCASSIA_batch_inputs,
    validate_runCASSIA_with_reference_inputs
)

# -----------------------------------------------------------------------------
# GENE ID CONVERSION (Ensembl/Entrez to symbols)
# -----------------------------------------------------------------------------
from .core.gene_id_converter import (
    convert_gene_ids,
    convert_dataframe_gene_ids,
    is_mygene_available
)

# -----------------------------------------------------------------------------
# SCORING AND EVALUATION
# -----------------------------------------------------------------------------
from .evaluation.scoring import (
    runCASSIA_score_batch,
    score_annotation_batch,
    prompt_creator_score,
    extract_score_and_reasoning,
    score_single_analysis,
    process_single_row
)

from .evaluation.cell_type_comparison import compareCelltypes

# -----------------------------------------------------------------------------
# MULTI-MODEL COMPARISON
# -----------------------------------------------------------------------------
try:
    from .comparison import symphonyCompare, SymphonyCompare
    _module_availability['comparison'] = True
except ImportError as e:
    _module_availability['comparison'] = False
    _import_errors['comparison'] = str(e)

# -----------------------------------------------------------------------------
# REPORT GENERATION
# -----------------------------------------------------------------------------
try:
    from .reports.generate_reports import (
        generate_analysis_html_report,
        generate_html_report,
        process_single_report,
        generate_score_index_page,
        runCASSIA_generate_score_report,
        generate_subclustering_report,
        calculate_evaluation_metrics
    )
    _module_availability['reports'] = True
except ImportError as e:
    _module_availability['reports'] = False
    _import_errors['reports'] = str(e)

try:
    from .reports.generate_batch_report import generate_batch_html_report_from_data
    _module_availability['batch_report'] = True
except ImportError as e:
    _module_availability['batch_report'] = False
    _import_errors['batch_report'] = str(e)

try:
    from .reports.generate_report_uncertainty import generate_uq_html_report
    _module_availability['uncertainty_report'] = True
except ImportError as e:
    _module_availability['uncertainty_report'] = False
    _import_errors['uncertainty_report'] = str(e)

# -----------------------------------------------------------------------------
# PIPELINE
# -----------------------------------------------------------------------------
try:
    from .pipeline.pipeline import runCASSIA_pipeline
    _module_availability['pipeline'] = True
except ImportError as e:
    _module_availability['pipeline'] = False
    _import_errors['pipeline'] = str(e)

# -----------------------------------------------------------------------------
# ANNOTATION BOOST AGENT
# -----------------------------------------------------------------------------
try:
    from .agents.annotation_boost.annotation_boost import (
        runCASSIA_annotationboost,
        runCASSIA_annotationboost_additional_task,
        iterative_marker_analysis
    )
    _module_availability['annotation_boost'] = True
except ImportError as e:
    _module_availability['annotation_boost'] = False
    _import_errors['annotation_boost'] = str(e)

# -----------------------------------------------------------------------------
# MERGING AGENT
# -----------------------------------------------------------------------------
try:
    from .agents.merging.merging_annotation import merge_annotations, merge_annotations_all
    from .agents.merging.merging_annotation import _create_annotation_prompt, _parse_llm_response
    _module_availability['merging'] = True
except ImportError as e:
    _module_availability['merging'] = False
    _import_errors['merging'] = str(e)

# -----------------------------------------------------------------------------
# UNCERTAINTY QUANTIFICATION AGENT
# -----------------------------------------------------------------------------
try:
    from .agents.uncertainty.Uncertainty_quantification import (
        runCASSIA_n_times,
        runCASSIA_batch_n_times,
        runCASSIA_similarity_score_batch,
        runCASSIA_n_times_similarity_score
    )
    _module_availability['uncertainty'] = True
except ImportError as e:
    _module_availability['uncertainty'] = False
    _import_errors['uncertainty'] = str(e)

# -----------------------------------------------------------------------------
# SUBCLUSTERING AGENT
# -----------------------------------------------------------------------------
try:
    from .agents.subclustering import (
        runCASSIA_subclusters,
        runCASSIA_subclustering,  # Alias from __init__.py
        runCASSIA_n_subcluster,
        annotate_subclusters
    )
    _module_availability['subclustering'] = True
except ImportError as e:
    _module_availability['subclustering'] = False
    _import_errors['subclustering'] = str(e)

# -----------------------------------------------------------------------------
# REFERENCE AGENT (intelligent reference retrieval)
# -----------------------------------------------------------------------------
try:
    from .agents.reference_agent import (
        ReferenceAgent,
        get_reference_content,
        format_reference_for_prompt,
        assess_complexity,
        select_references
    )
    from .agents.reference_agent.complexity_scorer import (
        assess_complexity_step1,
        select_references_step2
    )
    _module_availability['reference_agent'] = True
except ImportError as e:
    _module_availability['reference_agent'] = False
    _import_errors['reference_agent'] = str(e)

# -----------------------------------------------------------------------------
# HYPOTHESIS GENERATION
# -----------------------------------------------------------------------------
try:
    from .hypothesis import (
        generate_hypothesis,
        process_marker_file,
        run_multi_analysis,
        summarize_runs,
        # Aliases
        runCASSIA_hypothesis,
        run_hypothesis_analysis,
        summarize_hypothesis_runs
    )
    _module_availability['hypothesis'] = True
except ImportError as e:
    _module_availability['hypothesis'] = False
    _import_errors['hypothesis'] = str(e)

# -----------------------------------------------------------------------------
# IMAGE ANALYSIS
# -----------------------------------------------------------------------------
try:
    from .imaging.llm_image import analyze_image, call_llm_with_image
    _module_availability['imaging'] = True
except ImportError as e:
    _module_availability['imaging'] = False
    _import_errors['imaging'] = str(e)

# -----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------
from .core.utils import safe_get, natural_sort_key, clean_conversation_history, write_csv
from .core.progress_tracker import BatchProgressTracker

# -----------------------------------------------------------------------------
# ANNDATA INTEGRATION (optional Scanpy dependency)
# -----------------------------------------------------------------------------
try:
    from .core.anndata_utils import add_cassia_to_anndata, enhance_scanpy_markers
    _module_availability['anndata'] = True
except ImportError as e:
    _module_availability['anndata'] = False
    _import_errors['anndata'] = str(e)


# =============================================================================
# DIAGNOSTIC FUNCTIONS
# =============================================================================

def get_available_modules():
    """
    Get a dictionary of available optional modules.

    Returns:
        dict: Module names mapped to boolean availability status.

    Example:
        >>> from CASSIA import get_available_modules
        >>> modules = get_available_modules()
        >>> if modules.get('pipeline'):
        ...     from CASSIA import runCASSIA_pipeline
    """
    return _module_availability.copy()


def get_import_errors():
    """
    Get detailed import error messages for modules that failed to load.

    Returns:
        dict: Module names mapped to error messages for failed imports.

    Example:
        >>> from CASSIA import get_import_errors
        >>> errors = get_import_errors()
        >>> for module, error in errors.items():
        ...     print(f"{module}: {error}")
    """
    return _import_errors.copy()


def check_module_available(module_name):
    """
    Check if a specific optional module is available.

    Args:
        module_name: Name of the module to check (e.g., 'pipeline', 'comparison')

    Returns:
        bool: True if module is available, False otherwise.

    Example:
        >>> from CASSIA import check_module_available
        >>> if check_module_available('pipeline'):
        ...     print("Pipeline is available!")
    """
    return _module_availability.get(module_name, False)


# End of __init__.py
