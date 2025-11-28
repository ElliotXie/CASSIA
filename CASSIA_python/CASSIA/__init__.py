# CASSIA - Cell Annotation with Semantic Similarity for Intelligent Analysis
# Root module with backward-compatible exports from reorganized submodules

__version__ = "0.3.1.dev4"

# =============================================================================
# BACKWARD COMPATIBILITY LAYER
# =============================================================================
# This file re-exports all public functions from the new organized submodules
# to maintain backward compatibility with existing code that uses:
#   from CASSIA import runCASSIA_batch
#   from CASSIA import set_api_key
# etc.

# Import logging configuration (errors visible by default)
try:
    from .core.logging_config import get_logger, set_log_level, warn_user
except ImportError:
    from .logging_config import get_logger, set_log_level, warn_user

# -----------------------------------------------------------------------------
# CORE ANALYSIS FUNCTIONS (most commonly used)
# -----------------------------------------------------------------------------
try:
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
except ImportError:
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

# Main function code exports
try:
    from .engine.main_function_code import *
except ImportError:
    from .main_function_code import *

# -----------------------------------------------------------------------------
# MARKER UTILITIES
# -----------------------------------------------------------------------------
try:
    from .core.marker_utils import loadmarker, list_available_markers, split_markers, get_top_markers
except ImportError:
    from .marker_utils import loadmarker, list_available_markers, split_markers, get_top_markers

# -----------------------------------------------------------------------------
# LLM UTILITIES
# -----------------------------------------------------------------------------
try:
    from .core.llm_utils import call_llm
except ImportError:
    from .llm_utils import call_llm

# -----------------------------------------------------------------------------
# MODEL SETTINGS (with fuzzy alias support)
# -----------------------------------------------------------------------------
try:
    from .core.model_settings import (
        ModelSettings,
        resolve_model_name,
        get_recommended_model,
        get_available_aliases,
        print_available_models,
        print_available_aliases
    )
except ImportError:
    from .model_settings import (
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
try:
    from .core.exceptions import (
        CASSIAValidationError,
        MarkerValidationError,
        TemperatureValidationError,
        ProviderValidationError,
        ModelValidationError,
        TissueSpeciesValidationError,
        BatchParameterValidationError
    )
except ImportError:
    from .exceptions import (
        CASSIAValidationError,
        MarkerValidationError,
        TemperatureValidationError,
        ProviderValidationError,
        ModelValidationError,
        TissueSpeciesValidationError,
        BatchParameterValidationError
    )

# Alias for backward compatibility
CASSIAError = CASSIAValidationError

try:
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
except ImportError:
    from .validation import (
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
# SCORING AND EVALUATION
# -----------------------------------------------------------------------------
try:
    from .evaluation.scoring import (
        runCASSIA_score_batch,
        score_annotation_batch,
        prompt_creator_score,
        extract_score_and_reasoning,
        score_single_analysis,
        process_single_row
    )
except ImportError:
    from .scoring import (
        runCASSIA_score_batch,
        score_annotation_batch,
        prompt_creator_score,
        extract_score_and_reasoning,
        score_single_analysis,
        process_single_row
    )

try:
    from .evaluation.cell_type_comparison import compareCelltypes
except ImportError:
    from .cell_type_comparison import compareCelltypes

# -----------------------------------------------------------------------------
# MULTI-MODEL COMPARISON
# -----------------------------------------------------------------------------
try:
    from .comparison import symphonyCompare, SymphonyCompare
except ImportError:
    try:
        from .symphony_compare import symphonyCompare
        SymphonyCompare = symphonyCompare  # Alias
    except ImportError:
        pass  # Module may not be available

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
except ImportError:
    try:
        from .generate_reports import (
            generate_analysis_html_report,
            generate_html_report,
            process_single_report,
            generate_score_index_page,
            runCASSIA_generate_score_report,
            generate_subclustering_report,
            calculate_evaluation_metrics
        )
    except ImportError:
        pass  # Module may not be available

try:
    from .reports.generate_batch_report import generate_batch_html_report_from_data
except ImportError:
    try:
        from .generate_batch_report import generate_batch_html_report_from_data
    except ImportError:
        pass

# -----------------------------------------------------------------------------
# PIPELINE
# -----------------------------------------------------------------------------
try:
    from .pipeline.pipeline import runCASSIA_pipeline
except ImportError:
    try:
        from .pipeline import runCASSIA_pipeline
    except ImportError:
        pass  # Module may not be available

# -----------------------------------------------------------------------------
# ANNOTATION BOOST AGENT
# -----------------------------------------------------------------------------
try:
    from .agents.annotation_boost.annotation_boost import (
        runCASSIA_annotationboost,
        runCASSIA_annotationboost_additional_task,
        iterative_marker_analysis
    )
except ImportError:
    try:
        from .annotation_boost import (
            runCASSIA_annotationboost,
            runCASSIA_annotationboost_additional_task,
            iterative_marker_analysis
        )
    except ImportError:
        pass  # Module may not be available

# -----------------------------------------------------------------------------
# MERGING AGENT
# -----------------------------------------------------------------------------
try:
    from .agents.merging.merging_annotation import merge_annotations, merge_annotations_all
    from .agents.merging.merging_annotation import _create_annotation_prompt, _parse_llm_response
except ImportError:
    try:
        from .merging_annotation import merge_annotations, merge_annotations_all
        from .merging_annotation import _create_annotation_prompt, _parse_llm_response
    except ImportError:
        pass  # Module may not be available

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
except ImportError:
    try:
        from .Uncertainty_quantification import (
            runCASSIA_n_times,
            runCASSIA_batch_n_times,
            runCASSIA_similarity_score_batch,
            runCASSIA_n_times_similarity_score
        )
    except ImportError:
        pass  # Module may not be available

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
except ImportError:
    try:
        from .subclustering import (
            runCASSIA_subclusters,
            runCASSIA_n_subcluster,
            annotate_subclusters
        )
        runCASSIA_subclustering = runCASSIA_subclusters  # Alias
    except ImportError:
        pass  # Module may not be available

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
except ImportError:
    try:
        from .reference_agent import (
            ReferenceAgent,
            get_reference_content,
            format_reference_for_prompt,
            assess_complexity,
            select_references
        )
        from .reference_agent.complexity_scorer import (
            assess_complexity_step1,
            select_references_step2
        )
    except ImportError:
        pass  # Reference agent module may not be available

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
except ImportError:
    try:
        from .hypothesis_generation import generate_hypothesis, process_marker_file, run_multi_analysis
        from .summarize_hypothesis_runs import summarize_runs
        # Create aliases for backward compatibility
        runCASSIA_hypothesis = generate_hypothesis
        run_hypothesis_analysis = run_multi_analysis
        summarize_hypothesis_runs = summarize_runs
    except ImportError:
        pass  # Module may not be available

# -----------------------------------------------------------------------------
# IMAGE ANALYSIS
# -----------------------------------------------------------------------------
try:
    from .imaging.llm_image import analyze_image, call_llm_with_image
except ImportError:
    try:
        from .llm_image import analyze_image, call_llm_with_image
    except ImportError:
        pass  # Module may not be available

# -----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------
try:
    from .core.utils import safe_get, natural_sort_key, clean_conversation_history, write_csv
except ImportError:
    from .utils import safe_get, natural_sort_key, clean_conversation_history, write_csv

try:
    from .core.progress_tracker import BatchProgressTracker
except ImportError:
    from .progress_tracker import BatchProgressTracker

# End of __init__.py
