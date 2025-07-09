from .main_function_code import *
from .tools_function import *

# Explicitly import the functions we want to expose
from .tools_function import (
    loadmarker,
    list_available_markers,
    runCASSIA,
    runCASSIA_batch,
    runCASSIA_score_batch,
    runCASSIA_generate_score_report,
    runCASSIA_pipeline,
    set_openai_api_key,
    set_anthropic_api_key,
    set_openrouter_api_key,
    set_api_key
)

# Import cell type comparison from its own module
from .cell_type_comparison import compareCelltypes

# Import the new Symphony Compare function
from .symphony_compare import symphonyCompare

# Import the new LLM utilities
from .llm_utils import call_llm

# Import the annotation boost functionality
from .annotation_boost import iterative_marker_analysis
from .annotation_boost import runCASSIA_annotationboost, runCASSIA_annotationboost_additional_task

# Import main functions from extracted modules
try:
    from .merging_annotation import merge_annotations, merge_annotations_all
    # Internal utilities used by tutorial (not part of public API)
    from .merging_annotation import _create_annotation_prompt, _parse_llm_response
except ImportError:
    pass  # Module may not be available in all installations

try:
    from .Uncertainty_quantification import (
        runCASSIA_batch_n_times,
        runCASSIA_similarity_score_batch,
        runCASSIA_n_times_similarity_score
    )
except ImportError:
    pass  # Module may not be available in all installations

# Import subclustering if available
try:
    from .subclustering import (
        runCASSIA_subclusters,
        runCASSIA_n_subcluster,
        annotate_subclusters
    )
except ImportError:
    pass  # Module may not be available in all installations

__version__ = "0.2.21"
