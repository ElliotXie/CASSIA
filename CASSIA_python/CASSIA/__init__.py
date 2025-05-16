from .main_function_code import *
from .tools_function import *

# Explicitly import the functions we want to expose
from .tools_function import (
    loadmarker,
    list_available_markers,
    runCASSIA,
    runCASSIA_batch,
    runCASSIA_batch_n_times,
    runCASSIA_n_times,
    runCASSIA_similarity_score_batch,
    runCASSIA_n_times_similarity_score,
    runCASSIA_score_batch,
    runCASSIA_generate_score_report,
    runCASSIA_subclusters,
    runCASSIA_n_subcluster,
    runCASSIA_pipeline,
    compareCelltypes,
    set_openai_api_key,
    set_anthropic_api_key,
    set_openrouter_api_key,
    set_api_key
)

# Import the new LLM utilities
from .llm_utils import call_llm

# Import the annotation boost functionality
from .annotation_boost import iterative_marker_analysis
from .annotation_boost import runCASSIA_annotationboost, runCASSIA_annotationboost_additional_task

__version__ = "0.1.9"