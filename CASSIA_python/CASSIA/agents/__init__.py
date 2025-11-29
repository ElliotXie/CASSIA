# CASSIA Agents Module
# Specialized agent modules

# Import from submodules
from .annotation_boost.annotation_boost import (
    runCASSIA_annotationboost,
    runCASSIA_annotationboost_additional_task,
    iterative_marker_analysis
)
from .annotation_boost import runCASSIA_super_annotation_boost, runSuperAnnotationBoost

from .uncertainty.Uncertainty_quantification import (
    runCASSIA_n_times,
    runCASSIA_batch_n_times,
    runCASSIA_similarity_score_batch,
    runCASSIA_n_times_similarity_score
)

from .merging.merging_annotation import merge_annotations, merge_annotations_all

from .subclustering.subclustering import runCASSIA_subclusters, annotate_subclusters, runCASSIA_n_subcluster

# Alias for backward compatibility
runCASSIA_subclustering = runCASSIA_subclusters

# Reference agent
from .reference_agent import ReferenceAgent, get_reference_content, format_reference_for_prompt

__all__ = [
    # Annotation boost
    'runCASSIA_annotationboost',
    'runCASSIA_annotationboost_additional_task',
    'iterative_marker_analysis',
    'runCASSIA_super_annotation_boost',
    'runSuperAnnotationBoost',
    # Uncertainty
    'runCASSIA_n_times',
    'runCASSIA_batch_n_times',
    'runCASSIA_similarity_score_batch',
    'runCASSIA_n_times_similarity_score',
    # Merging
    'merge_annotations',
    'merge_annotations_all',
    # Subclustering
    'runCASSIA_subclusters',
    'runCASSIA_subclustering',  # Alias
    'annotate_subclusters',
    'runCASSIA_n_subcluster',
    # Reference agent
    'ReferenceAgent',
    'get_reference_content',
    'format_reference_for_prompt',
]
