# Annotation Boost Agent

from .annotation_boost import (
    runCASSIA_annotationboost,
    runCASSIA_annotationboost_additional_task,
    iterative_marker_analysis
)
from .super_annotation_boost import runSuperAnnotationBoost

# Alias for backward compatibility
runCASSIA_super_annotation_boost = runSuperAnnotationBoost

__all__ = [
    'runCASSIA_annotationboost',
    'runCASSIA_annotationboost_additional_task',
    'iterative_marker_analysis',
    'runSuperAnnotationBoost',
    'runCASSIA_super_annotation_boost',
]
