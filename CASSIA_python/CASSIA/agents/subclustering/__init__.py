# Subclustering Agent

from .subclustering import (
    runCASSIA_subclusters,
    annotate_subclusters,
    runCASSIA_n_subcluster,
    build_subcluster_reference_context,
)
from .auto_split import runCASSIA_subclusters_auto_split

# Alias for backward compatibility
runCASSIA_subclustering = runCASSIA_subclusters
runCASSIA_subclustering_auto_split = runCASSIA_subclusters_auto_split

__all__ = [
    'runCASSIA_subclusters',
    'runCASSIA_subclustering',  # Alias
    'runCASSIA_subclusters_auto_split',
    'runCASSIA_subclustering_auto_split',
    'annotate_subclusters',
    'runCASSIA_n_subcluster',
    'build_subcluster_reference_context',
]
