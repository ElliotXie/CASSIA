# CASSIA Comparison Module
# Multi-model comparison

from .symphony_compare import symphonyCompare

# Alias for backward compatibility (capital S was used in some imports)
SymphonyCompare = symphonyCompare

__all__ = [
    'symphonyCompare',
    'SymphonyCompare',
]
