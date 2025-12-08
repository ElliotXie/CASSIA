# CASSIA Imaging Module
# Image processing

from .llm_image import call_llm_with_image, call_llm_analyze_image

# Alias for backward compatibility
analyze_image = call_llm_analyze_image

__all__ = [
    'call_llm_with_image',
    'call_llm_analyze_image',
    'analyze_image',
]
