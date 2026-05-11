"""
Reference Agent Module for CASSIA

Provides intelligent reference document retrieval and context injection
for cell type annotation tasks. The agent analyzes marker genes to:
1. Determine preliminary cell type classification
2. Select relevant expert reference files from the router
3. Extract relevant reference content for prompt injection
"""

from .reference_agent import (
    ReferenceAgent,
    get_reference_content,
    get_subcluster_reference_brief,
    format_reference_for_prompt,
)
from .complexity_scorer import select_references_llm
from .reference_selector import select_references
from .section_extractor import extract_sections, parse_markdown

__all__ = [
    'ReferenceAgent',
    'get_reference_content',
    'get_subcluster_reference_brief',
    'format_reference_for_prompt',
    'select_references_llm',
    'select_references',
    'extract_sections',
    'parse_markdown',
]
