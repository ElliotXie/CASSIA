"""
Reference Agent Module for CASSIA

Provides intelligent reference document retrieval and context injection
for cell type annotation tasks. The agent analyzes marker genes to:
1. Determine preliminary cell type classification
2. Assess complexity/ambiguity using LLM
3. Decide whether reference retrieval is needed
4. Select and extract relevant reference content
"""

from .reference_agent import ReferenceAgent, get_reference_content
from .complexity_scorer import assess_complexity
from .reference_selector import select_references
from .section_extractor import extract_sections, parse_markdown

__all__ = [
    'ReferenceAgent',
    'get_reference_content',
    'assess_complexity',
    'select_references',
    'extract_sections',
    'parse_markdown',
]
