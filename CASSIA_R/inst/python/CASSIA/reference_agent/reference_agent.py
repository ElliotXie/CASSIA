"""
Reference Agent - Main Orchestrator for CASSIA.

Provides intelligent reference document retrieval and context injection
for cell type annotation tasks.

Uses a TWO-STEP ReAct workflow:
    Step 1: LLM assesses complexity and decides if reference is needed
    Step 2: If needed, LLM selects specific references from router
"""

from typing import Dict, List, Optional
from pathlib import Path

try:
    from .complexity_scorer import assess_complexity, quick_complexity_check, assess_complexity_step1, select_references_step2
    from .reference_selector import select_references, find_references_by_category
    from .section_extractor import extract_sections
    from .utils import (
        get_references_dir,
        load_reference_index,
        format_reference_content
    )
except ImportError:
    from complexity_scorer import assess_complexity, quick_complexity_check, assess_complexity_step1, select_references_step2
    from reference_selector import select_references, find_references_by_category
    from section_extractor import extract_sections
    from utils import (
        get_references_dir,
        load_reference_index,
        format_reference_content
    )


class ReferenceAgent:
    """
    Agent for intelligent reference document retrieval and context injection.

    Uses a TWO-STEP ReAct workflow:
    1. Step 1: LLM assesses marker complexity and decides if reference is needed
    2. Step 2: If needed, LLM sees router and selects specific references
    3. Load and extract content from selected references
    4. Return formatted content for prompt injection
    """

    def __init__(
        self,
        reference_dir: Optional[str] = None,
        provider: str = "openrouter",
        model: Optional[str] = None,
        api_key: Optional[str] = None
    ):
        """
        Initialize ReferenceAgent.

        Args:
            reference_dir: Optional custom path to references directory
            provider: LLM provider for complexity assessment
            model: Optional specific model for complexity assessment
            api_key: Optional API key for LLM provider
        """
        if reference_dir:
            self.reference_dir = Path(reference_dir)
        else:
            self.reference_dir = get_references_dir()

        self.provider = provider
        self.model = model
        self.api_key = api_key
        self._index = None

    @property
    def index(self) -> Dict:
        """Lazy load reference index."""
        if self._index is None:
            self._index = load_reference_index()
        return self._index

    def get_reference_for_markers(
        self,
        markers: List[str],
        tissue: Optional[str] = None,
        species: Optional[str] = None,
        depth: str = 'auto',
        max_content_length: Optional[int] = 8000,
        skip_llm: bool = False
    ) -> Dict:
        """
        Main entry point: Get reference content for a marker set.

        Uses TWO-STEP ReAct workflow:
        - Step 1: LLM assesses complexity → returns requires_reference
        - Step 2: If needed, LLM sees router → selects specific references

        Args:
            markers: List of marker genes (top 20 recommended)
            tissue: Tissue type for context
            species: Species for context
            depth: 'summary', 'detailed', or 'auto' (based on complexity)
            max_content_length: Maximum length of extracted content
            skip_llm: Skip LLM complexity assessment (use quick rules only)

        Returns:
            Dict with:
                - should_use_reference: bool
                - content: str (the reference content to inject)
                - references_used: List[str]
                - complexity_score: float
                - preliminary_cell_type: str
                - cell_type_range: List[str]
                - reasoning: str
        """
        # Ensure we have markers
        if not markers:
            return self._empty_result("No markers provided")

        markers = markers[:20]  # Limit to top 20

        # TWO-STEP ASSESSMENT
        if skip_llm:
            # Quick rule-based check (no LLM calls)
            quick_result = quick_complexity_check(markers)
            complexity_result = {
                'complexity_score': 70 if quick_result['likely_complex'] else 30,
                'preliminary_cell_type': quick_result.get('top_category', 'Unknown'),
                'cell_type_range': [],
                'requires_reference': quick_result['likely_complex'],
                'selected_references': [],  # No LLM selection in skip mode
                'reasoning': quick_result.get('reason', '')
            }
            # For skip_llm mode, use traditional reference selection
            if complexity_result['requires_reference']:
                selected_refs = select_references(
                    markers=markers,
                    preliminary_cell_type=complexity_result['preliminary_cell_type'],
                    reference_categories=quick_result.get('suggested_categories', []),
                    max_references=3
                )
                complexity_result['selected_references'] = [
                    ref.get('path', '').replace(str(self.reference_dir) + '/', '').replace(str(self.reference_dir) + '\\', '')
                    for ref in selected_refs
                ]
        else:
            # Full two-step LLM assessment
            # Step 1: Assess complexity (no router)
            # Step 2: If needed, select references (with router)
            complexity_result = assess_complexity(
                markers=markers,
                tissue=tissue,
                species=species,
                provider=self.provider,
                model=self.model,
                api_key=self.api_key
            )

        # Build result
        result = {
            'should_use_reference': False,
            'content': '',
            'references_used': [],
            'complexity_score': complexity_result.get('complexity_score', 50),
            'preliminary_cell_type': complexity_result.get('preliminary_cell_type', 'Unknown'),
            'cell_type_range': complexity_result.get('cell_type_range', []),
            'selected_references': complexity_result.get('selected_references', []),
            'reasoning': complexity_result.get('reasoning', '')
        }

        # Check if LLM decided reference is needed
        if not complexity_result.get('requires_reference', False):
            result['reasoning'] += " LLM determined reference not needed."
            return result

        result['should_use_reference'] = True

        # Get selected reference paths from LLM (Step 2 result)
        selected_paths = complexity_result.get('selected_references', [])

        if not selected_paths:
            result['should_use_reference'] = False
            result['reasoning'] += " No references selected by LLM."
            return result

        # Determine depth
        if depth == 'auto':
            depth = 'detailed' if result['complexity_score'] > 70 else 'summary'

        # Load and extract content from LLM-selected references
        all_content = []
        for ref_path_str in selected_paths:
            # Build full path
            ref_path = self.reference_dir / ref_path_str

            if not ref_path.exists():
                # Try without leading slash
                ref_path_str = ref_path_str.lstrip('/').lstrip('\\')
                ref_path = self.reference_dir / ref_path_str

            if not ref_path.exists():
                continue

            extraction = extract_sections(
                file_path=ref_path,
                cell_type_guess=result['preliminary_cell_type'],
                markers=markers,
                depth=depth
            )

            if extraction['full_content']:
                ref_name = ref_path_str.replace('/', '_').replace('\\', '_').replace('.md', '')
                all_content.append(f"### Reference: {ref_name}\n\n{extraction['full_content']}")
                result['references_used'].append(ref_path_str)

        # Combine and format content
        if all_content:
            combined = '\n\n---\n\n'.join(all_content)
            result['content'] = format_reference_content(combined, max_content_length)
        else:
            result['should_use_reference'] = False
            result['reasoning'] += " No content extracted from selected references."

        return result

    def _empty_result(self, reason: str) -> Dict:
        """Return an empty result with a reason."""
        return {
            'should_use_reference': False,
            'content': '',
            'references_used': [],
            'complexity_score': 0,
            'preliminary_cell_type': 'Unknown',
            'cell_type_range': [],
            'reference_categories': [],
            'reasoning': reason
        }

    def list_available_references(self, category: Optional[str] = None) -> List[Dict]:
        """
        List all available reference documents.

        Args:
            category: Optional category to filter by

        Returns:
            List of reference info dicts
        """
        if category:
            return find_references_by_category(category)

        # List all
        results = []
        if self.index and 'references' in self.index:
            for ref_id, ref_info in self.index['references'].items():
                results.append({
                    'id': ref_id,
                    'category': ref_info.get('category', ''),
                    'cell_types': ref_info.get('cell_types', []),
                    'path': ref_info.get('path', '')
                })

        return results

    def get_reference_content_direct(
        self,
        reference_id: str,
        section_path: Optional[str] = None
    ) -> Optional[str]:
        """
        Get content from a specific reference directly.

        Args:
            reference_id: The reference ID
            section_path: Optional path to specific section (e.g., "Subtypes/Tregs")

        Returns:
            Reference content or None if not found
        """
        if not self.index or 'references' not in self.index:
            return None

        ref_info = self.index['references'].get(reference_id)
        if not ref_info:
            return None

        ref_path = self.reference_dir / ref_info.get('path', '')
        if not ref_path.exists():
            return None

        if section_path:
            from .section_extractor import get_section_by_path
            return get_section_by_path(ref_path, section_path)

        # Return full content
        from .utils import load_markdown_file, parse_yaml_frontmatter
        content = load_markdown_file(ref_path)
        _, body = parse_yaml_frontmatter(content)
        return body


def get_reference_content(
    markers: List[str],
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    provider: str = "openrouter",
    model: Optional[str] = None,
    **kwargs
) -> Dict:
    """
    Convenience function to get reference content for markers.

    Uses TWO-STEP ReAct workflow:
    - Step 1: LLM assesses complexity and decides if reference is needed
    - Step 2: If needed, LLM sees router and selects specific references

    Args:
        markers: List of marker genes
        tissue: Optional tissue type
        species: Optional species
        provider: LLM provider for complexity assessment
        model: Optional specific model
        **kwargs: Additional arguments passed to ReferenceAgent

    Returns:
        Reference result dictionary
    """
    agent = ReferenceAgent(provider=provider, model=model)
    return agent.get_reference_for_markers(
        markers=markers,
        tissue=tissue,
        species=species,
        **kwargs
    )


def format_reference_for_prompt(reference_result: Dict) -> str:
    """
    Format reference result for injection into annotation prompts.

    Args:
        reference_result: Result from get_reference_for_markers

    Returns:
        Formatted string ready for prompt injection
    """
    if not reference_result.get('should_use_reference') or not reference_result.get('content'):
        return ""

    return f"""<expert_reference>
The following expert-curated reference information is provided to assist with cell type annotation.
Use this technical guidance to inform your analysis, particularly for subtype differentiation
and marker interpretation.

Preliminary assessment: {reference_result.get('preliminary_cell_type', 'Unknown')}
Complexity score: {reference_result.get('complexity_score', 'N/A')}/100

{reference_result['content']}
</expert_reference>"""
