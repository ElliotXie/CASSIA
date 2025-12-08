"""
Reference Selector for Reference Agent.

Selects and ranks relevant reference documents based on markers,
cell type classification, and reference metadata.
"""

import os
from typing import Dict, List, Optional, Tuple
from pathlib import Path

try:
    from .utils import (
        get_references_dir,
        load_reference_index,
        load_markdown_file,
        parse_yaml_frontmatter,
        normalize_marker_list,
        calculate_marker_overlap,
        list_reference_files
    )
except ImportError:
    from utils import (
        get_references_dir,
        load_reference_index,
        load_markdown_file,
        parse_yaml_frontmatter,
        normalize_marker_list,
        calculate_marker_overlap,
        list_reference_files
    )


class ReferenceCandidate:
    """Represents a candidate reference document with relevance scoring."""

    def __init__(self, ref_id: str, path: Path, metadata: Dict):
        self.id = ref_id
        self.path = path
        self.metadata = metadata
        self.score = 0.0
        self.match_reasons: List[str] = []

    def __repr__(self):
        return f"ReferenceCandidate(id='{self.id}', score={self.score:.2f})"


def select_references(
    markers: List[str],
    preliminary_cell_type: Optional[str] = None,
    reference_categories: Optional[List[str]] = None,
    specific_reference_paths: Optional[List[str]] = None,
    max_references: int = 3,
    min_score: float = 10.0
) -> List[Dict]:
    """
    Select relevant references based on markers and classification.

    Args:
        markers: List of marker genes
        preliminary_cell_type: Initial cell type guess from complexity assessment
        reference_categories: Suggested categories to search (e.g., ['t_cell', 'myeloid'])
        specific_reference_paths: Direct paths to references suggested by LLM (e.g., ['t_cell/cd4/_overview.md'])
        max_references: Maximum number of references to return
        min_score: Minimum relevance score to include

    Returns:
        List of dicts with:
            - id: Reference ID
            - path: Path to reference file
            - score: Relevance score
            - match_reasons: Why this reference was selected
    """
    candidates = []
    refs_dir = get_references_dir()

    # Priority 1: Use specific paths if provided by LLM
    if specific_reference_paths:
        for ref_path in specific_reference_paths:
            full_path = refs_dir / ref_path
            if full_path.exists():
                ref_id = ref_path.replace('/', '_').replace('.md', '')
                candidate = ReferenceCandidate(ref_id, full_path, {})
                candidate.score = 100  # High score for LLM-specified paths
                candidate.match_reasons.append("llm_specified_path")
                candidates.append(candidate)

    # Load index if available
    index = load_reference_index()

    if index and 'references' in index:
        # Score references from index
        candidates.extend(_score_indexed_references(
            index['references'],
            markers,
            preliminary_cell_type,
            reference_categories
        ))

    # Also scan filesystem for references not in index
    filesystem_candidates = _scan_filesystem_references(
        markers,
        preliminary_cell_type,
        reference_categories
    )

    # Merge candidates (prefer indexed ones)
    existing_ids = {c.id for c in candidates}
    for fc in filesystem_candidates:
        if fc.id not in existing_ids:
            candidates.append(fc)

    # Sort by score (descending) and filter
    candidates.sort(key=lambda x: -x.score)
    filtered = [c for c in candidates if c.score >= min_score]

    # Return top N
    return [
        {
            'id': c.id,
            'path': str(c.path),
            'score': c.score,
            'match_reasons': c.match_reasons,
            'metadata': c.metadata
        }
        for c in filtered[:max_references]
    ]


def _score_indexed_references(
    references: Dict,
    markers: List[str],
    preliminary_cell_type: Optional[str],
    reference_categories: Optional[List[str]]
) -> List[ReferenceCandidate]:
    """Score references from the index."""
    candidates = []
    refs_dir = get_references_dir()
    markers_normalized = set(normalize_marker_list(markers[:20]))

    for ref_id, ref_info in references.items():
        path = refs_dir / ref_info.get('path', '')
        if not path.exists():
            continue

        candidate = ReferenceCandidate(ref_id, path, ref_info)

        # Score based on trigger markers
        trigger_markers = set(normalize_marker_list(ref_info.get('trigger_markers', [])))
        trigger_overlap = len(markers_normalized & trigger_markers)
        if trigger_overlap > 0:
            candidate.score += trigger_overlap * 15
            candidate.match_reasons.append(f"trigger_markers_match:{trigger_overlap}")

        # Penalize exclusion marker matches
        exclusion_markers = set(normalize_marker_list(ref_info.get('exclusion_markers', [])))
        exclusion_overlap = len(markers_normalized & exclusion_markers)
        if exclusion_overlap > 0:
            candidate.score -= exclusion_overlap * 10
            candidate.match_reasons.append(f"exclusion_markers_penalty:{exclusion_overlap}")

        # Score based on cell type match
        ref_cell_types = [ct.lower() for ct in ref_info.get('cell_types', [])]
        if preliminary_cell_type:
            cell_type_lower = preliminary_cell_type.lower()
            for ref_ct in ref_cell_types:
                if cell_type_lower in ref_ct or ref_ct in cell_type_lower:
                    candidate.score += 25
                    candidate.match_reasons.append(f"cell_type_match:{ref_ct}")
                    break

        # Score based on category match
        ref_category = ref_info.get('category', '').lower()
        if reference_categories:
            for cat in reference_categories:
                if cat.lower() == ref_category or cat.lower() in ref_category:
                    candidate.score += 20
                    candidate.match_reasons.append(f"category_match:{cat}")
                    break

        # Bonus for priority
        priority = ref_info.get('priority', 5)
        candidate.score += max(0, 10 - priority * 2)

        if candidate.score > 0:
            candidates.append(candidate)

    return candidates


def _scan_filesystem_references(
    markers: List[str],
    preliminary_cell_type: Optional[str],
    reference_categories: Optional[List[str]]
) -> List[ReferenceCandidate]:
    """Scan filesystem for references and score them based on frontmatter."""
    candidates = []
    refs_dir = get_references_dir()

    if not refs_dir.exists():
        return candidates

    # If categories specified, only scan those directories
    if reference_categories:
        search_dirs = []
        for cat in reference_categories:
            cat_dir = refs_dir / cat
            if cat_dir.exists():
                search_dirs.append(cat_dir)
        if not search_dirs:
            search_dirs = [refs_dir]
    else:
        search_dirs = [refs_dir]

    markers_normalized = set(normalize_marker_list(markers[:20]))

    for search_dir in search_dirs:
        for md_file in search_dir.rglob("*.md"):
            # Generate ID from path
            rel_path = md_file.relative_to(refs_dir)
            ref_id = str(rel_path.with_suffix('')).replace(os.sep, '_')

            # Load and parse frontmatter
            content = load_markdown_file(md_file)
            if not content:
                continue

            frontmatter, _ = parse_yaml_frontmatter(content)

            candidate = ReferenceCandidate(ref_id, md_file, frontmatter)

            # Score based on frontmatter trigger markers
            trigger_markers = set(normalize_marker_list(frontmatter.get('trigger_markers', [])))
            trigger_overlap = len(markers_normalized & trigger_markers)
            if trigger_overlap > 0:
                candidate.score += trigger_overlap * 15
                candidate.match_reasons.append(f"trigger_markers_match:{trigger_overlap}")

            # Score based on cell types in frontmatter
            ref_cell_types = [ct.lower() for ct in frontmatter.get('cell_types', [])]
            if preliminary_cell_type:
                cell_type_lower = preliminary_cell_type.lower()
                for ref_ct in ref_cell_types:
                    if cell_type_lower in ref_ct or ref_ct in cell_type_lower:
                        candidate.score += 25
                        candidate.match_reasons.append(f"cell_type_match:{ref_ct}")
                        break

            # Score based on category in frontmatter
            ref_category = frontmatter.get('category', '').lower()
            if reference_categories:
                for cat in reference_categories:
                    if cat.lower() == ref_category:
                        candidate.score += 20
                        candidate.match_reasons.append(f"category_match:{cat}")
                        break

            # Score based on filename/path matching cell type
            path_str = str(md_file).lower()
            if preliminary_cell_type:
                # Check if cell type keywords appear in path
                cell_type_words = preliminary_cell_type.lower().split()
                for word in cell_type_words:
                    if len(word) > 2 and word in path_str:
                        candidate.score += 10
                        candidate.match_reasons.append(f"path_match:{word}")

            if candidate.score > 0:
                candidates.append(candidate)

    return candidates


def get_reference_hierarchy(reference_id: str) -> List[str]:
    """
    Get the hierarchy of references (parent chain) for a given reference.

    Args:
        reference_id: The reference ID to get hierarchy for

    Returns:
        List of reference IDs from root to the given reference
    """
    index = load_reference_index()

    if not index or 'references' not in index:
        return [reference_id]

    hierarchy = []
    current_id = reference_id

    while current_id:
        hierarchy.insert(0, current_id)
        ref_info = index['references'].get(current_id, {})
        current_id = ref_info.get('parent')

    return hierarchy


def get_child_references(reference_id: str) -> List[str]:
    """
    Get all child references of a given reference.

    Args:
        reference_id: The parent reference ID

    Returns:
        List of child reference IDs
    """
    index = load_reference_index()

    if not index or 'references' not in index:
        return []

    ref_info = index['references'].get(reference_id, {})
    return ref_info.get('children', [])


def find_references_by_category(category: str) -> List[Dict]:
    """
    Find all references in a given category.

    Args:
        category: Category name (e.g., 't_cell', 'myeloid')

    Returns:
        List of reference info dicts
    """
    results = []
    refs_dir = get_references_dir()

    # Check index
    index = load_reference_index()
    if index and 'references' in index:
        for ref_id, ref_info in index['references'].items():
            if ref_info.get('category', '').lower() == category.lower():
                results.append({
                    'id': ref_id,
                    'path': str(refs_dir / ref_info.get('path', '')),
                    'metadata': ref_info
                })

    # Also check filesystem
    cat_dir = refs_dir / category
    if cat_dir.exists():
        for md_file in cat_dir.rglob("*.md"):
            rel_path = md_file.relative_to(refs_dir)
            ref_id = str(rel_path.with_suffix('')).replace(os.sep, '_')

            # Skip if already in results
            if any(r['id'] == ref_id for r in results):
                continue

            content = load_markdown_file(md_file)
            frontmatter, _ = parse_yaml_frontmatter(content) if content else ({}, "")

            results.append({
                'id': ref_id,
                'path': str(md_file),
                'metadata': frontmatter
            })

    return results
