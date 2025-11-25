"""
LLM-Based Complexity Scorer for Reference Agent.

Uses an LLM to assess marker gene complexity and decide whether
reference retrieval is needed for cell type annotation.
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional

# Handle both package and direct imports
try:
    from ..llm_utils import call_llm
except ImportError:
    try:
        from CASSIA.llm_utils import call_llm
    except ImportError:
        # Fallback for direct script execution
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from llm_utils import call_llm


def _load_router_content() -> str:
    """Load the router markdown file content."""
    router_path = Path(__file__).parent / "references" / "_router.md"
    if router_path.exists():
        with open(router_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""


def _get_router_quick_reference() -> str:
    """Extract just the quick reference table from the router."""
    content = _load_router_content()
    if not content:
        return ""

    # Extract the Quick Reference Guide section
    start_marker = "## Quick Reference Guide"
    end_marker = "## Decision Tree"

    if start_marker in content:
        start_idx = content.find(start_marker)
        end_idx = content.find(end_marker) if end_marker in content else len(content)
        return content[start_idx:end_idx].strip()

    return ""


COMPLEXITY_ASSESSMENT_PROMPT = """You are an expert single-cell RNA-seq analyst. Analyze the following marker genes and assess the complexity of cell type annotation.

## Top 20 Marker Genes (ranked by expression):
{markers}

## Context:
- Tissue: {tissue}
- Species: {species}

## Available Reference Categories:
Use this guide to determine which reference folder would be helpful:

{router_guide}

## Your Task:
1. Identify the most likely cell type range (what general cell type(s) these markers suggest)
2. Rate the annotation complexity/ambiguity from 0-100:
   - 0-30: Clear, unambiguous markers pointing to a specific cell type
   - 31-60: Some ambiguity, multiple subtypes possible
   - 61-80: High complexity, markers suggest multiple lineages or rare cell types
   - 81-100: Very ambiguous, conflicting markers or unusual expression patterns
3. Determine if reference documentation would help (set needs_reference=true if complexity > 40 or subtypes need clarification)
4. Suggest which reference categories would be most useful based on the router guide above

## Response Format (JSON only):
```json
{{
    "preliminary_cell_type": "the most likely general cell type",
    "cell_type_range": ["list", "of", "possible", "cell", "types"],
    "complexity_score": 0-100,
    "needs_reference": true/false,
    "reference_categories": ["t_cell", "b_cell", "myeloid", etc.],
    "specific_reference_paths": ["t_cell/cd4/_overview.md", etc.],
    "reasoning": "brief explanation of your assessment"
}}
```

Respond ONLY with the JSON object, no additional text."""


def assess_complexity(
    markers: List[str],
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    provider: str = "openrouter",
    model: Optional[str] = None,
    temperature: float = 0,
    api_key: Optional[str] = None
) -> Dict:
    """
    Use LLM to assess marker complexity and decide if reference is needed.

    Args:
        markers: List of marker genes (top 20 recommended)
        tissue: Optional tissue type for context
        species: Optional species for context
        provider: LLM provider ("openai", "anthropic", "openrouter")
        model: Specific model to use (defaults to fast model)
        temperature: LLM temperature (0 for deterministic)
        api_key: Optional API key

    Returns:
        Dict with:
            - complexity_score: 0-100
            - preliminary_cell_type: str (best guess)
            - cell_type_range: List[str] (possible cell types)
            - needs_reference: bool
            - reference_categories: List[str] (which reference folders to check)
            - reasoning: str
    """
    # Use fast model by default for cost efficiency
    if model is None:
        default_models = {
            "openrouter": "google/gemini-2.5-flash",
            "openai": "gpt-4o-mini",
            "anthropic": "claude-3-haiku-20240307"
        }
        model = default_models.get(provider, "google/gemini-2.5-flash")

    # Format markers
    markers_str = ", ".join(markers[:20])

    # Get router guide for reference categories
    router_guide = _get_router_quick_reference()
    if not router_guide:
        # Fallback if router not available
        router_guide = """| Key Markers | Reference Path |
|-------------|----------------|
| CD3D, CD3E, CD4 | t_cell/cd4/ |
| CD3D, CD3E, CD8A | t_cell/cd8/ |
| CD19, CD79A | b_cell/ |
| CD14, CD68 | myeloid/ |
| NCAM1, NKG7 | nk_cell/ |"""

    # Build prompt
    prompt = COMPLEXITY_ASSESSMENT_PROMPT.format(
        markers=markers_str,
        tissue=tissue or "Unknown",
        species=species or "Unknown",
        router_guide=router_guide
    )

    try:
        # Call LLM
        response = call_llm(
            prompt=prompt,
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=1024,
            api_key=api_key
        )

        # Parse JSON response
        result = _parse_complexity_response(response)
        return result

    except Exception as e:
        # Return a safe default on error
        return {
            "complexity_score": 50,
            "preliminary_cell_type": "Unknown",
            "cell_type_range": [],
            "needs_reference": True,  # Default to using reference when uncertain
            "reference_categories": [],
            "specific_reference_paths": [],
            "reasoning": f"Error during complexity assessment: {str(e)}",
            "error": str(e)
        }


def _parse_complexity_response(response: str) -> Dict:
    """
    Parse the LLM response into a structured dictionary.

    Args:
        response: Raw LLM response string

    Returns:
        Parsed dictionary with complexity assessment
    """
    # Try to extract JSON from response
    json_match = re.search(r'\{[\s\S]*\}', response)

    if json_match:
        try:
            result = json.loads(json_match.group())

            # Validate and normalize fields
            return {
                "complexity_score": _normalize_score(result.get("complexity_score", 50)),
                "preliminary_cell_type": result.get("preliminary_cell_type", "Unknown"),
                "cell_type_range": result.get("cell_type_range", []),
                "needs_reference": result.get("needs_reference", True),
                "reference_categories": _normalize_categories(result.get("reference_categories", [])),
                "specific_reference_paths": result.get("specific_reference_paths", []),
                "reasoning": result.get("reasoning", "")
            }
        except json.JSONDecodeError:
            pass

    # Fallback parsing if JSON extraction fails
    return _fallback_parse(response)


def _normalize_score(score) -> int:
    """Normalize complexity score to 0-100 range."""
    try:
        score = int(score)
        return max(0, min(100, score))
    except (TypeError, ValueError):
        return 50


def _normalize_categories(categories) -> List[str]:
    """Normalize reference category names."""
    if not categories:
        return []

    normalized = []
    category_mapping = {
        "t cell": "t_cell",
        "t-cell": "t_cell",
        "tcell": "t_cell",
        "b cell": "b_cell",
        "b-cell": "b_cell",
        "bcell": "b_cell",
        "myeloid": "myeloid",
        "monocyte": "myeloid",
        "macrophage": "myeloid",
        "dendritic": "myeloid",
        "nk cell": "nk_cell",
        "nk": "nk_cell",
        "natural killer": "nk_cell",
        "epithelial": "epithelial",
        "endothelial": "endothelial",
        "stromal": "stromal",
        "fibroblast": "stromal",
    }

    for cat in categories:
        cat_lower = cat.lower().strip()
        if cat_lower in category_mapping:
            normalized.append(category_mapping[cat_lower])
        else:
            # Use as-is if not in mapping (with underscore normalization)
            normalized.append(cat_lower.replace(" ", "_").replace("-", "_"))

    return list(set(normalized))  # Remove duplicates


def _fallback_parse(response: str) -> Dict:
    """
    Fallback parsing when JSON extraction fails.
    Try to extract information from natural language response.
    """
    response_lower = response.lower()

    # Try to determine complexity from keywords
    complexity_score = 50
    if "clear" in response_lower or "unambiguous" in response_lower:
        complexity_score = 25
    elif "ambiguous" in response_lower or "complex" in response_lower:
        complexity_score = 70
    elif "very ambiguous" in response_lower or "conflicting" in response_lower:
        complexity_score = 85

    # Try to extract cell type
    cell_type = "Unknown"
    cell_type_patterns = [
        r"(?:likely|probably|appears to be)\s+(?:a\s+)?([A-Za-z\s]+cell)",
        r"([A-Za-z]+\s+cell)s?\s+(?:are|is)",
        r"suggests?\s+([A-Za-z\s]+cell)",
    ]
    for pattern in cell_type_patterns:
        match = re.search(pattern, response, re.IGNORECASE)
        if match:
            cell_type = match.group(1).strip()
            break

    # Determine reference categories from content
    categories = []
    if "t cell" in response_lower or "cd4" in response_lower or "cd8" in response_lower:
        categories.append("t_cell")
    if "b cell" in response_lower or "plasma" in response_lower:
        categories.append("b_cell")
    if "monocyte" in response_lower or "macrophage" in response_lower or "myeloid" in response_lower:
        categories.append("myeloid")

    return {
        "complexity_score": complexity_score,
        "preliminary_cell_type": cell_type,
        "cell_type_range": [],
        "needs_reference": complexity_score > 40,
        "reference_categories": categories,
        "specific_reference_paths": [],
        "reasoning": "Parsed from natural language response (JSON extraction failed)"
    }


def quick_complexity_check(markers: List[str]) -> Dict:
    """
    Quick rule-based complexity check without LLM call.
    Useful for fast pre-screening before LLM assessment.

    Args:
        markers: List of marker genes

    Returns:
        Dict with basic complexity indicators
    """
    # Known clear cell type markers
    clear_markers = {
        "t_cell": {"CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"},
        "b_cell": {"CD19", "CD79A", "CD79B", "MS4A1", "PAX5"},
        "myeloid": {"CD14", "CD68", "CSF1R", "ITGAM", "CD33"},
        "nk_cell": {"NCAM1", "NKG7", "KLRB1", "KLRD1", "GNLY"},
        "epithelial": {"EPCAM", "KRT18", "KRT19", "CDH1"},
        "endothelial": {"PECAM1", "VWF", "CDH5", "KDR"},
    }

    markers_upper = {m.upper() for m in markers[:20]}

    # Count overlaps with each category
    category_scores = {}
    for category, cat_markers in clear_markers.items():
        overlap = len(markers_upper & cat_markers)
        if overlap > 0:
            category_scores[category] = overlap

    # Determine complexity
    if not category_scores:
        return {
            "likely_complex": True,
            "suggested_categories": [],
            "reason": "No clear lineage markers detected"
        }

    # Sort by score
    sorted_categories = sorted(category_scores.items(), key=lambda x: -x[1])
    top_category, top_score = sorted_categories[0]

    # Check for ambiguity (multiple categories with similar scores)
    is_ambiguous = len(sorted_categories) > 1 and sorted_categories[1][1] >= top_score * 0.5

    return {
        "likely_complex": is_ambiguous or top_score < 2,
        "suggested_categories": [cat for cat, _ in sorted_categories],
        "top_category": top_category,
        "top_score": top_score,
        "reason": "Multiple lineages detected" if is_ambiguous else "Clear lineage markers"
    }
