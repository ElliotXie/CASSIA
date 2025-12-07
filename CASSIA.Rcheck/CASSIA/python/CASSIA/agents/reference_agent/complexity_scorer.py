"""
Two-Step LLM-Based Complexity Scorer for Reference Agent.

Step 1: Assess marker complexity and decide if reference is needed (no router)
Step 2: If reference needed, select specific references from router (ReAct pattern)
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
    router_path = Path(__file__).parent / "references_brain" / "_router.md"
    if router_path.exists():
        with open(router_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""


# =============================================================================
# STEP 1: Complexity Assessment (Does NOT see router)
# =============================================================================

STEP1_COMPLEXITY_PROMPT = """You are an expert single-cell RNA-seq analyst. Analyze the following marker genes and assess the complexity of cell type annotation.

## Top 20 Marker Genes (ranked by expression):
{markers}

## Context:
- Tissue: {tissue}
- Species: {species}

## Your Task:
1. Identify the most likely cell type based on these markers
2. List possible cell types this could be
3. Rate the annotation complexity/ambiguity from 0-100:
   - 0-30: Clear, unambiguous markers pointing to a specific cell type
   - 31-60: Some ambiguity, multiple subtypes possible
   - 61-80: High complexity, markers suggest multiple lineages or rare cell types
   - 81-100: Very ambiguous, conflicting markers or unusual expression patterns
4. Decide if you would benefit from expert reference documentation to make a more accurate annotation

## Response Format (JSON only):
```json
{{
    "preliminary_cell_type": "the most likely general cell type",
    "cell_type_range": ["list", "of", "possible", "cell", "types"],
    "complexity_score": 0-100,
    "requires_reference": true/false,
    "reasoning": "brief explanation of your assessment"
}}
```

Respond ONLY with the JSON object, no additional text."""


def assess_complexity_step1(
    markers: List[str],
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    provider: str = "openrouter",
    model: Optional[str] = None,
    temperature: float = 0,
    api_key: Optional[str] = None
) -> Dict:
    """
    STEP 1: Assess marker complexity and decide if reference is needed.

    This is the FIRST LLM call. It does NOT see the router yet.
    The LLM decides if it needs reference documentation based on marker complexity.

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
            - requires_reference: bool (LLM's decision)
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

    # Build prompt (NO router in Step 1)
    prompt = STEP1_COMPLEXITY_PROMPT.format(
        markers=markers_str,
        tissue=tissue or "Unknown",
        species=species or "Unknown"
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
        result = _parse_step1_response(response)
        return result

    except Exception as e:
        # Return a safe default on error
        return {
            "complexity_score": 50,
            "preliminary_cell_type": "Unknown",
            "cell_type_range": [],
            "requires_reference": True,  # Default to using reference when uncertain
            "reasoning": f"Error during Step 1 complexity assessment: {str(e)}",
            "error": str(e)
        }


def _parse_step1_response(response: str) -> Dict:
    """Parse the Step 1 LLM response."""
    json_match = re.search(r'\{[\s\S]*\}', response)

    if json_match:
        try:
            result = json.loads(json_match.group())
            return {
                "complexity_score": _normalize_score(result.get("complexity_score", 50)),
                "preliminary_cell_type": result.get("preliminary_cell_type", "Unknown"),
                "cell_type_range": result.get("cell_type_range", []),
                "requires_reference": result.get("requires_reference", True),
                "reasoning": result.get("reasoning", "")
            }
        except json.JSONDecodeError:
            pass

    # Fallback parsing
    return _fallback_parse_step1(response)


def _fallback_parse_step1(response: str) -> Dict:
    """Fallback parsing for Step 1 when JSON extraction fails."""
    response_lower = response.lower()

    # Determine complexity from keywords
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
    ]
    for pattern in cell_type_patterns:
        match = re.search(pattern, response, re.IGNORECASE)
        if match:
            cell_type = match.group(1).strip()
            break

    # Check if reference is needed based on keywords
    requires_reference = (
        complexity_score > 40 or
        "need" in response_lower or
        "reference" in response_lower or
        "uncertain" in response_lower
    )

    return {
        "complexity_score": complexity_score,
        "preliminary_cell_type": cell_type,
        "cell_type_range": [],
        "requires_reference": requires_reference,
        "reasoning": "Parsed from natural language response (JSON extraction failed)"
    }


# =============================================================================
# STEP 2: Reference Selection (NOW sees router)
# =============================================================================

STEP2_SELECTION_PROMPT = """You are selecting expert reference documents to help with cell type annotation.

## Annotation Context:
- Preliminary cell type: {preliminary_cell_type}
- Possible cell types: {cell_type_range}
- Top markers: {markers}
- Tissue: {tissue}
- Species: {species}

## Available Reference Library Structure:
{router_content}

## Your Task:
Select which reference file(s) would be most helpful for annotating this cell type accurately.
Choose 1-3 files that are most relevant. Prefer specific files over overview files when appropriate.

## Response Format (JSON only):
```json
{{
    "selected_references": ["t_cell/cd4/treg.md", "t_cell/cd4/_overview.md"],
    "reasoning": "brief explanation of why these references were selected"
}}
```

Respond ONLY with the JSON object, no additional text."""


def select_references_step2(
    markers: List[str],
    preliminary_cell_type: str,
    cell_type_range: List[str],
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    provider: str = "openrouter",
    model: Optional[str] = None,
    temperature: float = 0,
    api_key: Optional[str] = None
) -> Dict:
    """
    STEP 2: Select references from router.

    This is the SECOND LLM call. It NOW sees the router structure.
    Only called if Step 1 said requires_reference=True.

    Args:
        markers: List of marker genes
        preliminary_cell_type: Best guess from Step 1
        cell_type_range: Possible cell types from Step 1
        tissue: Optional tissue type
        species: Optional species
        provider: LLM provider
        model: Specific model to use
        temperature: LLM temperature
        api_key: Optional API key

    Returns:
        Dict with:
            - selected_references: List[str] (file paths)
            - reasoning: str
    """
    # Use fast model by default
    if model is None:
        default_models = {
            "openrouter": "google/gemini-2.5-flash",
            "openai": "gpt-4o-mini",
            "anthropic": "claude-3-haiku-20240307"
        }
        model = default_models.get(provider, "google/gemini-2.5-flash")

    # Load router content
    router_content = _load_router_content()
    if not router_content:
        return {
            "selected_references": [],
            "reasoning": "Router file not found",
            "error": "Router not available"
        }

    # Format inputs
    markers_str = ", ".join(markers[:20])
    cell_type_range_str = ", ".join(cell_type_range) if cell_type_range else "Unknown"

    # Build prompt (NOW includes router)
    prompt = STEP2_SELECTION_PROMPT.format(
        markers=markers_str,
        preliminary_cell_type=preliminary_cell_type or "Unknown",
        cell_type_range=cell_type_range_str,
        tissue=tissue or "Unknown",
        species=species or "Unknown",
        router_content=router_content
    )

    try:
        # Call LLM
        response = call_llm(
            prompt=prompt,
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=512,
            api_key=api_key
        )

        # Parse JSON response
        result = _parse_step2_response(response)
        return result

    except Exception as e:
        return {
            "selected_references": [],
            "reasoning": f"Error during Step 2 reference selection: {str(e)}",
            "error": str(e)
        }


def _parse_step2_response(response: str) -> Dict:
    """Parse the Step 2 LLM response."""
    json_match = re.search(r'\{[\s\S]*\}', response)

    if json_match:
        try:
            result = json.loads(json_match.group())
            selected = result.get("selected_references", [])

            # Normalize paths (remove leading slashes, ensure consistent format)
            normalized = []
            for path in selected:
                path = path.strip().lstrip("/").lstrip("\\")
                # Remove "references/" prefix if present
                if path.startswith("references/"):
                    path = path[11:]
                normalized.append(path)

            return {
                "selected_references": normalized,
                "reasoning": result.get("reasoning", "")
            }
        except json.JSONDecodeError:
            pass

    # Fallback: try to extract paths from response
    return _fallback_parse_step2(response)


def _fallback_parse_step2(response: str) -> Dict:
    """Fallback parsing for Step 2 when JSON extraction fails."""
    # Try to extract file paths from response
    path_pattern = r'([a-z_]+/[a-z_/]+\.md)'
    matches = re.findall(path_pattern, response.lower())

    return {
        "selected_references": list(set(matches))[:3],  # Dedupe and limit to 3
        "reasoning": "Parsed from natural language response (JSON extraction failed)"
    }


# =============================================================================
# Main Entry Point (Two-Step Orchestration)
# =============================================================================

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
    Two-step complexity assessment with ReAct-style reference selection.

    Step 1: Assess complexity (no router)
    Step 2: Select references (with router) - only if Step 1 says requires_reference=True

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
            - preliminary_cell_type: str
            - cell_type_range: List[str]
            - requires_reference: bool (from Step 1)
            - selected_references: List[str] (from Step 2, empty if not required)
            - reasoning: str (combined reasoning)
    """
    # STEP 1: Assess complexity
    step1_result = assess_complexity_step1(
        markers=markers,
        tissue=tissue,
        species=species,
        provider=provider,
        model=model,
        temperature=temperature,
        api_key=api_key
    )

    # Build combined result
    result = {
        "complexity_score": step1_result.get("complexity_score", 50),
        "preliminary_cell_type": step1_result.get("preliminary_cell_type", "Unknown"),
        "cell_type_range": step1_result.get("cell_type_range", []),
        "requires_reference": step1_result.get("requires_reference", False),
        "selected_references": [],
        "reasoning": step1_result.get("reasoning", ""),
        # Legacy fields for backward compatibility
        "needs_reference": step1_result.get("requires_reference", False),
        "reference_categories": [],
        "specific_reference_paths": []
    }

    # STEP 2: Select references (only if Step 1 says it's needed)
    if step1_result.get("requires_reference", False):
        step2_result = select_references_step2(
            markers=markers,
            preliminary_cell_type=step1_result.get("preliminary_cell_type"),
            cell_type_range=step1_result.get("cell_type_range", []),
            tissue=tissue,
            species=species,
            provider=provider,
            model=model,
            temperature=temperature,
            api_key=api_key
        )

        result["selected_references"] = step2_result.get("selected_references", [])
        result["specific_reference_paths"] = step2_result.get("selected_references", [])
        result["reasoning"] += f" | Reference selection: {step2_result.get('reasoning', '')}"

    return result


# =============================================================================
# Utility Functions
# =============================================================================

def _normalize_score(score) -> int:
    """Normalize complexity score to 0-100 range."""
    try:
        score = int(score)
        return max(0, min(100, score))
    except (TypeError, ValueError):
        return 50


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
