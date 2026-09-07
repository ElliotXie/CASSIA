"""
Single-step LLM reference selector for Reference Agent.

One LLM call that, given markers + optional context, simultaneously:
  - Infers a preliminary cell type / candidate cell types
  - Selects the most relevant reference files from the router
  - Returns a short reasoning string

There is no separate "do we need a reference?" gating step: the caller
opts into reference mode explicitly (use_reference=True), so the agent
always proceeds to selection. If the library has nothing relevant, the
LLM is free to return an empty selected_references list.
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional

try:
    from CASSIA.core.llm_utils import call_llm
except ImportError:
    try:
        from ...core.llm_utils import call_llm
    except ImportError:
        from llm_utils import call_llm


def _load_router_content() -> str:
    """Load the router markdown file content."""
    router_path = Path(__file__).parent / "references_brain" / "_router.md"
    if router_path.exists():
        with open(router_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""


SELECTION_PROMPT = """You are helping annotate a single-cell RNA-seq cluster. Given its top marker genes, infer the likely cell type and pick the most relevant expert reference files from the library.

## Top Marker Genes (ranked by expression):
{markers}

## Context:
- Tissue: {tissue}
- Species: {species}
- Cell type hint (from caller, may be empty): {cell_type_hint}

## Available Reference Library Structure:
{router_content}

## Your Task:
1. Infer the most likely cell type from the markers (one short label, e.g. "Macrophage", "CD8 T cell").
2. List the plausible alternatives the reference should help distinguish between.
3. Pick 1-3 reference files from the library that are most relevant. Prefer specific subtype files over overview files when the markers point to a subtype. If nothing in the library is a good match, return an empty list.

## Response Format (JSON only):
```json
{{
    "preliminary_cell_type": "...",
    "cell_type_range": ["...", "..."],
    "selected_references": ["myeloid/macrophage/tam_spp1.md", "myeloid/macrophage/_overview.md"],
    "reasoning": "brief explanation"
}}
```

Respond ONLY with the JSON object, no additional text."""


def select_references_llm(
    markers: List[str],
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    cell_type_hint: Optional[str] = None,
    provider: str = "openrouter",
    model: Optional[str] = None,
    temperature: float = 0,
    api_key: Optional[str] = None,
) -> Dict:
    """Single LLM call: infer cell type + select references from the router.

    Returns a dict with:
        - preliminary_cell_type: str
        - cell_type_range: List[str]
        - selected_references: List[str]  (paths relative to references_brain/)
        - reasoning: str
    """
    if model is None:
        try:
            from ...core.model_settings import get_model_settings
            model = (
                get_model_settings()
                .settings.get("providers", {})
                .get(provider, {})
                .get("fast", "google/gemini-3.8-flash")
            )
        except Exception:
            model = "google/gemini-3.8-flash"

    router_content = _load_router_content()
    if not router_content:
        return {
            "preliminary_cell_type": "Unknown",
            "cell_type_range": [],
            "selected_references": [],
            "reasoning": "Router file not found",
            "error": "Router not available",
        }

    prompt = SELECTION_PROMPT.format(
        markers=", ".join(markers[:20]),
        tissue=tissue or "Unknown",
        species=species or "Unknown",
        cell_type_hint=cell_type_hint or "",
        router_content=router_content,
    )

    try:
        response = call_llm(
            prompt=prompt,
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=768,
            api_key=api_key,
        )
        return _parse_response(response)
    except Exception as e:
        return {
            "preliminary_cell_type": "Unknown",
            "cell_type_range": [],
            "selected_references": [],
            "reasoning": f"Error during reference selection: {e}",
            "error": str(e),
        }


def _parse_response(response: str) -> Dict:
    """Parse the LLM's JSON response, with a regex fallback."""
    json_match = re.search(r"\{[\s\S]*\}", response)
    if json_match:
        try:
            result = json.loads(json_match.group())
            return {
                "preliminary_cell_type": result.get("preliminary_cell_type", "Unknown"),
                "cell_type_range": result.get("cell_type_range", []) or [],
                "selected_references": _normalize_paths(
                    result.get("selected_references", []) or []
                ),
                "reasoning": result.get("reasoning", ""),
            }
        except json.JSONDecodeError:
            pass
    return _fallback_parse(response)


def _normalize_paths(paths: List[str]) -> List[str]:
    normalized = []
    for path in paths:
        path = (path or "").strip().lstrip("/").lstrip("\\")
        if path.startswith("references/"):
            path = path[len("references/"):]
        if path.startswith("references_brain/"):
            path = path[len("references_brain/"):]
        if path:
            normalized.append(path)
    return normalized


def _fallback_parse(response: str) -> Dict:
    """Best-effort fallback when JSON extraction fails."""
    path_pattern = r"([a-z_]+/[a-z_/]+\.md)"
    matches = re.findall(path_pattern, response.lower())
    return {
        "preliminary_cell_type": "Unknown",
        "cell_type_range": [],
        "selected_references": list(dict.fromkeys(matches))[:3],
        "reasoning": "Parsed from natural-language response (JSON extraction failed)",
    }
