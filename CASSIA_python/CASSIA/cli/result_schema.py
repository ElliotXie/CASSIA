"""Stable, deterministic result normalization for CASSIA CLI workflows.

The LLM-facing prompts intentionally remain workflow-specific.  This module is
the compatibility boundary after an LLM response has been parsed: every public
annotation artifact receives the same core fields without requiring another
LLM formatting call.
"""

from __future__ import annotations

import json
from typing import Any, Dict, Iterable, List, Optional, Sequence


ANNOTATION_SCHEMA_VERSION = "cassia.annotation.v1"


def _first_present(payload: Dict[str, Any], keys: Sequence[str]) -> Any:
    for key in keys:
        value = payload.get(key)
        if value is None:
            continue
        if isinstance(value, str) and not value.strip():
            continue
        if isinstance(value, (list, tuple, dict, set)) and not value:
            continue
        return value
    return None


def _text(value: Any, field: str, required: bool = False) -> str:
    if value is None:
        value = ""
    if isinstance(value, str):
        normalized = value.strip()
    elif isinstance(value, (dict, list, tuple)):
        normalized = json.dumps(value, ensure_ascii=False)
    else:
        normalized = str(value).strip()
    if required and not normalized:
        raise ValueError(f"Annotation JSON is missing required field '{field}'")
    return normalized


def _string_list(value: Any) -> List[str]:
    if value is None or value == "":
        return []
    if isinstance(value, str):
        candidates: Iterable[Any] = value.split(",")
    elif isinstance(value, (list, tuple, set)):
        candidates = value
    else:
        candidates = [value]

    normalized: List[str] = []
    seen = set()
    for item in candidates:
        text = _text(item, "list item")
        if text and text.casefold() not in seen:
            normalized.append(text)
            seen.add(text.casefold())
    return normalized


def _normalize_changed(value: Any) -> Optional[bool]:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"null", "none", "n/a", "na", "not applicable"}:
            return None
        if lowered in {"true", "yes", "y", "1"}:
            return True
        if lowered in {"false", "no", "n", "0"}:
            return False
    raise ValueError("Annotation field 'changed_from_original' must be a boolean")


def normalize_annotation_payload(
    payload: Dict[str, Any],
    *,
    cluster_id: Optional[str] = None,
    markers: Optional[Sequence[str]] = None,
    annotation_mode: str = "standard",
) -> Dict[str, Any]:
    """Return a canonical CASSIA annotation while preserving extra fields."""
    if not isinstance(payload, dict):
        raise ValueError("Annotation result must be a JSON object")

    normalized = dict(payload)
    main_cell_type = _text(
        _first_present(normalized, ("main_cell_type", "final_cell_type")),
        "main_cell_type",
        required=True,
    )
    sub_cell_types = _string_list(
        _first_present(
            normalized,
            ("sub_cell_types", "ranked_sub_cell_types", "final_sub_cell_type"),
        )
    )
    possible_mixed = _string_list(normalized.get("possible_mixed_cell_types"))

    normalized.update({
        "schema_version": ANNOTATION_SCHEMA_VERSION,
        "annotation_mode": annotation_mode,
        "main_cell_type": main_cell_type,
        "sub_cell_types": sub_cell_types,
        "possible_mixed_cell_types": possible_mixed,
        "confidence": (
            normalized.get("confidence")
            if isinstance(normalized.get("confidence"), (int, float))
            and not isinstance(normalized.get("confidence"), bool)
            else _text(normalized.get("confidence"), "confidence")
        ),
        "evidence": _text(normalized.get("evidence"), "evidence"),
    })
    if cluster_id is not None:
        normalized["cluster_id"] = str(cluster_id)
    if markers is not None:
        marker_list = _string_list(list(markers))
        normalized["num_markers"] = len(marker_list)
        normalized["marker_list"] = marker_list
    return normalized


def normalize_fused_boost_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize Fused Boost output and expose the canonical annotation fields."""
    normalized = normalize_annotation_payload(payload, annotation_mode="fused_boost")
    normalized["final_cell_type"] = normalized["main_cell_type"]

    final_sub_cell_type = _text(
        _first_present(normalized, ("final_sub_cell_type",)),
        "final_sub_cell_type",
    )
    ranked = _string_list(normalized.get("ranked_sub_cell_types"))
    if final_sub_cell_type:
        ranked = [final_sub_cell_type] + [
            item for item in ranked if item.casefold() != final_sub_cell_type.casefold()
        ]
    elif ranked:
        final_sub_cell_type = ranked[0]
    elif normalized["sub_cell_types"]:
        final_sub_cell_type = normalized["sub_cell_types"][0]
        ranked = list(normalized["sub_cell_types"])

    normalized["final_sub_cell_type"] = final_sub_cell_type
    normalized["ranked_sub_cell_types"] = ranked
    normalized["sub_cell_types"] = ranked or (
        [final_sub_cell_type] if final_sub_cell_type else []
    )
    raw_changed = normalized.get("changed_from_original")
    try:
        normalized["changed_from_original"] = _normalize_changed(raw_changed)
    except ValueError:
        # Direct Fused Boost has no trusted prior annotation. Some models put
        # the name of the displaced preliminary candidate in this legacy field
        # even though the prompt requests null. Preserve that useful audit text
        # separately while keeping the public field schema-correct.
        normalized["changed_from_candidate"] = _text(
            raw_changed,
            "changed_from_candidate",
        )
        normalized["changed_from_original"] = None
    for key in (
        "checked_genes",
        "supporting_markers",
        "refuting_markers",
        "alternatives",
    ):
        normalized[key] = _string_list(normalized.get(key))
    normalized["recommended_next_steps"] = _text(
        normalized.get("recommended_next_steps"),
        "recommended_next_steps",
    )
    return normalized


__all__ = [
    "ANNOTATION_SCHEMA_VERSION",
    "normalize_annotation_payload",
    "normalize_fused_boost_payload",
]
