"""Topology policies for integrated clustering and annotation runs."""

from __future__ import annotations

from typing import Any, Dict


AUTOMATION_STRATEGIES = ("fixed", "conservative", "adaptive")
DEFAULT_MIN_CHILD_CELLS = 5
DEFAULT_MAX_CHILDREN_PER_SPLIT = 4


def strategy_policy(
    strategy: str,
    max_topology_edits: int,
    min_child_cells: int = DEFAULT_MIN_CHILD_CELLS,
    max_children_per_split: int = DEFAULT_MAX_CHILDREN_PER_SPLIT,
) -> Dict[str, Any]:
    """Validate user limits and resolve the selected strategy contract."""
    if strategy not in AUTOMATION_STRATEGIES:
        raise ValueError(
            f"Unknown automation strategy '{strategy}'. Choose: "
            f"{', '.join(AUTOMATION_STRATEGIES)}"
        )
    if max_topology_edits < 0:
        raise ValueError("max_topology_edits must be non-negative")
    if min_child_cells < 1:
        raise ValueError("min_child_cells must be at least 1")
    if max_children_per_split < 2:
        raise ValueError("max_children_per_split must be at least 2")

    if strategy == "fixed":
        return {
            "allowed_topology_actions": [],
            "max_topology_edits": 0,
            "min_child_cells": min_child_cells,
            "max_children_per_split": max_children_per_split,
            "description": "Keep the supplied partition fixed; investigate and label only.",
        }
    if strategy == "conservative":
        return {
            "allowed_topology_actions": ["merge", "subcluster"],
            "max_topology_edits": min(max_topology_edits, 2),
            "min_child_cells": min_child_cells,
            "max_children_per_split": min(max_children_per_split, 3),
            "description": (
                "Permit at most two strongly evidenced local topology edits; "
                "prefer leaving an ambiguous cluster unresolved."
            ),
        }
    return {
        "allowed_topology_actions": ["merge", "subcluster"],
        "max_topology_edits": max_topology_edits,
        "min_child_cells": min_child_cells,
        "max_children_per_split": max_children_per_split,
        "description": (
            "Permit bounded evidence-led local refinement while preserving the "
            "global partition and exact cell membership provenance."
        ),
    }


__all__ = [
    "AUTOMATION_STRATEGIES",
    "DEFAULT_MAX_CHILDREN_PER_SPLIT",
    "DEFAULT_MIN_CHILD_CELLS",
    "strategy_policy",
]
