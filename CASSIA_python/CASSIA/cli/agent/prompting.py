"""Prompt construction for the integrated coding-agent controller."""

from __future__ import annotations

import shlex
import sys
from pathlib import Path
from typing import Optional

from .policy import strategy_policy


def default_command_prefix() -> str:
    """Return the current interpreter's installed CASSIA CLI prefix."""
    return shlex.join([sys.executable, "-m", "CASSIA.cli"])


def build_automation_prompt(
    *,
    workdir: Path,
    strategy: str,
    species: str,
    tissue: str,
    additional_info: Optional[str],
    max_topology_edits: int,
    max_agent_commands: int,
    min_child_cells: int = 5,
    max_children_per_split: int = 4,
    command_prefix: Optional[str] = None,
) -> str:
    """Build a reproducible task prompt for the coding-agent controller."""
    if max_agent_commands < 1:
        raise ValueError("max_agent_commands must be at least 1")
    policy = strategy_policy(
        strategy, max_topology_edits, min_child_cells, max_children_per_split
    )
    prefix = command_prefix or default_command_prefix()
    extra = additional_info.strip() if additional_info else "None"
    allowed = ", ".join(policy["allowed_topology_actions"]) or "none"
    system_prompt = (
        Path(__file__).with_name("system_prompt.md").read_text(encoding="utf-8").strip()
    )

    return f"""{system_prompt}

## Automated run contract

You are operating an already initialized, persisted CASSIA run. Complete one
clustering-and-annotation run without editing source code or the input RDS.

Execution context:
- Workdir: {workdir}
- Species: {species}
- Tissue: {tissue}
- Additional dataset context: {extra}
- Strategy: {strategy}
- Strategy meaning: {policy['description']}
- Allowed topology actions: {allowed}
- Enforced topology-edit budget: {policy['max_topology_edits']}
- Minimum cells in every accepted split child: {policy['min_child_cells']}
- Maximum children in one accepted split: {policy['max_children_per_split']}
- Agent command budget: {max_agent_commands}

Use this exact command prefix instead of the shorter `cassia` examples above:

  {prefix} agent

Every stateful command must include:

  --workdir {shlex.quote(str(workdir.parent))}

Required workflow:
1. Run `status`, `clusters`, and `markers all --top 15` to orient and populate
   the evidence cache.
2. Use batched `gene` queries and pairwise `markers A --vs B` only when they
   test a concrete biological ambiguity.
3. Follow the selected strategy. Never run global `preprocess`; it is locked.
   For an evidence-supported split, use `subcluster UUID --auto-resolution`
   so the CLI, not you, searches the conservative resolution ladder and commits
   only the first split that satisfies the child-size/fragmentation policy.
4. After every merge or subcluster, refresh markers for every new alive cluster
   before labeling it.
5. Label every resolved alive cluster with marker-backed evidence. Explicitly
   `--skip` ambiguous, low-quality, or mixed clusters with a concise reason.
6. Run `qa`. Fix actionable errors without `--force`.
7. Run `finalize` only after QA passes. Do not use `finalize --force`.
8. End with a short summary containing the final QA result, labels/skips,
   topology edits, and output paths.

Hard boundaries:
- Use only the CASSIA agent commands above for scientific state changes.
- Do not invoke `cassia agent auto` from inside this run.
- Do not access the network, install packages, publish, commit, or push.
- Do not claim an absent top marker is absent expression; say it was not found
  in the cached differential-marker evidence.
- Stop when the command budget is reached and report what remains unresolved.
"""


__all__ = ["build_automation_prompt", "default_command_prefix"]
