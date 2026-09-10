"""Audit parsing and deterministic comparison for saved automation runs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence


def read_audit_entries(path: Path) -> List[Dict[str, Any]]:
    """Read valid JSON objects from an append-only CASSIA audit log."""
    if not path.exists():
        return []
    entries: List[Dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            try:
                value = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(value, dict):
                entries.append(value)
    return entries


def summarize_automation_audit(
    entries: Iterable[Mapping[str, Any]],
    *,
    start_index: int = 0,
) -> Dict[str, Any]:
    """Return deterministic policy/operation metrics for one agent attempt."""
    values = list(entries)[start_index:]
    ops = [entry for entry in values if entry.get("event") == "op"]
    successful = [entry for entry in ops if entry.get("ok") is True]
    failed = [entry for entry in ops if entry.get("ok") is False]
    committed_ids = {
        str(entry.get("id"))
        for entry in values
        if entry.get("event") == "transaction_commit"
    }
    topology = [
        entry
        for entry in successful
        if entry.get("op") in {"merge", "subcluster"}
        and (
            entry.get("execution_mode") != "one-shot"
            or str(entry.get("id")) in committed_ids
        )
    ]
    return {
        "command_count": len(ops),
        "successful_command_count": len(successful),
        "failed_command_count": len(failed),
        "topology_edit_count": len(topology),
        "topology_operations": [str(entry.get("op")) for entry in topology],
        "operation_counts": {
            name: sum(1 for entry in ops if entry.get("op") == name)
            for name in sorted({str(entry.get("op")) for entry in ops})
        },
    }


def summarize_automation_run(run_dir: Path) -> Dict[str, Any]:
    """Summarize one completed automation directory without rerunning science."""
    run_dir = Path(run_dir).expanduser().resolve()
    manifest_path = run_dir / "automation_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Automation manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    attempts = manifest.get("attempts") or []
    attempt = attempts[-1] if attempts else {}
    audit = attempt.get("audit") or {}
    artifacts = manifest.get("artifacts") or {}
    baseline = artifacts.get("baseline_partition") or {}
    final = artifacts.get("final_partition") or {}
    outputs = artifacts.get("final_outputs") or {}
    qa = attempt.get("qa") or outputs.get("qa") or {}
    policy = manifest.get("policy") or {}

    annotation_path_value = outputs.get("out_tsv")
    local_candidate = run_dir / "outputs" / "annotation.tsv"
    if local_candidate.exists():
        annotation_path: Optional[Path] = local_candidate
    elif annotation_path_value:
        annotation_path = Path(annotation_path_value)
    else:
        annotation_path = None

    cluster_sizes: List[int] = []
    labeled_cells = 0
    skipped_cells = 0
    if annotation_path is not None and annotation_path.exists():
        with annotation_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    n_cells = int(row.get("n_cells") or 0)
                except ValueError:
                    n_cells = 0
                cluster_sizes.append(n_cells)
                if str(row.get("skip", "")).strip().lower() in {"true", "1", "yes"}:
                    skipped_cells += n_cells
                else:
                    labeled_cells += n_cells

    total_cells = labeled_cells + skipped_cells
    if total_cells == 0:
        total_cells = int(final.get("n_cells") or baseline.get("n_cells") or 0)
    # Older manifests predate the split guard; use today's conservative
    # threshold so historical experiments remain comparable.
    min_child_cells = int(policy.get("min_child_cells") or 5)
    final_clusters = int(final.get("n_clusters") or qa.get("n_alive") or len(cluster_sizes))
    baseline_clusters = int(baseline.get("n_clusters") or 0)
    return {
        "name": run_dir.name,
        "run_dir": str(run_dir),
        "status": manifest.get("status"),
        "strategy": manifest.get("strategy"),
        "model": manifest.get("model"),
        "reasoning_effort": manifest.get("reasoning_effort"),
        "qa_pass": bool(qa.get("pass")),
        "commands": int(audit.get("command_count") or 0),
        "failed_commands": int(audit.get("failed_command_count") or 0),
        "topology_edits": int(audit.get("topology_edit_count") or 0),
        "baseline_clusters": baseline_clusters,
        "final_clusters": final_clusters,
        "labeled_clusters": int(qa.get("n_labeled") or outputs.get("n_labeled") or 0),
        "skipped_clusters": int(qa.get("n_skipped") or outputs.get("n_skipped") or 0),
        "total_cells": total_cells,
        "labeled_cells": labeled_cells,
        "skipped_cells": skipped_cells,
        "cell_coverage": round(labeled_cells / total_cells, 6) if total_cells else 0.0,
        "smallest_cluster": min(cluster_sizes) if cluster_sizes else None,
        "tiny_clusters": sum(1 for size in cluster_sizes if size < min_child_cells),
        "fragmentation_ratio": (
            round(final_clusters / baseline_clusters, 6) if baseline_clusters else None
        ),
    }


def render_automation_comparison(run_dirs: Sequence[Path]) -> str:
    """Render a deterministic Markdown comparison of saved runs."""
    summaries = [summarize_automation_run(Path(path)) for path in run_dirs]
    lines = [
        "# CASSIA integrated-agent comparison",
        "",
        "| run | strategy | QA | commands | edits | clusters | labeled clusters | cell coverage | smallest cluster | tiny clusters |",
        "|---|---:|:---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in summaries:
        coverage = f"{100 * item['cell_coverage']:.1f}%"
        cluster_change = f"{item['baseline_clusters']}→{item['final_clusters']}"
        smallest = item["smallest_cluster"] if item["smallest_cluster"] is not None else "—"
        lines.append(
            "| {name} | {strategy} | {qa} | {commands} | {edits} | {clusters} | "
            "{labeled} | {coverage} | {smallest} | {tiny} |".format(
                name=item["name"],
                strategy=item["strategy"],
                qa="PASS" if item["qa_pass"] else "FAIL",
                commands=item["commands"],
                edits=item["topology_edits"],
                clusters=cluster_change,
                labeled=item["labeled_clusters"],
                coverage=coverage,
                smallest=smallest,
                tiny=item["tiny_clusters"],
            )
        )
    lines.extend([
        "",
        "`tiny clusters` means final clusters below the run's configured `min_child_cells`; legacy manifests use 5.",
        "Command and edit counts come from the transaction audit, not from the model's prose.",
        "",
    ])
    return "\n".join(lines)


__all__ = [
    "read_audit_entries",
    "render_automation_comparison",
    "summarize_automation_audit",
    "summarize_automation_run",
]
