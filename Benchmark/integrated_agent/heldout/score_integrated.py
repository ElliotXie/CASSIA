#!/usr/bin/env python3
"""Score a completed CASSIA integrated-agent run against held-out cell truth.

This module is intentionally benchmark-only. The agent run is completed before
this program receives the truth file, so scoring fields never need to be in the
Seurat object or automation prompt.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SKIP_VALUES = {"1", "true", "yes", "y"}


def _read_csv(path: Path, *, delimiter: str = ",") -> List[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def _require_columns(rows: Sequence[Mapping[str, str]], columns: Iterable[str], name: str) -> None:
    available = set(rows[0]) if rows else set()
    missing = sorted(set(columns) - available)
    if missing:
        raise ValueError(f"{name} is missing required columns: {', '.join(missing)}")


def _index_unique(rows: Sequence[Mapping[str, str]], key: str, name: str) -> Dict[str, Mapping[str, str]]:
    result: Dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = str(row.get(key, "")).strip()
        if not value:
            raise ValueError(f"{name} contains an empty {key}")
        if value in result:
            raise ValueError(f"{name} contains duplicate {key}: {value}")
        result[value] = row
    return result


def normalize_label(value: str) -> str:
    """Conservative normalization for exact scoring without an adjudicated map."""
    return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()


def _is_within(path: Path, directory: Path) -> bool:
    path = path.resolve()
    directory = directory.resolve()
    return path == directory or directory in path.parents


def _comb2(value: int) -> float:
    return value * (value - 1) / 2


def clustering_metrics(truth: Sequence[str], predicted: Sequence[str]) -> Dict[str, float]:
    """Return dependency-free ARI/NMI and information-theoretic metrics."""
    if len(truth) != len(predicted) or not truth:
        raise ValueError("truth and predicted clusters must have equal non-zero length")

    n = len(truth)
    truth_counts = Counter(truth)
    predicted_counts = Counter(predicted)
    contingency = Counter(zip(truth, predicted))
    pair_total = _comb2(n)
    sum_cells = sum(_comb2(count) for count in contingency.values())
    sum_truth = sum(_comb2(count) for count in truth_counts.values())
    sum_predicted = sum(_comb2(count) for count in predicted_counts.values())
    expected = (sum_truth * sum_predicted / pair_total) if pair_total else 0.0
    denominator = 0.5 * (sum_truth + sum_predicted) - expected
    ari = (sum_cells - expected) / denominator if denominator else 1.0

    def entropy(counts: Iterable[int]) -> float:
        return -sum((count / n) * math.log(count / n) for count in counts if count)

    truth_entropy = entropy(truth_counts.values())
    predicted_entropy = entropy(predicted_counts.values())
    mutual_information = sum(
        (count / n) * math.log((count * n) / (truth_counts[t] * predicted_counts[p]))
        for (t, p), count in contingency.items()
        if count
    )
    nmi_denominator = 0.5 * (truth_entropy + predicted_entropy)
    nmi = mutual_information / nmi_denominator if nmi_denominator else 1.0
    homogeneity = mutual_information / truth_entropy if truth_entropy else 1.0
    completeness = mutual_information / predicted_entropy if predicted_entropy else 1.0
    v_measure = (
        2 * homogeneity * completeness / (homogeneity + completeness)
        if homogeneity + completeness
        else 0.0
    )
    purity = sum(
        max((contingency.get((truth_label, cluster), 0) for truth_label in truth_counts), default=0)
        for cluster in predicted_counts
    ) / n

    return {
        "adjusted_rand_index": ari,
        "normalized_mutual_information": nmi,
        "homogeneity": homogeneity,
        "completeness": completeness,
        "v_measure": v_measure,
        "purity": purity,
    }


def _classification_metrics(truth: Sequence[str], predicted: Sequence[str]) -> Dict[str, float]:
    if len(truth) != len(predicted) or not truth:
        raise ValueError("truth and predicted labels must have equal non-zero length")
    labels = sorted(set(truth) | set(predicted))
    f1_values: List[float] = []
    recalls: List[float] = []
    for label in labels:
        tp = sum(t == label and p == label for t, p in zip(truth, predicted))
        fp = sum(t != label and p == label for t, p in zip(truth, predicted))
        fn = sum(t == label and p != label for t, p in zip(truth, predicted))
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1_values.append(2 * precision * recall / (precision + recall) if precision + recall else 0.0)
        if any(t == label for t in truth):
            recalls.append(recall)
    return {
        "accuracy": sum(t == p for t, p in zip(truth, predicted)) / len(truth),
        "macro_f1": sum(f1_values) / len(f1_values),
        "balanced_accuracy": sum(recalls) / len(recalls),
    }


def _load_label_map(path: Optional[Path]) -> Dict[str, Tuple[str, str]]:
    if path is None:
        return {}
    rows = _read_csv(path)
    _require_columns(
        rows,
        ["predicted_label", "canonical_fine_label", "canonical_broad_label"],
        "label map",
    )
    result: Dict[str, Tuple[str, str]] = {}
    for row in rows:
        key = normalize_label(row["predicted_label"])
        if not key or key in result:
            raise ValueError(f"label map has an empty or duplicate predicted_label: {row['predicted_label']}")
        result[key] = (
            normalize_label(row["canonical_fine_label"]),
            normalize_label(row["canonical_broad_label"]),
        )
    return result


def score_run(
    *,
    truth_path: Path,
    run_dir: Path,
    label_map_path: Optional[Path] = None,
    min_cluster_cells: int = 5,
) -> Dict[str, object]:
    """Score one run using its exported memberships and annotation table."""
    if _is_within(truth_path, run_dir):
        raise ValueError("held-out truth must be outside the agent run directory")
    if label_map_path is not None and _is_within(label_map_path, run_dir):
        raise ValueError("label map must be outside the agent run directory")
    truth_rows = _read_csv(truth_path)
    membership_rows = _read_csv(run_dir / "provenance" / "final_memberships.csv")
    annotation_rows = _read_csv(run_dir / "outputs" / "annotation.tsv", delimiter="\t")
    _require_columns(truth_rows, ["cell_id", "truth_fine_label"], "truth")
    _require_columns(membership_rows, ["cell_id", "cluster_id"], "memberships")
    _require_columns(annotation_rows, ["cluster_uuid", "label", "skip"], "annotations")

    truth_by_cell = _index_unique(truth_rows, "cell_id", "truth")
    membership_by_cell = _index_unique(membership_rows, "cell_id", "memberships")
    annotations = _index_unique(annotation_rows, "cluster_uuid", "annotations")
    truth_cells = set(truth_by_cell)
    membership_cells = set(membership_by_cell)
    if truth_cells != membership_cells:
        missing = sorted(truth_cells - membership_cells)[:5]
        extra = sorted(membership_cells - truth_cells)[:5]
        raise ValueError(
            "truth/membership cell IDs differ; "
            f"missing_from_run={missing}, missing_from_truth={extra}"
        )

    cells = sorted(truth_cells)
    truth_fine = [normalize_label(truth_by_cell[cell]["truth_fine_label"]) for cell in cells]
    has_broad = "truth_broad_label" in truth_rows[0]
    truth_broad = [normalize_label(truth_by_cell[cell].get("truth_broad_label", "")) for cell in cells]
    predicted_clusters = [membership_by_cell[cell]["cluster_id"] for cell in cells]
    topology = clustering_metrics(truth_fine, predicted_clusters)
    topology_broad = clustering_metrics(truth_broad, predicted_clusters) if has_broad else None

    label_map = _load_label_map(label_map_path)
    prediction_fine: List[str] = []
    prediction_broad: List[str] = []
    covered: List[bool] = []
    mapped: List[bool] = []
    for cell in cells:
        cluster = membership_by_cell[cell]["cluster_id"]
        annotation = annotations.get(cluster)
        is_covered = bool(annotation) and str(annotation.get("skip", "")).strip().casefold() not in SKIP_VALUES
        normalized = (
            normalize_label(str(annotation.get("label", "")))
            if annotation and is_covered
            else ""
        )
        canonical = label_map.get(normalized) if normalized else None
        covered.append(is_covered)
        mapped.append(bool(canonical) if label_map else bool(normalized))
        prediction_fine.append(canonical[0] if canonical else (normalized or "__unlabeled__"))
        prediction_broad.append(canonical[1] if canonical else (normalized or "__unlabeled__"))

    fine = _classification_metrics(truth_fine, prediction_fine)
    broad = _classification_metrics(truth_broad, prediction_broad) if has_broad else None
    cluster_sizes = Counter(predicted_clusters)
    tiny_clusters = [cluster for cluster, size in cluster_sizes.items() if size < min_cluster_cells]
    tiny_cells = sum(cluster_sizes[cluster] for cluster in tiny_clusters)
    n_cells = len(cells)
    covered_count = sum(covered)
    execution: Optional[Dict[str, object]] = None
    manifest_path = run_dir / "automation_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        attempts = manifest.get("attempts") or []
        latest = attempts[-1] if attempts else {}
        audit = latest.get("audit") or {}
        execution = {
            "strategy": manifest.get("strategy"),
            "model": manifest.get("model"),
            "reasoning_effort": manifest.get("reasoning_effort"),
            "command_count": int(audit.get("command_count") or 0),
            "failed_command_count": int(audit.get("failed_command_count") or 0),
            "topology_edit_count": int(audit.get("topology_edit_count") or 0),
        }

    return {
        "schema_version": "cassia.integrated-heldout.v1",
        "run_dir": str(run_dir.resolve()),
        "truth_path": str(truth_path.resolve()),
        "n_cells": n_cells,
        "n_truth_fine_labels": len(set(truth_fine)),
        "n_predicted_clusters": len(cluster_sizes),
        "topology": topology,
        "topology_broad": topology_broad,
        "annotation": {
            "fine": fine,
            "broad": broad,
            "cell_coverage": covered_count / n_cells,
            "mapped_cell_fraction": sum(mapped) / n_cells,
            "uncovered_cells": n_cells - covered_count,
            "scoring_mode": "adjudicated_label_map" if label_map else "normalized_exact_label",
        },
        "fragmentation": {
            "cluster_count_ratio": len(cluster_sizes) / len(set(truth_fine)),
            "smallest_cluster_cells": min(cluster_sizes.values()),
            "tiny_cluster_count": len(tiny_clusters),
            "tiny_cell_fraction": tiny_cells / n_cells,
            "min_cluster_cells": min_cluster_cells,
        },
        "execution": execution,
    }


def render_markdown(result: Mapping[str, object]) -> str:
    topology = result["topology"]
    annotation = result["annotation"]
    fragmentation = result["fragmentation"]
    assert isinstance(topology, Mapping)
    assert isinstance(annotation, Mapping)
    assert isinstance(fragmentation, Mapping)
    fine = annotation["fine"]
    broad = annotation.get("broad")
    assert isinstance(fine, Mapping)
    lines = [
        "# Integrated held-out evaluation",
        "",
        f"- Cells: {result['n_cells']}",
        f"- Truth labels / predicted clusters: {result['n_truth_fine_labels']} / {result['n_predicted_clusters']}",
        f"- ARI / NMI: {topology['adjusted_rand_index']:.3f} / {topology['normalized_mutual_information']:.3f}",
        f"- Homogeneity / completeness: {topology['homogeneity']:.3f} / {topology['completeness']:.3f}",
        f"- Fine accuracy / macro-F1: {fine['accuracy']:.3f} / {fine['macro_f1']:.3f}",
        f"- Cell coverage: {annotation['cell_coverage']:.3f}",
        f"- Tiny clusters / tiny-cell fraction: {fragmentation['tiny_cluster_count']} / {fragmentation['tiny_cell_fraction']:.3f}",
    ]
    if isinstance(broad, Mapping):
        lines.insert(-2, f"- Broad accuracy / macro-F1: {broad['accuracy']:.3f} / {broad['macro_f1']:.3f}")
    lines.extend(["", f"Scoring mode: `{annotation['scoring_mode']}`.", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--label-map", type=Path)
    parser.add_argument("--min-cluster-cells", type=int, default=5)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.min_cluster_cells < 1:
        parser.error("--min-cluster-cells must be at least 1")
    result = score_run(
        truth_path=args.truth,
        run_dir=args.run,
        label_map_path=args.label_map,
        min_cluster_cells=args.min_cluster_cells,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "metrics.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (args.out / "report.md").write_text(render_markdown(result), encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
