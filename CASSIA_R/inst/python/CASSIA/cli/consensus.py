"""Deterministic consensus utilities for CASSIA CLI outputs."""

from __future__ import annotations

import ast
import csv
import glob as globlib
import html
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


CLUSTER_COLUMNS = [
    "Cluster ID",
    "Result ID",
    "Cell Type",
    "True Cell Type",
    "cluster",
    "Cluster",
    "cluster_id",
    "id",
    "group",
    "leiden",
    "seurat_clusters",
]

MAIN_COLUMNS = [
    "Predicted General Cell Type",
    "General Cell Type LLM",
    "General Cell Type Oncology",
    "main_cell_type",
    "final_cell_type",
    "boost_cell_type",
    "original_cell_type",
    "cell_type",
    "celltype1",
    "cell_type_1",
    "primary_cell_type",
    "consensus_main_cell_type",
]

SUB_COLUMNS = [
    "Predicted Detailed Cell Type",
    "Predicted Sub Cell Types",
    "Sub Cell Type LLM",
    "Sub Cell Type Oncology",
    "sub_cell_type",
    "sub_cell_types",
    "final_sub_cell_type",
    "boost_sub_cell_type",
    "subtype",
    "celltype2",
    "cell_type_2",
    "secondary_cell_type",
    "consensus_sub_cell_type",
]

CONFIDENCE_COLUMNS = [
    "Confidence",
    "confidence",
    "boost_confidence",
    "original_confidence",
]

EVIDENCE_COLUMNS = [
    "Evidence",
    "evidence",
    "reason",
    "selection_reasons",
    "error",
]

GENERIC_OUTPUT_NAMES = {
    "summary.csv",
    "subcluster_results.csv",
    "auto_summary.csv",
    "consensus.csv",
}

CONSENSUS_COLUMNS = [
    "cluster_id",
    "consensus_main_cell_type",
    "consensus_sub_cell_type",
    "status",
    "n_votes",
    "n_inputs",
    "main_agreement",
    "sub_agreement",
    "combined_agreement",
    "main_vote_counts",
    "sub_vote_counts",
    "combined_vote_counts",
    "sources",
    "source_annotations",
]


@dataclass(frozen=True)
class AnnotationVote:
    """One annotation vote loaded from a CASSIA output table."""

    source: str
    path: str
    cluster_id: str
    main_label: str
    sub_label: str
    confidence: str = ""
    evidence: str = ""


def _field_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def _label_key(value: str) -> str:
    return re.sub(r"[\s_-]+", " ", value.casefold()).strip()


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text:
        return ""
    return re.sub(r"\s+", " ", text)


def _as_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [_clean_text(item) for item in value if _clean_text(item)]

    text = _clean_text(value)
    if not text:
        return []

    if text[0] in "[(" and text[-1] in "])":
        for parser in (json.loads, ast.literal_eval):
            try:
                parsed = parser(text)
            except Exception:
                continue
            if isinstance(parsed, (list, tuple, set)):
                return [_clean_text(item) for item in parsed if _clean_text(item)]

    normalized = text.replace(";", ",").replace("|", ",")
    return [item.strip() for item in normalized.split(",") if item.strip()]


def _primary_label(value: Any) -> str:
    labels = _as_list(value)
    return labels[0] if labels else _clean_text(value)


def _choose_column(
    headers: Sequence[str],
    candidates: Sequence[str],
    override: Optional[str] = None,
    required: bool = False,
    role: str = "column",
) -> Optional[str]:
    header_map = {_field_key(header): header for header in headers}
    if override:
        key = _field_key(override)
        if key not in header_map:
            raise ValueError(f"Could not find {role} '{override}'. Available columns: {', '.join(headers)}")
        return header_map[key]

    for candidate in candidates:
        key = _field_key(candidate)
        if key in header_map:
            return header_map[key]

    if required:
        raise ValueError(
            f"Could not auto-detect {role}. Available columns: {', '.join(headers)}"
        )
    return None


def _source_name(path: Path) -> str:
    if path.name in GENERIC_OUTPUT_NAMES and path.parent.name:
        return path.parent.name
    return path.stem


def _resolve_directory_input(path: Path) -> Optional[Path]:
    candidates = [
        path / "summary.csv",
        path / "subcluster_results.csv",
        path / "auto_summary.csv",
        path / "boost" / "_auto" / "auto_summary.csv",
    ]
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate
    return None


def expand_consensus_inputs(inputs: Optional[Iterable[str]], patterns: Optional[Iterable[str]]) -> List[Path]:
    """Expand explicit file/directory paths and glob patterns into input CSV paths."""
    paths: List[Path] = []
    for raw_path in inputs or []:
        path = Path(raw_path).expanduser()
        if path.is_dir():
            resolved = _resolve_directory_input(path)
            if resolved is None:
                raise ValueError(f"No supported CASSIA output CSV found in directory: {path}")
            paths.append(resolved)
        else:
            paths.append(path)

    for pattern in patterns or []:
        for match in sorted(globlib.glob(str(Path(pattern).expanduser()), recursive=True)):
            path = Path(match)
            if path.is_dir():
                resolved = _resolve_directory_input(path)
                if resolved is not None:
                    paths.append(resolved)
            else:
                paths.append(path)

    deduped: List[Path] = []
    seen = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        if not resolved.exists():
            raise ValueError(f"Consensus input does not exist: {path}")
        if not resolved.is_file():
            raise ValueError(f"Consensus input is not a file: {path}")
        seen.add(resolved)
        deduped.append(resolved)

    if not deduped:
        raise ValueError("Provide at least one input CSV with --inputs or --glob")
    return deduped


def load_annotation_votes(
    path: Path,
    cluster_column: Optional[str] = None,
    main_column: Optional[str] = None,
    sub_column: Optional[str] = None,
) -> List[AnnotationVote]:
    """Load one CASSIA summary-like CSV into normalized annotation votes."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        if not headers:
            raise ValueError(f"Input CSV has no header: {path}")

        cluster_col = _choose_column(headers, CLUSTER_COLUMNS, cluster_column, required=True, role="cluster column")
        main_col = _choose_column(headers, MAIN_COLUMNS, main_column, required=False, role="main annotation column")
        sub_col = _choose_column(headers, SUB_COLUMNS, sub_column, required=False, role="sub annotation column")
        confidence_col = _choose_column(headers, CONFIDENCE_COLUMNS, required=False, role="confidence column")
        evidence_col = _choose_column(headers, EVIDENCE_COLUMNS, required=False, role="evidence column")

        if main_col is None and sub_col is None:
            raise ValueError(
                f"Could not auto-detect annotation columns in {path}. "
                "Use --main-column and/or --sub-column."
            )

        source = _source_name(path)
        votes: List[AnnotationVote] = []
        for row in reader:
            cluster_id = _clean_text(row.get(cluster_col, ""))
            if not cluster_id:
                continue
            main_label = _primary_label(row.get(main_col, "")) if main_col else ""
            sub_label = _primary_label(row.get(sub_col, "")) if sub_col else ""
            if not main_label and sub_label:
                main_label = sub_label
            if not sub_label and main_label:
                sub_label = main_label
            if not main_label and not sub_label:
                continue
            votes.append(
                AnnotationVote(
                    source=source,
                    path=str(path),
                    cluster_id=cluster_id,
                    main_label=main_label,
                    sub_label=sub_label,
                    confidence=_clean_text(row.get(confidence_col, "")) if confidence_col else "",
                    evidence=_clean_text(row.get(evidence_col, "")) if evidence_col else "",
                )
            )
    return votes


def _natural_sort_key(value: str) -> Tuple[Any, ...]:
    pieces = re.split(r"(\d+)", str(value))
    return tuple(int(piece) if piece.isdigit() else piece.casefold() for piece in pieces)


def _vote_summary(labels: Sequence[str]) -> Dict[str, Any]:
    key_counts: Counter[str] = Counter()
    displays: Dict[str, Counter[str]] = defaultdict(Counter)
    for label in labels:
        clean = _clean_text(label)
        if not clean:
            continue
        key = _label_key(clean)
        key_counts[key] += 1
        displays[key][clean] += 1

    total = sum(key_counts.values())
    if total == 0:
        return {
            "label": "",
            "count": 0,
            "agreement": 0.0,
            "counts": {},
            "tie": False,
        }

    def display_for_key(key: str) -> str:
        return sorted(
            displays[key].items(),
            key=lambda item: (-item[1], item[0].casefold(), item[0]),
        )[0][0]

    def sort_item(item: Tuple[str, int]) -> Tuple[int, str, str]:
        key, count = item
        display = display_for_key(key)
        return (-count, display.casefold(), display)

    sorted_counts = sorted(key_counts.items(), key=sort_item)
    top_key, top_count = sorted_counts[0]
    top_display = display_for_key(top_key)
    top_count_value = top_count
    tie = sum(1 for count in key_counts.values() if count == top_count_value) > 1
    display_counts = {
        display_for_key(key): count
        for key, count in sorted_counts
    }
    return {
        "label": top_display,
        "count": top_count,
        "agreement": top_count / total,
        "counts": display_counts,
        "tie": tie,
    }


def _pair_summary(votes: Sequence[AnnotationVote]) -> Dict[str, Any]:
    pair_counts: Counter[Tuple[str, str]] = Counter()
    pair_displays: Dict[Tuple[str, str], Counter[Tuple[str, str]]] = defaultdict(Counter)
    for vote in votes:
        main = _clean_text(vote.main_label)
        sub = _clean_text(vote.sub_label)
        if not main and not sub:
            continue
        key = (_label_key(main), _label_key(sub))
        pair_counts[key] += 1
        pair_displays[key][(main, sub)] += 1

    total = sum(pair_counts.values())
    if total == 0:
        return {
            "main": "",
            "sub": "",
            "count": 0,
            "agreement": 0.0,
            "counts": {},
            "tie": False,
        }

    def display_for_pair(key: Tuple[str, str]) -> Tuple[str, str, str]:
        main, sub = sorted(
            pair_displays[key].items(),
            key=lambda item: (
                -item[1],
                " / ".join(part for part in item[0] if part).casefold(),
                " / ".join(part for part in item[0] if part),
            ),
        )[0][0]
        label = " / ".join(part for part in (main, sub) if part)
        return main, sub, label

    def sort_pair(item: Tuple[Tuple[str, str], int]) -> Tuple[int, str, str]:
        key, count = item
        _, _, label = display_for_pair(key)
        return (-count, label.casefold(), label)

    sorted_counts = sorted(pair_counts.items(), key=sort_pair)
    top_key, top_count = sorted_counts[0]
    top_main, top_sub, _ = display_for_pair(top_key)
    tie = sum(1 for count in pair_counts.values() if count == top_count) > 1
    display_counts = {}
    for key, count in sorted_counts:
        _, _, label = display_for_pair(key)
        display_counts[label] = count
    return {
        "main": top_main,
        "sub": top_sub,
        "count": top_count,
        "agreement": top_count / total,
        "counts": display_counts,
        "tie": tie,
    }


def _format_fraction(value: float) -> str:
    return f"{value:.3f}"


def _status_for_vote(
    vote_count: int,
    threshold: float,
    main_summary: Dict[str, Any],
    pair_summary: Dict[str, Any],
) -> str:
    if vote_count == 0:
        return "no_votes"
    if vote_count == 1:
        return "single_vote"
    if not pair_summary["tie"] and pair_summary["agreement"] >= threshold:
        return "consensus"
    if not main_summary["tie"] and main_summary["agreement"] >= threshold:
        return "partial_consensus"
    return "conflict"


def build_consensus_rows(votes: Sequence[AnnotationVote], n_inputs: int, threshold: float) -> List[Dict[str, str]]:
    """Build consensus rows from normalized annotation votes."""
    by_cluster: Dict[str, List[AnnotationVote]] = defaultdict(list)
    for vote in votes:
        by_cluster[vote.cluster_id].append(vote)

    rows: List[Dict[str, str]] = []
    for cluster_id in sorted(by_cluster, key=_natural_sort_key):
        cluster_votes = by_cluster[cluster_id]
        main_summary = _vote_summary([vote.main_label for vote in cluster_votes])
        sub_summary = _vote_summary([vote.sub_label for vote in cluster_votes])
        pair = _pair_summary(cluster_votes)
        status = _status_for_vote(len(cluster_votes), threshold, main_summary, pair)
        source_annotations = [
            {
                "source": vote.source,
                "path": vote.path,
                "main_cell_type": vote.main_label,
                "sub_cell_type": vote.sub_label,
                "confidence": vote.confidence,
                "evidence": vote.evidence,
            }
            for vote in cluster_votes
        ]
        rows.append(
            {
                "cluster_id": cluster_id,
                "consensus_main_cell_type": main_summary["label"] or pair["main"],
                "consensus_sub_cell_type": sub_summary["label"] or pair["sub"],
                "status": status,
                "n_votes": str(len(cluster_votes)),
                "n_inputs": str(n_inputs),
                "main_agreement": _format_fraction(main_summary["agreement"]),
                "sub_agreement": _format_fraction(sub_summary["agreement"]),
                "combined_agreement": _format_fraction(pair["agreement"]),
                "main_vote_counts": json.dumps(main_summary["counts"], ensure_ascii=False, sort_keys=True),
                "sub_vote_counts": json.dumps(sub_summary["counts"], ensure_ascii=False, sort_keys=True),
                "combined_vote_counts": json.dumps(pair["counts"], ensure_ascii=False, sort_keys=True),
                "sources": "; ".join(vote.source for vote in cluster_votes),
                "source_annotations": json.dumps(source_annotations, ensure_ascii=False),
            }
        )
    return rows


def write_consensus_csv(rows: Sequence[Dict[str, str]], out_path: Path) -> Path:
    """Write consensus rows to a CSV file."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONSENSUS_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return out_path


def _status_badge(status: str) -> str:
    palette = {
        "consensus": ("#065f46", "#d1fae5"),
        "partial_consensus": ("#92400e", "#fef3c7"),
        "conflict": ("#991b1b", "#fee2e2"),
        "single_vote": ("#374151", "#f3f4f6"),
        "no_votes": ("#374151", "#f3f4f6"),
    }
    fg, bg = palette.get(status, ("#374151", "#f3f4f6"))
    return (
        f'<span style="display:inline-block;padding:2px 7px;border-radius:999px;'
        f'color:{fg};background:{bg};font-size:12px;font-weight:600;">'
        f"{html.escape(status)}</span>"
    )


def write_consensus_html(
    rows: Sequence[Dict[str, str]],
    input_paths: Sequence[Path],
    html_path: Path,
    threshold: float,
) -> Path:
    """Write a compact HTML report for consensus rows."""
    html_path.parent.mkdir(parents=True, exist_ok=True)
    status_counts = Counter(row["status"] for row in rows)
    input_items = "\n".join(
        f"<li><code>{html.escape(str(path))}</code></li>"
        for path in input_paths
    )
    table_rows = []
    for row in rows:
        source_annotations = json.loads(row["source_annotations"] or "[]")
        source_lines = []
        for item in source_annotations:
            label = item.get("main_cell_type", "")
            subtype = item.get("sub_cell_type", "")
            if subtype and subtype != label:
                label = f"{label} / {subtype}" if label else subtype
            confidence = item.get("confidence", "")
            suffix = f" ({confidence})" if confidence else ""
            source_lines.append(f"{item.get('source', '')}: {label}{suffix}")
        table_rows.append(
            "<tr>"
            f"<td>{html.escape(row['cluster_id'])}</td>"
            f"<td>{html.escape(row['consensus_main_cell_type'])}</td>"
            f"<td>{html.escape(row['consensus_sub_cell_type'])}</td>"
            f"<td>{_status_badge(row['status'])}</td>"
            f"<td>{html.escape(row['combined_agreement'])}</td>"
            f"<td>{html.escape(row['main_vote_counts'])}</td>"
            f"<td>{html.escape(row['sub_vote_counts'])}</td>"
            f"<td>{'<br>'.join(html.escape(line) for line in source_lines)}</td>"
            "</tr>"
        )

    html_doc = """<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CASSIA Consensus Report</title>
  <style>
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 32px; color: #1f2937; }
    h1 { margin-bottom: 6px; }
    .meta { color: #4b5563; margin-top: 0; }
    .summary { display: flex; flex-wrap: wrap; gap: 12px; margin: 18px 0 22px; }
    .metric { border: 1px solid #d1d5db; border-radius: 6px; padding: 10px 12px; min-width: 135px; }
    .metric strong { display: block; font-size: 20px; }
    table { border-collapse: collapse; width: 100%; font-size: 14px; }
    th, td { border: 1px solid #d1d5db; padding: 8px; vertical-align: top; }
    th { background: #f3f4f6; text-align: left; }
    code { background: #f3f4f6; padding: 2px 4px; border-radius: 4px; }
    ul { margin-top: 6px; }
  </style>
</head>
<body>
  <h1>CASSIA Consensus Report</h1>
  <p class="meta">Deterministic vote threshold: <code>__THRESHOLD__</code>. Inputs: __INPUT_COUNT__.</p>
  <div class="summary">
    <div class="metric"><strong>__CLUSTERS__</strong>Clusters</div>
    <div class="metric"><strong>__CONSENSUS__</strong>Consensus</div>
    <div class="metric"><strong>__PARTIAL__</strong>Partial</div>
    <div class="metric"><strong>__CONFLICT__</strong>Conflict</div>
  </div>
  <h2>Inputs</h2>
  <ul>
    __INPUTS__
  </ul>
  <h2>Cluster Consensus</h2>
  <table>
    <thead>
      <tr>
        <th>Cluster</th><th>Consensus Main</th><th>Consensus Subtype</th><th>Status</th>
        <th>Combined Agreement</th><th>Main Votes</th><th>Subtype Votes</th><th>Source Annotations</th>
      </tr>
    </thead>
    <tbody>
      __ROWS__
    </tbody>
  </table>
</body>
</html>
"""
    html_doc = (
        html_doc
        .replace("__THRESHOLD__", f"{threshold:.3f}")
        .replace("__INPUT_COUNT__", str(len(input_paths)))
        .replace("__CLUSTERS__", str(len(rows)))
        .replace("__CONSENSUS__", str(status_counts.get("consensus", 0)))
        .replace("__PARTIAL__", str(status_counts.get("partial_consensus", 0)))
        .replace("__CONFLICT__", str(status_counts.get("conflict", 0)))
        .replace("__INPUTS__", input_items)
        .replace("__ROWS__", "\n".join(table_rows))
    )
    html_path.write_text(html_doc, encoding="utf-8")
    return html_path


def run_consensus(args: Any) -> int:
    """Run deterministic consensus voting from CLI arguments."""
    threshold = float(args.threshold)
    if threshold <= 0 or threshold > 1:
        raise ValueError("--threshold must be greater than 0 and at most 1")

    input_paths = expand_consensus_inputs(args.inputs, args.glob_patterns)
    all_votes: List[AnnotationVote] = []
    for path in input_paths:
        all_votes.extend(
            load_annotation_votes(
                path,
                cluster_column=args.cluster_column,
                main_column=args.main_column,
                sub_column=args.sub_column,
            )
        )
    if not all_votes:
        raise ValueError("No annotation votes were found in the provided inputs")

    rows = build_consensus_rows(all_votes, n_inputs=len(input_paths), threshold=threshold)
    out_path = Path(args.out)
    write_consensus_csv(rows, out_path)
    print(f"Wrote {out_path}")

    if not args.no_html:
        html_path = Path(args.html) if args.html else out_path.with_suffix(".html")
        write_consensus_html(rows, input_paths, html_path, threshold=threshold)
        print(f"Wrote {html_path}")
    return 0
