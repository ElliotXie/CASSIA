"""Stable, blinded LLM-as-judge evaluation for CASSIA benchmark outputs."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import random
import re
import string
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from .backends import AGENT_BACKENDS, AgentCLIBackend
from .runner import extract_json_object, utc_now


JUDGE_PROTOCOL_VERSION = "stable-judge-v1.1"
DEFAULT_JUDGE_BACKENDS = ("cursor-agent",)
DEFAULT_JUDGE_MODELS = {
    "cursor-agent": "composer-2.5",
    "claude-cli": "claude-sonnet-4-6",
    "codex-cli": "gpt-5.6-luna",
}
DEFAULT_BATCH_SIZE = 4
VALID_VERDICTS = {"correct", "partial", "wrong", "missing"}
VALID_MAJOR_ERRORS = {
    "none", "missing_output", "wrong_lineage", "wrong_subtype", "wrong_state",
    "unsupported_markers", "overclaim", "too_generic",
}

CASE_ID_CANDIDATES = ["case_id", "cluster_id", "Result ID", "Cluster ID", "cluster", "id"]
EXPECTED_LABEL_CANDIDATES = [
    "expected_label", "expected_cell_type", "label", "ground_truth", "Expected Label", "cell_type",
]
EXPECTED_TERMS_CANDIDATES = ["expected_terms", "rubric", "terms"]
SOURCE_CANDIDATES = ["source_sheet", "source", "paper_cluster", "cluster_name"]
MARKER_CANDIDATES = ["markers", "marker_list", "marker_genes", "Marker List", "key_markers"]
MAIN_TYPE_CANDIDATES = [
    "main_cell_type", "Predicted General Cell Type", "final_cell_type", "general_cell_type",
]
TOP1_CANDIDATES = [
    "top1_subtype", "primary_subtype", "Predicted Detailed Cell Type",
    "final_sub_cell_type", "sub_cell_type", "subtype", "cell_type", "annotation",
]
TOP2_CANDIDATES = ["top2_subtype", "second_subtype"]
TOP3_CANDIDATES = ["top3_subtype", "third_subtype"]
RANKED_CANDIDATES = [
    "sub_cell_types", "all_sub_cell_types", "ranked_subtypes", "top_three_subtypes", "alternatives",
]


@dataclass(frozen=True)
class PredictionSpec:
    mode: str
    path: Path


def default_judge_dir() -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("cassia_runs") / f"judge_{stamp}"


def default_cache_dir() -> Path:
    configured = os.environ.get("CASSIA_JUDGE_CACHE_DIR")
    return Path(configured).expanduser() if configured else Path.home() / ".cache" / "cassia" / JUDGE_PROTOCOL_VERSION


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _find_column(columns: Sequence[str], candidates: Sequence[str], override: Optional[str] = None) -> str:
    mapping = {str(column).casefold(): str(column) for column in columns}
    if override:
        found = mapping.get(str(override).casefold())
        if found:
            return found
        raise ValueError(f"Column '{override}' not found. Available: {', '.join(map(str, columns))}")
    for candidate in candidates:
        found = mapping.get(candidate.casefold())
        if found:
            return found
    raise ValueError(f"None of the expected columns were found: {', '.join(candidates)}")


def _optional_column(columns: Sequence[str], candidates: Sequence[str], override: Optional[str] = None) -> Optional[str]:
    if override:
        return _find_column(columns, candidates, override)
    mapping = {str(column).casefold(): str(column) for column in columns}
    for candidate in candidates:
        found = mapping.get(candidate.casefold())
        if found:
            return found
    return None


def _cell(row: pd.Series, column: Optional[str]) -> str:
    if not column:
        return ""
    value = row.get(column, "")
    return "" if value is None or pd.isna(value) else str(value).strip()


def _compact_text(value: object, max_chars: int) -> str:
    text = " ".join(str(value or "").split())
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 20].rstrip() + " ...[truncated]"


def _first_value(row: pd.Series, candidates: Sequence[str]) -> str:
    mapping = {str(column).casefold(): str(column) for column in row.index}
    for candidate in candidates:
        column = mapping.get(candidate.casefold())
        if column:
            value = _cell(row, column)
            if value:
                return value
    return ""


def _ranked_values(value: Any) -> List[str]:
    if value is None or (not isinstance(value, (list, tuple)) and pd.isna(value)):
        return []
    if isinstance(value, (list, tuple)):
        raw = list(value)
    else:
        text = str(value).strip()
        raw = []
        if text.startswith("["):
            try:
                parsed = json.loads(text)
                if isinstance(parsed, list):
                    raw = parsed
            except json.JSONDecodeError:
                pass
        if not raw:
            raw = re.split(r"\s*;\s*|\s*\|\s*", text)
    values: List[str] = []
    seen = set()
    for item in raw:
        cleaned = str(item).strip().strip('"\'')
        key = cleaned.casefold()
        if cleaned and key not in seen and key != "nan":
            values.append(cleaned)
            seen.add(key)
    return values


def _canonical_prediction(row: Optional[pd.Series]) -> Dict[str, Any]:
    if row is None:
        return {
            "missing_output": True,
            "main_cell_type": "",
            "top1_subtype": "",
            "top2_subtype": "",
            "top3_subtype": "",
        }
    main = _first_value(row, MAIN_TYPE_CANDIDATES)
    explicit = [
        _first_value(row, TOP1_CANDIDATES),
        _first_value(row, TOP2_CANDIDATES),
        _first_value(row, TOP3_CANDIDATES),
    ]
    ranked: List[str] = []
    mapping = {str(column).casefold(): str(column) for column in row.index}
    for candidate in RANKED_CANDIDATES:
        column = mapping.get(candidate.casefold())
        if column:
            ranked = _ranked_values(row.get(column))
            if ranked:
                break
    combined: List[str] = []
    seen = set()
    for value in [explicit[0], *ranked, explicit[1], explicit[2]]:
        key = value.casefold() if value else ""
        if value and key not in seen:
            combined.append(value)
            seen.add(key)
    top = (combined + ["", "", ""])[:3]
    return {
        "missing_output": not bool(main or top[0]),
        "main_cell_type": _compact_text(main, 240),
        "top1_subtype": _compact_text(top[0], 240),
        "top2_subtype": _compact_text(top[1], 240),
        "top3_subtype": _compact_text(top[2], 240),
    }


def parse_prediction_spec(value: str) -> PredictionSpec:
    if ":" not in value:
        raise ValueError("Prediction must be MODE:CSV_PATH")
    mode, path = value.split(":", 1)
    if not mode.strip() or not path.strip():
        raise ValueError("Prediction must be MODE:CSV_PATH")
    return PredictionSpec(mode.strip(), Path(path.strip()))


def load_truth(args: Any) -> pd.DataFrame:
    path = Path(args.truth)
    frame = pd.read_csv(path)
    case_col = _find_column(frame.columns, CASE_ID_CANDIDATES, args.truth_case_column)
    label_col = _find_column(frame.columns, EXPECTED_LABEL_CANDIDATES, args.expected_label_column)
    terms_col = _optional_column(frame.columns, EXPECTED_TERMS_CANDIDATES, args.expected_terms_column)
    source_col = _optional_column(frame.columns, SOURCE_CANDIDATES, args.source_column)
    marker_col = _optional_column(frame.columns, MARKER_CANDIDATES, args.truth_marker_column)
    rows = []
    for _, row in frame.iterrows():
        case_id = _cell(row, case_col)
        if case_id:
            rows.append({
                "case_id": case_id,
                "expected_label": _cell(row, label_col),
                "expected_terms": _cell(row, terms_col),
                "source": _cell(row, source_col),
                "markers": _cell(row, marker_col),
            })
    if not rows:
        raise ValueError(f"No truth cases loaded from {path}")
    result = pd.DataFrame(rows)
    if result["case_id"].duplicated().any():
        raise ValueError("Truth CSV contains duplicate case IDs")
    return result


def load_markers(args: Any) -> Dict[str, str]:
    if not args.markers:
        return {}
    frame = pd.read_csv(args.markers)
    case_col = _find_column(frame.columns, CASE_ID_CANDIDATES, args.marker_case_column)
    marker_col = _find_column(frame.columns, MARKER_CANDIDATES, args.marker_column)
    return {
        _cell(row, case_col): _cell(row, marker_col)
        for _, row in frame.iterrows()
        if _cell(row, case_col)
    }


def load_prediction(path: Path, case_column: Optional[str] = None) -> pd.DataFrame:
    frame = pd.read_csv(path)
    case_col = _find_column(frame.columns, CASE_ID_CANDIDATES, case_column)
    frame = frame.copy()
    frame["__case_id__"] = frame[case_col].astype(str)
    if frame["__case_id__"].duplicated().any():
        raise ValueError(f"Prediction CSV contains duplicate case IDs: {path}")
    return frame.set_index("__case_id__", drop=False)


def build_judge_items(
    truth: pd.DataFrame,
    markers: Mapping[str, str],
    predictions: Mapping[str, pd.DataFrame],
    case_ids: Sequence[str],
) -> List[Dict[str, Any]]:
    truth_map = {str(row.case_id): row for row in truth.itertuples(index=False)}
    items: List[Dict[str, Any]] = []
    for case_id in case_ids:
        truth_row = truth_map[case_id]
        marker_text = markers.get(case_id) or getattr(truth_row, "markers") or ""
        for mode, prediction_frame in predictions.items():
            pred_row = None if case_id not in prediction_frame.index else prediction_frame.loc[case_id]
            items.append({
                "case_id": case_id,
                "mode": mode,
                "source": getattr(truth_row, "source"),
                "expected_label": getattr(truth_row, "expected_label"),
                "expected_terms": getattr(truth_row, "expected_terms"),
                "marker_genes": _compact_text(marker_text, 1400),
                "prediction": _canonical_prediction(pred_row),
            })
    return items


def _blind_items(
    items: Sequence[Mapping[str, Any]], seed: int,
) -> Tuple[List[Dict[str, Any]], Dict[Tuple[str, str], str]]:
    rng = random.Random(seed)
    by_case: Dict[str, List[Dict[str, Any]]] = {}
    for item in items:
        by_case.setdefault(str(item["case_id"]), []).append(dict(item))
    decoder: Dict[Tuple[str, str], str] = {}
    blinded: List[Dict[str, Any]] = []
    for case_id, case_items in by_case.items():
        labels = [f"pred_{letter}" for letter in string.ascii_uppercase[: len(case_items)]]
        rng.shuffle(labels)
        for item, label in zip(case_items, labels):
            decoder[(case_id, label)] = str(item["mode"])
            item["mode"] = label
            blinded.append(item)
    rng.shuffle(blinded)
    return blinded, decoder


def _group_prompt_items(items: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[str, Dict[str, Any]] = {}
    for item in items:
        case_id = str(item["case_id"])
        case = grouped.setdefault(case_id, {
            "case_id": case_id,
            "source": item.get("source", ""),
            "expected_label": item.get("expected_label", ""),
            "expected_terms": item.get("expected_terms", ""),
            "marker_genes": item.get("marker_genes", ""),
            "predictions": [],
        })
        case["predictions"].append({
            "mode": item["mode"],
            **dict(item["prediction"]),
        })
    return list(grouped.values())


def build_judge_prompt(items: Sequence[Mapping[str, Any]], dataset_context: str = "") -> str:
    cases = json.dumps(_group_prompt_items(items), indent=2, ensure_ascii=False)
    expected_count = len(items)
    required_pairs = json.dumps([
        {"case_id": str(item["case_id"]), "mode": str(item["mode"])}
        for item in items
    ], indent=2, ensure_ascii=False)
    return f"""You are an independent expert judge for single-cell RNA-seq annotation.

Protocol: {JUDGE_PROTOCOL_VERSION}

Evaluate each anonymized prediction independently against the held-out label and the
fixed marker evidence. Never infer an annotator identity from pred_A/pred_B labels.

Primary-call rule (mandatory):
- main_cell_type is the primary input for lineage_score.
- top1_subtype is the ONLY ranked subtype allowed to affect top1_verdict,
  lineage_score, subtype_score, state_score, marker_support_score, or pass status.
- top2_subtype and top3_subtype are diagnostic alternatives only. They may determine
  matched_rank, but a correct label at rank 2 or 3 gives ZERO credit to all core scores.
- Judge whether top1 itself is exact, broad/near, or wrong. Do not let a lower-ranked
  exact answer turn a wrong top1 into partial.

Four axes, each 0-2 when applicable:
1. lineage_score: 2 correct broad lineage; 1 related but unresolved broad lineage;
   0 wrong/contradictory lineage or missing.
2. subtype_score: 2 top1 recovers the expected subtype or a fully compatible finer
   subtype; 1 top1 itself is a defensible broader/near subtype; 0 wrong/generic/missing.
3. state_score: 2 expected functional state/stage is recovered by top1; 1 partially
   recovered; 0 wrong/absent when the truth requires a state distinction; null when
   the ground truth contains no meaningful state/stage distinction. Apply this
   mechanically: an explicitly named related but non-equivalent state/stage is 1;
   omission or contradiction of a required state/stage is 0.
4. marker_support_score: 2 main+top1 are strongly supported by the supplied markers;
   1 broadly defensible but incomplete/ambiguous; 0 contradicted or unsupported.
   Score 2 only when markers support the distinguishing subtype/state in top1, not
   merely its broad lineage. If markers support lineage but do not resolve top1 from
   plausible siblings, score 1. Score 0 for contradiction or absent support;
   null only when marker evidence is absent. Judge marker compatibility, not writing
   length: no rationale is required from the annotator.

top1_verdict:
- correct: top1 is equivalent to the truth, a standard synonym, or fully compatible finer subtype.
- partial: top1 itself is correct-lineage but too broad, or a defensible near-neighbor.
- wrong: top1 is an incorrect sibling/state/lineage.
- missing: no usable main/top1 output.

matched_rank is the semantic rank of the ground-truth-equivalent label among top1/top2/top3,
or null if absent. It is diagnostic and never changes the core score.

Dataset context:
{dataset_context.strip() or 'Not provided.'}

Cases:
{cases}

Required case_id/mode pairs (return each pair exactly once and copy both strings verbatim):
{required_pairs}

Return only JSON with exactly {expected_count} judgments:
{{
  "protocol_version": "{JUDGE_PROTOCOL_VERSION}",
  "judgments": [
    {{
      "case_id": "exact case_id",
      "mode": "pred_A",
      "top1_verdict": "correct|partial|wrong|missing",
      "matched_rank": 1,
      "lineage_score": 0,
      "subtype_score": 0,
      "state_score": null,
      "marker_support_score": 0,
      "major_error": "none|missing_output|wrong_lineage|wrong_subtype|wrong_state|unsupported_markers|overclaim|too_generic",
      "rationale": "short reason based on main and top1; mention lower-rank recovery only diagnostically"
    }}
  ]
}}

Preserve every case_id and blinded mode exactly. No markdown or commentary.
"""


def _extract_judgments(payload: Mapping[str, Any]) -> List[Mapping[str, Any]]:
    judgments = payload.get("judgments")
    if not isinstance(judgments, list):
        raise ValueError("Judge JSON must contain a judgments list")
    return [item for item in judgments if isinstance(item, Mapping)]


def validate_judge_payload(payload: Mapping[str, Any], expected_items: Sequence[Mapping[str, Any]]) -> None:
    judgments = _extract_judgments(payload)
    expected = {(str(item["case_id"]), str(item["mode"])) for item in expected_items}
    returned = [(str(item.get("case_id", "")), str(item.get("mode", ""))) for item in judgments]
    if len(judgments) != len(expected_items) or set(returned) != expected or len(set(returned)) != len(returned):
        raise ValueError("Judge must return exactly one judgment for every case_id/mode")
    for item in judgments:
        verdict = str(item.get("top1_verdict", ""))
        if verdict not in VALID_VERDICTS:
            raise ValueError(f"Invalid top1_verdict: {verdict}")
        for field in ("lineage_score", "subtype_score"):
            if item.get(field) not in (0, 1, 2):
                raise ValueError(f"{field} must be 0, 1, or 2")
        for field in ("state_score", "marker_support_score"):
            if item.get(field) not in (None, 0, 1, 2):
                raise ValueError(f"{field} must be null, 0, 1, or 2")
        if item.get("matched_rank") not in (None, 1, 2, 3):
            raise ValueError("matched_rank must be null, 1, 2, or 3")


def _int_score(value: Any) -> int:
    return max(0, min(2, int(value)))


def normalize_judgments(payload: Mapping[str, Any], expected_items: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    validate_judge_payload(payload, expected_items)
    raw_map = {
        (str(item["case_id"]), str(item["mode"])): item
        for item in _extract_judgments(payload)
    }
    rows = []
    for expected in expected_items:
        key = (str(expected["case_id"]), str(expected["mode"]))
        item = raw_map[key]
        missing = bool(expected["prediction"].get("missing_output"))
        verdict = "missing" if missing else str(item["top1_verdict"])
        lineage = 0 if missing else _int_score(item["lineage_score"])
        subtype = 0 if missing else _int_score(item["subtype_score"])
        state = None if item.get("state_score") is None else _int_score(item["state_score"])
        marker = None if item.get("marker_support_score") is None else _int_score(item["marker_support_score"])
        if missing:
            state = 0 if state is not None else None
            marker = 0 if marker is not None else None
        applicable = [lineage, subtype] + [score for score in (state, marker) if score is not None]
        total = int(sum(applicable))
        maximum = len(applicable) * 2
        normalized = total / maximum if maximum else 0.0
        major_error = "missing_output" if missing else str(item.get("major_error") or "none")
        if major_error not in VALID_MAJOR_ERRORS:
            major_error = "none"
        matched_rank = item.get("matched_rank")
        quality_pass = bool(
            verdict in {"correct", "partial"}
            and normalized >= 0.75
            and lineage >= 1
            and subtype >= 1
            and marker != 0
            and major_error not in {"missing_output", "wrong_lineage"}
        )
        strict_pass = bool(
            verdict == "correct"
            and lineage == 2
            and subtype == 2
            and (state is None or state >= 1)
            and (marker is None or marker >= 1)
            and major_error == "none"
        )
        rows.append({
            "case_id": key[0],
            "mode": key[1],
            "source": expected.get("source", ""),
            "expected_label": expected.get("expected_label", ""),
            "expected_terms": expected.get("expected_terms", ""),
            "marker_genes": expected.get("marker_genes", ""),
            **{f"prediction_{name}": value for name, value in expected["prediction"].items()},
            "top1_verdict": verdict,
            "matched_rank": matched_rank,
            "lower_rank_only_recovery": matched_rank in {2, 3},
            "lineage_score": lineage,
            "subtype_score": subtype,
            "state_score": state,
            "marker_support_score": marker,
            "judge_total": total,
            "judge_max": maximum,
            "normalized_score": round(normalized, 6),
            "judge_pass": quality_pass,
            "strict_pass": strict_pass,
            "major_error": major_error,
            "rationale": str(item.get("rationale") or ""),
        })
    return pd.DataFrame(rows)


def _mean_nullable(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce")
    return float(values.mean()) if values.notna().any() else float("nan")


def summarize_scores(scored: pd.DataFrame, by: Sequence[str] = ("mode",)) -> pd.DataFrame:
    rows = []
    for keys, group in scored.groupby(list(by), sort=False, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        verdicts = group["top1_verdict"].value_counts()
        rows.append({
            **{column: key for column, key in zip(by, keys)},
            "judged_rows": int(len(group)),
            "unique_cases": int(group["case_id"].nunique()),
            "mean_normalized_score": float(group["normalized_score"].mean()),
            "correct_count": int(verdicts.get("correct", 0)),
            "partial_count": int(verdicts.get("partial", 0)),
            "wrong_count": int(verdicts.get("wrong", 0)),
            "missing_count": int(verdicts.get("missing", 0)),
            "lower_rank_only_count": int(group["lower_rank_only_recovery"].sum()),
            "judge_pass_count": int(group["judge_pass"].sum()),
            "strict_pass_count": int(group["strict_pass"].sum()),
            "mean_lineage_score": float(group["lineage_score"].mean()),
            "mean_subtype_score": float(group["subtype_score"].mean()),
            "mean_state_score": _mean_nullable(group["state_score"]),
            "mean_marker_support_score": _mean_nullable(group["marker_support_score"]),
        })
    return pd.DataFrame(rows)


def build_consensus_scores(scored: pd.DataFrame) -> pd.DataFrame:
    """Keep per-case inter-judge agreement explicit instead of hiding it in a mean."""
    rows = []
    for (case_id, mode), group in scored.groupby(["case_id", "mode"], sort=False):
        verdicts = [str(value) for value in group["top1_verdict"]]
        ranks = [None if pd.isna(value) else int(value) for value in group["matched_rank"]]
        unique_verdicts = sorted(set(verdicts))
        unique_ranks = sorted(set(ranks), key=lambda value: (-1 if value is None else value))
        row = {
            "case_id": case_id,
            "mode": mode,
            "expected_label": group.iloc[0]["expected_label"],
            "judge_count": int(len(group)),
            "judge_backends": ";".join(group["judge_backend"].astype(str)),
            "top1_verdict_consensus": unique_verdicts[0] if len(unique_verdicts) == 1 else "disagreement",
            "top1_verdicts": ";".join(verdicts),
            "verdict_agreement": len(unique_verdicts) == 1,
            "matched_rank_consensus": unique_ranks[0] if len(unique_ranks) == 1 else None,
            "matched_ranks": ";".join("NA" if value is None else str(value) for value in ranks),
            "matched_rank_agreement": len(unique_ranks) == 1,
            "mean_normalized_score": float(group["normalized_score"].mean()),
            "normalized_score_range": float(group["normalized_score"].max() - group["normalized_score"].min()),
            "all_judges_pass": bool(group["judge_pass"].all()),
            "all_judges_strict": bool(group["strict_pass"].all()),
        }
        for axis in ("lineage_score", "subtype_score", "state_score", "marker_support_score"):
            values = pd.to_numeric(group[axis], errors="coerce")
            row[f"mean_{axis}"] = float(values.mean()) if values.notna().any() else None
            row[f"{axis}_agreement"] = bool(values.dropna().nunique() <= 1 and values.notna().nunique() <= 1)
        rows.append(row)
    return pd.DataFrame(rows)


def _safe_path_part(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value or "default")


def score_cache_key(
    item: Mapping[str, Any],
    context: str,
    backend: str,
    model: str,
    reasoning_effort: str = "",
) -> str:
    payload = {
        "protocol": JUDGE_PROTOCOL_VERSION,
        "backend": backend,
        "model": model,
        "dataset_context": context.strip(),
        "case": {
            "case_id": item.get("case_id"),
            "source": item.get("source"),
            "expected_label": item.get("expected_label"),
            "expected_terms": item.get("expected_terms"),
            "marker_genes": item.get("marker_genes"),
            "prediction": item.get("prediction"),
        },
    }
    # Preserve stable-judge-v1.1 cache identity for backends that do not expose
    # an effort setting. Only an explicitly pinned effort creates a new arm.
    if reasoning_effort:
        payload["reasoning_effort"] = reasoning_effort
    canonical = json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _empty_effort_compat_cache_key(
    item: Mapping[str, Any], context: str, backend: str, model: str,
) -> str:
    """Read keys produced briefly when an empty effort was hashed explicitly."""
    payload = {
        "protocol": JUDGE_PROTOCOL_VERSION,
        "backend": backend,
        "model": model,
        "reasoning_effort": "",
        "dataset_context": context.strip(),
        "case": {
            "case_id": item.get("case_id"),
            "source": item.get("source"),
            "expected_label": item.get("expected_label"),
            "expected_terms": item.get("expected_terms"),
            "marker_genes": item.get("marker_genes"),
            "prediction": item.get("prediction"),
        },
    }
    canonical = json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _cache_path(
    cache_dir: Path,
    backend: str,
    model: str,
    reasoning_effort: str,
    key: str,
) -> Path:
    identity = f"{model}@{reasoning_effort}" if reasoning_effort else model
    return cache_dir / _safe_path_part(backend) / _safe_path_part(identity) / f"{key}.json"


def _load_cached_row(path: Path, item: Mapping[str, Any]) -> Optional[pd.DataFrame]:
    try:
        cached = json.loads(path.read_text(encoding="utf-8"))
        row = dict(cached["row"])
    except (OSError, KeyError, json.JSONDecodeError, TypeError):
        return None
    row["case_id"] = item["case_id"]
    row["mode"] = item["mode"]
    return pd.DataFrame([row])


def _write_cache_row(
    path: Path,
    row: Mapping[str, Any],
    backend: str,
    model: str,
    reasoning_effort: str,
    key: str,
) -> None:
    serializable = {}
    for name, value in row.items():
        if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
            serializable[name] = None
        elif hasattr(value, "item"):
            serializable[name] = value.item()
        else:
            serializable[name] = value
    _write_json(path, {
        "protocol_version": JUDGE_PROTOCOL_VERSION,
        "backend": backend,
        "model": model,
        "reasoning_effort": reasoning_effort,
        "cache_key": key,
        "created_at": utc_now(),
        "row": serializable,
    })


def _decode_scores(scored: pd.DataFrame, decoder: Mapping[Tuple[str, str], str]) -> pd.DataFrame:
    scored = scored.copy()
    scored["mode"] = [
        decoder[(str(case_id), str(mode))]
        for case_id, mode in zip(scored["case_id"], scored["mode"])
    ]
    return scored


def _judge_batch(
    backend: AgentCLIBackend,
    backend_name: str,
    model: str,
    reasoning_effort: str,
    items: Sequence[Mapping[str, Any]],
    batch_index: int,
    out_dir: Path,
    context: str,
    seed: int,
    max_attempts: int,
) -> pd.DataFrame:
    blinded, decoder = _blind_items(items, seed=seed * 1000 + batch_index)
    prompt = build_judge_prompt(blinded, dataset_context=context)
    prompt_path = out_dir / "prompts" / f"judge_batch_{batch_index:04d}__{_safe_path_part(backend_name)}.md"
    decoder_path = out_dir / "prompts" / f"judge_batch_{batch_index:04d}.decoder.json"
    prompt_path.parent.mkdir(parents=True, exist_ok=True)
    prompt_path.write_text(prompt, encoding="utf-8")
    _write_json(decoder_path, [
        {"case_id": case_id, "blind_label": blind, "mode": mode}
        for (case_id, blind), mode in decoder.items()
    ])
    last_error: Optional[Exception] = None
    for attempt in range(1, max_attempts + 1):
        raw_path = out_dir / "raw" / (
            f"judge_batch_{batch_index:04d}__{_safe_path_part(backend_name)}__attempt{attempt}.txt"
        )
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            workspace = (
                out_dir
                / "workspace"
                / _safe_path_part(backend_name)
                / f"batch_{batch_index:04d}"
            )
            workspace.mkdir(parents=True, exist_ok=True)
            output = backend.run(
                prompt=prompt,
                prompt_file=prompt_path,
                cwd=workspace,
                context={
                    "input": "",
                    "out": str(out_dir),
                    "cluster": f"judge_batch_{batch_index:04d}",
                    "agent_output_file": str(raw_path.resolve()),
                },
            )
            raw_path.write_text(output + ("\n" if not output.endswith("\n") else ""), encoding="utf-8")
            payload = extract_json_object(output)
            scored = normalize_judgments(payload, blinded)
            scored = _decode_scores(scored, decoder)
            scored["batch"] = batch_index
            scored["judge_backend"] = backend_name
            scored["judge_model"] = model
            scored["judge_reasoning_effort"] = reasoning_effort
            scored["protocol_version"] = JUDGE_PROTOCOL_VERSION
            return scored
        except Exception as exc:
            last_error = exc
            prompt = build_judge_prompt(blinded, dataset_context=context) + (
                f"\nPrevious attempt was invalid: {_compact_text(exc, 300)}\n"
                "Return a complete corrected JSON object only."
            )
            prompt_path.write_text(prompt, encoding="utf-8")
    raise RuntimeError(f"Judge batch {batch_index} failed after {max_attempts} attempts: {last_error}")


def _chunk_pending_cases(items: Sequence[Mapping[str, Any]], size: int) -> Iterable[List[Mapping[str, Any]]]:
    by_case: Dict[str, List[Mapping[str, Any]]] = {}
    for item in items:
        by_case.setdefault(str(item["case_id"]), []).append(item)
    case_ids = list(by_case)
    for start in range(0, len(case_ids), size):
        batch_ids = case_ids[start:start + size]
        yield [item for case_id in batch_ids for item in by_case[case_id]]


def write_scores_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def write_judge_html(scored: pd.DataFrame, summary: pd.DataFrame, out_dir: Path) -> Path:
    cards = []
    for _, row in summary.iterrows():
        cards.append(
            f"<article><h3>{html.escape(str(row['mode']))}</h3>"
            f"<strong>{row['mean_normalized_score']:.1%}</strong>"
            f"<p>correct {int(row['correct_count'])} · partial {int(row['partial_count'])} · "
            f"strict {int(row['strict_pass_count'])}/{int(row['judged_rows'])}</p></article>"
        )
    rows = []
    for row in scored.itertuples(index=False):
        state = "NA" if pd.isna(row.state_score) else int(row.state_score)
        marker = "NA" if pd.isna(row.marker_support_score) else int(row.marker_support_score)
        rows.append(
            f"<tr><td>{html.escape(str(row.case_id))}</td><td>{html.escape(str(row.mode))}</td>"
            f"<td>{html.escape(str(row.expected_label))}</td>"
            f"<td><b>{html.escape(str(row.top1_verdict))}</b> · rank {row.matched_rank}</td>"
            f"<td>L{int(row.lineage_score)} S{int(row.subtype_score)} T{state} M{marker}</td>"
            f"<td>{row.normalized_score:.1%}</td><td>{html.escape(str(row.rationale))}</td></tr>"
        )
    document = f"""<!doctype html><html><head><meta charset="utf-8"><title>CASSIA stable judge</title>
<style>body{{font:15px system-ui;margin:32px;color:#17201b}}.cards{{display:flex;gap:16px;flex-wrap:wrap}}article{{padding:16px;border:1px solid #ccd5ce;border-radius:12px}}strong{{font-size:28px}}table{{border-collapse:collapse;width:100%;margin-top:24px}}th,td{{padding:8px;border-bottom:1px solid #ddd;text-align:left;vertical-align:top}}</style></head>
<body><h1>CASSIA {JUDGE_PROTOCOL_VERSION}</h1><p>Main type determines lineage; top1 is the only scored subtype. Top2/3 are diagnostic only.</p>
<div class="cards">{''.join(cards)}</div><table><thead><tr><th>case</th><th>mode</th><th>truth</th><th>top1</th><th>axes</th><th>normalized</th><th>rationale</th></tr></thead><tbody>{''.join(rows)}</tbody></table></body></html>"""
    path = out_dir / "judge_report.html"
    path.write_text(document, encoding="utf-8")
    return path


def _model_map(args: Any, backends: Sequence[str]) -> Dict[str, str]:
    overrides: Dict[str, str] = {}
    for value in args.judge_model or []:
        if "=" not in value:
            if len(backends) != 1:
                raise ValueError("Use --judge-model BACKEND=MODEL with multiple backends")
            overrides[backends[0]] = value
        else:
            backend, model = value.split("=", 1)
            overrides[backend.strip()] = model.strip()
    result = {}
    for backend in backends:
        model = overrides.get(backend) or DEFAULT_JUDGE_MODELS.get(backend)
        if not model:
            raise ValueError(f"Pin a model for {backend} with --judge-model {backend}=MODEL")
        result[backend] = model
    return result


def _reasoning_effort_map(args: Any, backends: Sequence[str]) -> Dict[str, str]:
    """Resolve optional per-backend reasoning-effort overrides."""
    overrides: Dict[str, str] = {}
    for value in args.judge_reasoning_effort or []:
        if "=" not in value:
            if len(backends) != 1:
                raise ValueError(
                    "Use --judge-reasoning-effort BACKEND=EFFORT with multiple backends"
                )
            overrides[backends[0]] = value.strip()
        else:
            backend, effort = value.split("=", 1)
            overrides[backend.strip()] = effort.strip()
    unknown = sorted(set(overrides) - set(backends))
    if unknown:
        raise ValueError(f"Reasoning effort supplied for unselected backend: {', '.join(unknown)}")
    for backend, effort in overrides.items():
        if backend != "codex-cli":
            raise ValueError("--judge-reasoning-effort is currently supported only for codex-cli")
        if effort not in AgentCLIBackend._CODEX_REASONING_EFFORTS:
            raise ValueError(f"Unsupported Codex reasoning effort: {effort}")
    return {backend: overrides.get(backend, "") for backend in backends}


def _resolve_backends(requested: Optional[Sequence[str]]) -> List[str]:
    names = list(requested or DEFAULT_JUDGE_BACKENDS)
    for name in names:
        backend = AGENT_BACKENDS.get(name)
        if backend is None:
            raise ValueError(f"Unsupported judge backend: {name}")
        if not backend.available:
            raise RuntimeError(f"Judge backend executable is unavailable: {name}")
    return names


def run_judge(args: Any) -> int:
    if args.batch_size <= 0:
        raise ValueError("--batch-size must be greater than 0")
    if args.max_workers <= 0:
        raise ValueError("--max-workers must be greater than 0")
    out_dir = Path(args.out) if args.out else default_judge_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    specs = [parse_prediction_spec(value) for value in args.prediction]
    truth = load_truth(args)
    markers = load_markers(args)
    predictions = {spec.mode: load_prediction(spec.path, args.prediction_case_column) for spec in specs}
    case_ids = truth["case_id"].astype(str).tolist()
    if args.limit_cases is not None:
        case_ids = case_ids[: args.limit_cases]
    items = build_judge_items(truth, markers, predictions, case_ids)
    backends = _resolve_backends(args.backend)
    models = _model_map(args, backends)
    reasoning_efforts = _reasoning_effort_map(args, backends)
    context = args.context or ""
    cache_dir = Path(args.cache_dir).expanduser() if args.cache_dir else default_cache_dir()
    seed = args.blind_seed

    manifest = {
        "protocol_version": JUDGE_PROTOCOL_VERSION,
        "created_at": utc_now(),
        "status": "running",
        "truth": str(Path(args.truth).resolve()),
        "markers": str(Path(args.markers).resolve()) if args.markers else "",
        "predictions": {spec.mode: str(spec.path.resolve()) for spec in specs},
        "judge_backends": backends,
        "judge_models": models,
        "judge_reasoning_efforts": reasoning_efforts,
        "batch_size_cases": args.batch_size,
        "max_workers": args.max_workers,
        "case_count": len(case_ids),
        "prediction_count": len(items),
        "blind_seed": seed,
        "cache_dir": str(cache_dir.resolve()),
        "top2_top3_affect_core_score": False,
        "parameters": {key: value for key, value in vars(args).items() if key != "func"},
    }
    _write_json(out_dir / "judge_manifest.json", manifest)

    if args.dry_run:
        blinded, decoder = _blind_items(items, seed)
        prompt = build_judge_prompt(blinded, context)
        (out_dir / "prompts").mkdir(parents=True, exist_ok=True)
        (out_dir / "prompts" / "dry_run.md").write_text(prompt, encoding="utf-8")
        _write_json(out_dir / "prompts" / "dry_run.decoder.json", [
            {"case_id": case_id, "blind_label": blind, "mode": mode}
            for (case_id, blind), mode in decoder.items()
        ])
        manifest.update({"status": "dry-run", "updated_at": utc_now()})
        _write_json(out_dir / "judge_manifest.json", manifest)
        return 0

    all_scores: List[pd.DataFrame] = []
    cache_hits = 0
    new_scores = 0
    new_unique_scores = 0
    for backend_name in backends:
        model = models[backend_name]
        reasoning_effort = reasoning_efforts[backend_name]
        backend = AgentCLIBackend(
            backend_name,
            timeout_seconds=args.timeout,
            model=model,
            agent_mode=args.agent_mode,
            reasoning_effort=reasoning_effort or None,
        )
        pending: List[Mapping[str, Any]] = []
        backend_frames: List[pd.DataFrame] = []
        key_by_identity: Dict[Tuple[str, str], str] = {}
        items_by_key: Dict[str, List[Mapping[str, Any]]] = {}
        for item in items:
            key = score_cache_key(item, context, backend_name, model, reasoning_effort)
            key_by_identity[(str(item["case_id"]), str(item["mode"]))] = key
            items_by_key.setdefault(key, []).append(item)
        for key, equivalent_items in items_by_key.items():
            item = equivalent_items[0]
            path = _cache_path(cache_dir, backend_name, model, reasoning_effort, key)
            cached = None if args.force or not args.resume else _load_cached_row(path, item)
            if cached is None and args.resume and not args.force and not reasoning_effort:
                compat_key = _empty_effort_compat_cache_key(item, context, backend_name, model)
                compat_path = _cache_path(
                    cache_dir, backend_name, model, reasoning_effort, compat_key,
                )
                cached = _load_cached_row(compat_path, item)
                if cached is not None:
                    _write_cache_row(
                        path,
                        cached.iloc[0].to_dict(),
                        backend_name,
                        model,
                        reasoning_effort,
                        key,
                    )
            if cached is None:
                pending.append(item)
            else:
                cached_rows = []
                base_row = cached.iloc[0].to_dict()
                for equivalent in equivalent_items:
                    row = dict(base_row)
                    row["case_id"] = equivalent["case_id"]
                    row["mode"] = equivalent["mode"]
                    row["cache_hit"] = True
                    cached_rows.append(row)
                backend_frames.append(pd.DataFrame(cached_rows))
                cache_hits += len(equivalent_items)
        batches = list(enumerate(_chunk_pending_cases(pending, args.batch_size), start=1))

        def record_completed_batch(batch_index: int, scored: pd.DataFrame) -> None:
            nonlocal new_scores, new_unique_scores
            expanded_rows = []
            score_path = out_dir / "batches" / (
                f"judge_batch_{batch_index:04d}__{_safe_path_part(backend_name)}.csv"
            )
            for _, row in scored.iterrows():
                identity = (str(row["case_id"]), str(row["mode"]))
                key = key_by_identity[identity]
                _write_cache_row(
                    _cache_path(cache_dir, backend_name, model, reasoning_effort, key),
                    row.to_dict(), backend_name, model, reasoning_effort, key,
                )
                for equivalent in items_by_key[key]:
                    expanded = row.to_dict()
                    expanded["case_id"] = equivalent["case_id"]
                    expanded["mode"] = equivalent["mode"]
                    expanded["cache_hit"] = False
                    expanded_rows.append(expanded)
                    new_scores += 1
                new_unique_scores += 1
            expanded_frame = pd.DataFrame(expanded_rows)
            write_scores_csv(score_path, expanded_frame)
            backend_frames.append(expanded_frame)

        if args.max_workers == 1 or len(batches) <= 1:
            for batch_index, batch in batches:
                record_completed_batch(
                    batch_index,
                    _judge_batch(
                        backend, backend_name, model, reasoning_effort, batch,
                        batch_index, out_dir, context, seed, args.max_attempts,
                    ),
                )
        else:
            with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                future_map = {
                    executor.submit(
                        _judge_batch,
                        backend,
                        backend_name,
                        model,
                        reasoning_effort,
                        batch,
                        batch_index,
                        out_dir,
                        context,
                        seed,
                        args.max_attempts,
                    ): batch_index
                    for batch_index, batch in batches
                }
                failed_batches = []
                for future in as_completed(future_map):
                    batch_index = future_map[future]
                    try:
                        record_completed_batch(batch_index, future.result())
                    except Exception as exc:
                        failed_batches.append((batch_index, exc))
                        print(f"judge batch {batch_index} failed: {exc}", flush=True)
                if failed_batches:
                    details = "; ".join(
                        f"{batch_index}: {_compact_text(exc, 200)}"
                        for batch_index, exc in failed_batches
                    )
                    raise RuntimeError(
                        f"{len(failed_batches)} judge batch(es) failed after all other "
                        f"completed batches were cached: {details}"
                    )
        if backend_frames:
            records = [
                record
                for frame in backend_frames
                for record in frame.to_dict(orient="records")
            ]
            all_scores.append(pd.DataFrame(records))

    if not all_scores:
        raise RuntimeError("No judge scores were produced")
    scored = pd.concat(all_scores, ignore_index=True)
    scored = scored.sort_values(["case_id", "mode", "judge_backend"]).reset_index(drop=True)
    summary = summarize_scores(scored, by=("mode",))
    per_judge = summarize_scores(
        scored,
        by=("mode", "judge_backend", "judge_model", "judge_reasoning_effort"),
    )
    consensus = build_consensus_scores(scored)
    write_scores_csv(out_dir / "judge_scores.csv", scored)
    write_scores_csv(out_dir / "judge_summary.csv", summary)
    write_scores_csv(out_dir / "judge_summary_by_judge.csv", per_judge)
    write_scores_csv(out_dir / "judge_consensus.csv", consensus)
    _write_json(out_dir / "judge_summary.json", summary.to_dict(orient="records"))
    html_path = write_judge_html(scored, summary, out_dir)
    manifest.update({
        "status": "completed",
        "updated_at": utc_now(),
        "cache_hits": cache_hits,
        "new_scores": new_scores,
        "new_unique_scores": new_unique_scores,
        "scores_csv": str((out_dir / "judge_scores.csv").resolve()),
        "summary_csv": str((out_dir / "judge_summary.csv").resolve()),
        "consensus_csv": str((out_dir / "judge_consensus.csv").resolve()),
        "html_report": str(html_path.resolve()),
    })
    _write_json(out_dir / "judge_manifest.json", manifest)
    print(f"stable judge complete: {new_scores} new, {cache_hits} cached")
    print(f"wrote {out_dir / 'judge_scores.csv'}")
    return 0


def add_judge_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the stable Judge protocol on a parser and return it."""
    parser.add_argument("--truth", required=True, help="Truth CSV with case_id and expected label.")
    parser.add_argument(
        "--prediction", action="append", required=True, metavar="MODE:CSV",
        help="Prediction CSV; repeat for multiple annotation modes.",
    )
    parser.add_argument("--markers", help="Optional case_id-to-marker-list CSV.")
    parser.add_argument("--out", help="Output directory.")
    parser.add_argument("--backend", action="append", help="Judge backend; repeat for an ensemble.")
    parser.add_argument(
        "--judge-model", action="append",
        help="Pinned model, either MODEL for one backend or BACKEND=MODEL.",
    )
    parser.add_argument(
        "--judge-reasoning-effort",
        action="append",
        metavar="[BACKEND=]EFFORT",
        help=(
            "Optional Codex reasoning effort (for example high or max). "
            "Use BACKEND=EFFORT when multiple judge backends are selected."
        ),
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="Cases per judge call.")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="Concurrent Judge calls. Batch contents and blind labels remain deterministic.",
    )
    parser.add_argument("--blind-seed", type=int, default=17)
    parser.add_argument("--context", default="")
    parser.add_argument("--cache-dir", help="Versioned score-cache directory.")
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--force", action="store_true", help="Ignore and overwrite matching cache entries.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit-cases", type=int)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--agent-mode", default="ask")
    parser.add_argument("--truth-case-column")
    parser.add_argument("--expected-label-column")
    parser.add_argument("--expected-terms-column")
    parser.add_argument("--source-column")
    parser.add_argument("--truth-marker-column")
    parser.add_argument("--marker-case-column")
    parser.add_argument("--marker-column")
    parser.add_argument("--prediction-case-column")
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    add_judge_arguments(parser)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    return run_judge(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
