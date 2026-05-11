"""Held-out Zhang CRC T-cell subtype benchmark.

This benchmark uses Zhang et al. Nature 2018 Supplementary Table 5 signature
markers. The model input intentionally hides the paper's CD4/CD8 cluster names:
only neutral cluster IDs and marker genes are passed to CASSIA.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
HELDOUT_DIR = SUITE_DIR / "heldout" / "zhang_crc_2018"
RESULTS_DIR = SUITE_DIR / "results"
KEY_FILE = SUITE_DIR / ".openrouter_key"
INPUT_CSV = HELDOUT_DIR / "zhang_crc_tcell_marker_inputs.csv"
GROUND_TRUTH_CSV = HELDOUT_DIR / "zhang_crc_tcell_ground_truth.csv"

sys.path.insert(0, str(ROOT / "Test" / "shared" / "python"))
from test_utils import setup_cassia_imports  # noqa: E402

setup_cassia_imports()

MODEL = "moonshotai/kimi-k2.6"
PROVIDER = "openrouter"


def load_key_from_local_file() -> None:
    if os.environ.get("OPENROUTER_API_KEY"):
        return
    if KEY_FILE.exists():
        key = KEY_FILE.read_text(encoding="utf-8").strip()
        if key:
            os.environ["OPENROUTER_API_KEY"] = key


def load_cases() -> pd.DataFrame:
    inputs = pd.read_csv(INPUT_CSV)
    truth = pd.read_csv(GROUND_TRUTH_CSV)
    cases = inputs.merge(truth, on="cluster_id", how="inner")
    if len(cases) != len(inputs) or len(cases) != len(truth):
        raise RuntimeError("held-out input/ground-truth row mismatch")
    return cases


def marker_dataframe() -> pd.DataFrame:
    cases = load_cases()
    return pd.DataFrame({
        "cluster": cases["cluster_id"].astype(str),
        "markers": cases["marker_genes"].astype(str),
    })


def parse_expected_terms(value: str) -> List[List[str]]:
    groups: List[List[str]] = []
    for group in str(value).split(";"):
        terms = [term.strip() for term in group.split("/") if term.strip()]
        if terms:
            groups.append(terms)
    return groups


def score_text(text: str, expected_terms: List[List[str]]) -> int:
    text_lower = str(text).lower()
    score = 0
    for group in expected_terms:
        if any(term.lower() in text_lower for term in group):
            score += 1
    return score


def combine_output_text(row: pd.Series) -> str:
    cols = ["main_cell_type", "sub_cell_type", "reason", "cell_type", "annotation", "evidence"]
    return " ".join(str(row.get(col, "")) for col in cols)


def run_mode(outdir: Path, mode: str, use_reference: bool) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    from CASSIA import get_llm_usage_summary, reset_llm_usage_log, runCASSIA_subclusters

    output_base = str(outdir / f"zhang_crc_{mode}")
    reset_llm_usage_log()
    kwargs: Dict[str, Any] = {}
    if use_reference:
        kwargs.update({
            "reference_model": MODEL,
            "reference_cell_type_hint": "t cell",
            "reference_max_content_length": 18000,
            "reference_max_context_length": 18000,
        })

    runCASSIA_subclusters(
        marker=marker_dataframe(),
        major_cluster_info="human colorectal cancer T cells; CD4/CD8 paper labels are hidden from input",
        output_name=output_base,
        model=MODEL,
        temperature=0,
        provider=PROVIDER,
        n_genes=30,
        tissue="colorectal cancer tumor and adjacent tissue",
        species="human",
        use_reference=use_reference,
        **kwargs,
    )
    usage = get_llm_usage_summary(reset=True)
    result = pd.read_csv(f"{output_base}.csv")
    result["Result ID"] = result["Result ID"].astype(str)
    return result.set_index("Result ID"), usage


def score_outputs(outdir: Path) -> pd.DataFrame:
    baseline_path = outdir / "zhang_crc_baseline.csv"
    reference_path = outdir / "zhang_crc_reference.csv"
    if not baseline_path.exists() or not reference_path.exists():
        raise FileNotFoundError(f"missing outputs: {baseline_path} or {reference_path}")

    baseline = pd.read_csv(baseline_path)
    reference = pd.read_csv(reference_path)
    baseline["Result ID"] = baseline["Result ID"].astype(str)
    reference["Result ID"] = reference["Result ID"].astype(str)
    baseline = baseline.set_index("Result ID")
    reference = reference.set_index("Result ID")

    rows = []
    cases = load_cases()
    for _, case in cases.iterrows():
        case_id = str(case["cluster_id"])
        expected_terms = parse_expected_terms(case["expected_terms"])
        base_missing = case_id not in baseline.index
        ref_missing = case_id not in reference.index
        base_text = "[MISSING OUTPUT ROW]" if base_missing else combine_output_text(baseline.loc[case_id])
        ref_text = "[MISSING OUTPUT ROW]" if ref_missing else combine_output_text(reference.loc[case_id])
        base_score = score_text(base_text, expected_terms)
        ref_score = score_text(ref_text, expected_terms)
        rows.append({
            "case_id": case_id,
            "source_sheet": case["source_sheet"],
            "expected_label": case["expected_label"],
            "expected_terms": case["expected_terms"],
            "baseline_score": base_score,
            "baseline_pass": base_score >= 2,
            "baseline_missing": base_missing,
            "baseline_text": base_text,
            "reference_score": ref_score,
            "reference_pass": ref_score >= 2,
            "reference_missing": ref_missing,
            "reference_text": ref_text,
            "reference_minus_baseline": ref_score - base_score,
        })

    scores = pd.DataFrame(rows)
    scores.to_csv(outdir / "heldout_scores.csv", index=False)
    write_summary(outdir, scores)
    return scores


def write_summary(outdir: Path, scores: pd.DataFrame) -> None:
    baseline_total = int(scores["baseline_score"].sum())
    reference_total = int(scores["reference_score"].sum())
    baseline_pass = int(scores["baseline_pass"].sum())
    reference_pass = int(scores["reference_pass"].sum())
    worsened = scores[scores["reference_minus_baseline"] < 0]

    lines = [
        "# Zhang CRC held-out T-cell benchmark evaluation",
        "",
        f"- Cases: {len(scores)}",
        f"- Baseline total score: {baseline_total}",
        f"- Reference total score: {reference_total}",
        f"- Baseline pass count: {baseline_pass}/{len(scores)}",
        f"- Reference pass count: {reference_pass}/{len(scores)}",
        f"- Reference worsened cases: {len(worsened)}",
        "",
        "## Per-case scores",
        "",
    ]
    for _, row in scores.iterrows():
        lines.extend([
            f"### {row['case_id']} ({row['source_sheet']})",
            f"- Expected: {row['expected_label']}",
            f"- Baseline score/pass: {row['baseline_score']} / {row['baseline_pass']}",
            f"- Reference score/pass: {row['reference_score']} / {row['reference_pass']}",
            f"- Delta: {row['reference_minus_baseline']}",
            "",
        ])
    (outdir / "manual_evaluation.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-existing", type=Path, help="Score an existing result directory without new LLM calls")
    parser.add_argument("--outdir", type=Path, help="Output directory for a new run")
    args = parser.parse_args()

    load_key_from_local_file()
    if args.score_existing:
        scores = score_outputs(args.score_existing)
        print(json.dumps({
            "outdir": str(args.score_existing),
            "baseline_total": int(scores["baseline_score"].sum()),
            "reference_total": int(scores["reference_score"].sum()),
            "baseline_pass": int(scores["baseline_pass"].sum()),
            "reference_pass": int(scores["reference_pass"].sum()),
        }, indent=2))
        return 0

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    outdir = args.outdir or RESULTS_DIR / f"{timestamp}_heldout_zhang_crc_tcell"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "metadata.json").write_text(json.dumps({
        "dataset": "Zhang CRC T-cell held-out",
        "input_csv": str(INPUT_CSV),
        "ground_truth_csv": str(GROUND_TRUTH_CSV),
        "model": MODEL,
        "provider": PROVIDER,
    }, indent=2), encoding="utf-8")

    run_mode(outdir, "baseline", use_reference=False)
    run_mode(outdir, "reference", use_reference=True)
    scores = score_outputs(outdir)
    print(json.dumps({
        "outdir": str(outdir),
        "baseline_total": int(scores["baseline_score"].sum()),
        "reference_total": int(scores["reference_score"].sum()),
        "baseline_pass": int(scores["baseline_pass"].sum()),
        "reference_pass": int(scores["reference_pass"].sum()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
