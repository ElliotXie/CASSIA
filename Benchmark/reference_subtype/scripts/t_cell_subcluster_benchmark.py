"""
T-cell subtype reference benchmark.

Compares Kimi K2.6 CASSIA subclustering output with and without the T-cell
reference agent on the bundled real subcluster marker fixture:
``CASSIA_python/CASSIA/data/subcluster_results.csv``.

The fixture is useful because it is not a clean all-T-cell toy set: it contains
T-cell-like subclusters plus epithelial, NK/myeloid-like, and neural/stromal
contaminating programs from a parent CD8-positive alpha-beta T-cell cluster.
"""

from __future__ import annotations

import os
import sys
import time
import argparse
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = SUITE_DIR / "results"
INPUT_MARKERS = ROOT / "CASSIA_python" / "CASSIA" / "data" / "subcluster_results.csv"
KEY_FILE = SUITE_DIR / ".openrouter_key"

sys.path.insert(0, str(ROOT / "Test" / "shared" / "python"))
from test_utils import setup_cassia_imports  # noqa: E402

setup_cassia_imports()


MODEL = "moonshotai/kimi-k2.6"
PROVIDER = "openrouter"

CASES: List[Dict[str, Any]] = [
    {
        "id": "0",
        "expected_label": "IL7R/RORA memory or helper-like T cell",
        "expected_terms": [["T cell", "T-cell"], ["IL7R", "memory"], ["RORA", "helper", "Th17"]],
        "min_score": 2,
    },
    {
        "id": "1",
        "expected_label": "TRDC/KLRC2 innate-like or gamma-delta/NKT-like T cell",
        "expected_terms": [["TRDC", "gamma-delta", "gamma delta", "gdT"], ["innate-like", "NKT", "NK-like"], ["TIGIT", "ENTPD1", "LAYN", "exhausted", "regulatory"]],
        "min_score": 2,
    },
    {
        "id": "2",
        "expected_label": "GZMK/GZMH cytotoxic effector-memory CD8 T cell",
        "expected_terms": [["CD8", "T cell", "T-cell"], ["GZMK", "GZMH"], ["cytotoxic", "effector", "memory"]],
        "min_score": 2,
    },
    {
        "id": "3",
        "expected_label": "EPCAM/KRT epithelial contaminant",
        "expected_terms": [["epithelial", "epithelium"], ["EPCAM", "KRT", "keratin"], ["contaminant", "non-T", "non T", "not T"]],
        "min_score": 2,
    },
    {
        "id": "4",
        "expected_label": "KLRF1/GNLY/TYROBP NK or NK-like contaminant",
        "expected_terms": [["NK", "natural killer"], ["KLRF1", "GNLY", "NKG7"], ["contaminant", "non-T", "non T", "myeloid", "TYROBP", "FCER1G"]],
        "min_score": 2,
    },
    {
        "id": "5",
        "expected_label": "PLP1/PMP22/HAND2 neural crest or Schwann/stromal contaminant",
        "expected_terms": [["neural", "Schwann", "glial"], ["PLP1", "PMP22", "HAND2"], ["contaminant", "non-T", "non T", "stromal"]],
        "min_score": 2,
    },
]


def load_key_from_local_file() -> None:
    if os.environ.get("OPENROUTER_API_KEY"):
        return
    if KEY_FILE.exists():
        key = KEY_FILE.read_text(encoding="utf-8").strip()
        if key:
            os.environ["OPENROUTER_API_KEY"] = key


def marker_preview(marker_df: pd.DataFrame, cluster_id: str, n: int = 20) -> str:
    rows = marker_df[marker_df["cluster"].astype(str) == str(cluster_id)]
    return ", ".join(rows["gene"].astype(str).head(n).tolist())


def score_text(text: str, expected_terms: List[List[str]]) -> int:
    text_lower = str(text).lower()
    score = 0
    for group in expected_terms:
        if any(term.lower() in text_lower for term in group):
            score += 1
    return score


def combine_output_text(row: pd.Series) -> str:
    return " ".join(str(row.get(col, "")) for col in ["main_cell_type", "sub_cell_type", "reason"])


def run_mode(outdir: Path, mode: str, use_reference: bool) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    from CASSIA import get_llm_usage_summary, reset_llm_usage_log, runCASSIA_subclusters

    output_base = str(outdir / mode)
    reset_llm_usage_log()

    kwargs: Dict[str, Any] = {}
    if use_reference:
        kwargs.update({
            "reference_model": MODEL,
            "reference_cell_type_hint": "t cell",
            "reference_max_content_length": 14000,
            "reference_max_context_length": 14000,
        })

    runCASSIA_subclusters(
        marker=pd.read_csv(INPUT_MARKERS),
        major_cluster_info="CD8-positive, alpha-beta T cell",
        output_name=output_base,
        model=MODEL,
        temperature=0,
        provider=PROVIDER,
        n_genes=40,
        tissue="human intestine",
        species="human",
        use_reference=use_reference,
        **kwargs,
    )
    usage = get_llm_usage_summary(reset=True)
    result = pd.read_csv(f"{output_base}.csv").set_index("Result ID")
    return result, usage


def score_existing(outdir: Path) -> int:
    baseline_path = outdir / "t_cell_baseline.csv"
    reference_path = outdir / "t_cell_reference.csv"
    if not baseline_path.exists() or not reference_path.exists():
        print(f"ERROR: missing benchmark outputs in {outdir}")
        return 2

    baseline = pd.read_csv(baseline_path)
    reference = pd.read_csv(reference_path)
    baseline["Result ID"] = baseline["Result ID"].astype(str)
    reference["Result ID"] = reference["Result ID"].astype(str)
    baseline = baseline.set_index("Result ID")
    reference = reference.set_index("Result ID")

    score_rows = []
    baseline_total = 0
    reference_total = 0
    reference_passed = 0
    improved_or_equal = 0

    for case in CASES:
        case_id = case["id"]
        base_text = combine_output_text(baseline.loc[case_id])
        ref_text = combine_output_text(reference.loc[case_id])
        base_score = score_text(base_text, case["expected_terms"])
        ref_score = score_text(ref_text, case["expected_terms"])
        baseline_total += base_score
        reference_total += ref_score
        reference_passed += int(ref_score >= case["min_score"])
        improved_or_equal += int(ref_score >= base_score)
        score_rows.append({
            "case_id": case_id,
            "expected_label": case["expected_label"],
            "expected_terms": ";".join("/".join(group) for group in case["expected_terms"]),
            "min_score": case["min_score"],
            "baseline_score": base_score,
            "baseline_pass": base_score >= case["min_score"],
            "baseline_text": base_text,
            "reference_score": ref_score,
            "reference_pass": ref_score >= case["min_score"],
            "reference_text": ref_text,
            "improved_or_equal": ref_score >= base_score,
        })

    scores = pd.DataFrame(score_rows)
    scores.to_csv(outdir / "scores.csv", index=False)

    print(f"RESULT_DIR={outdir}")
    print(f"TOTAL baseline={baseline_total} reference={reference_total}")
    print(f"REFERENCE_PASSED={reference_passed}/{len(CASES)}")
    print(f"REFERENCE_IMPROVED_OR_EQUAL={improved_or_equal}/{len(CASES)}")
    print("\nSCORES")
    print(scores[[
        "case_id",
        "expected_label",
        "baseline_score",
        "reference_score",
        "improved_or_equal",
    ]].to_string(index=False))
    print("\nOUTPUT_TEXT")
    for row in score_rows:
        print(f"\n[{row['case_id']}] expected: {row['expected_label']}")
        print(f"baseline({row['baseline_score']}): {row['baseline_text']}")
        print(f"reference({row['reference_score']}): {row['reference_text']}")

    return 0 if reference_passed >= 5 and reference_total >= baseline_total else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-existing", help="Existing result directory to score without making API calls.")
    args = parser.parse_args()
    if args.score_existing:
        return score_existing(Path(args.score_existing))

    load_key_from_local_file()
    if not os.environ.get("OPENROUTER_API_KEY"):
        print("ERROR: OPENROUTER_API_KEY is not set and local key file is missing.")
        return 2
    if not INPUT_MARKERS.exists():
        print(f"ERROR: input marker file not found: {INPUT_MARKERS}")
        return 2

    run_id = time.strftime("%Y%m%d_%H%M%S") + "_tcell"
    outdir = RESULTS_DIR / run_id
    outdir.mkdir(parents=True, exist_ok=True)

    marker_df = pd.read_csv(INPUT_MARKERS)
    case_rows = []
    for case in CASES:
        case_rows.append({
            "case_id": case["id"],
            "expected_label": case["expected_label"],
            "expected_terms": ";".join("/".join(group) for group in case["expected_terms"]),
            "min_score": case["min_score"],
            "markers": marker_preview(marker_df, case["id"]),
        })
    pd.DataFrame(case_rows).to_csv(outdir / "cases.csv", index=False)

    baseline, baseline_usage = run_mode(outdir, "t_cell_baseline", use_reference=False)
    reference, reference_usage = run_mode(outdir, "t_cell_reference", use_reference=True)
    baseline = baseline.reset_index()
    reference = reference.reset_index()
    baseline.to_csv(outdir / "t_cell_baseline.csv", index=False)
    reference.to_csv(outdir / "t_cell_reference.csv", index=False)

    usage = pd.DataFrame([
        {"mode": "t_cell_baseline", **baseline_usage},
        {"mode": "t_cell_reference", **reference_usage},
    ]).drop(columns=["by_model"], errors="ignore")
    usage.to_csv(outdir / "usage.csv", index=False)
    print("\nUSAGE")
    print(usage.to_string(index=False))
    return score_existing(outdir)


if __name__ == "__main__":
    raise SystemExit(main())
