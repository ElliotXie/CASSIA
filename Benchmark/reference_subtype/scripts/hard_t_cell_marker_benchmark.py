"""
Hard T-cell subtype marker benchmark.

This benchmark uses closely related T-cell subtype marker panels adapted from
the pan-cancer tumor-infiltrating T-cell literature and a supplementary marker
table. The goal is harder than the bundled CD8 subcluster fixture: most cases
are all T cells, many share checkpoint, cytotoxic, memory, or IFN genes, and
some subtype pairs are deliberately close.

Source anchors:
- Zheng et al. Science 2021, doi:10.1126/science.abe6474
- Supplementary marker table listing Treg, exhausted CD8, proliferative,
  helper, cytotoxic CD8, pre-exhausted CD8, naive, effector-memory CD8,
  naive-memory CD4, Th17, transitional-memory CD4, and recently activated CD4
  marker panels.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = SUITE_DIR / "results"
KEY_FILE = SUITE_DIR / ".openrouter_key"

sys.path.insert(0, str(ROOT / "Test" / "shared" / "python"))
from test_utils import setup_cassia_imports  # noqa: E402

setup_cassia_imports()


MODEL = "moonshotai/kimi-k2.6"
PROVIDER = "openrouter"

CASES: List[Dict[str, Any]] = [
    {
        "id": "treg",
        "label": "Regulatory T cell",
        "markers": "FOXP3, TNFRSF18, IL2RA, TIGIT, CTLA4, IKZF2",
        "expected_terms": [["regulatory", "Treg"], ["FOXP3"], ["CTLA4", "IL2RA", "TIGIT"]],
        "min_score": 2,
    },
    {
        "id": "terminal_exhausted_cd8",
        "label": "Terminally exhausted CD8 T cell",
        "markers": "CXCL13, LAG3, GZMB, CCL5, NKG7, IFNG, GZMA, HAVCR2, GNLY, PDCD1, TIGIT, TNFRSF9, ENTPD1, CTLA4, PRF1, TOX, GZMH, GZMK",
        "expected_terms": [["exhausted", "dysfunctional"], ["CD8"], ["CXCL13", "TOX", "PDCD1", "HAVCR2", "LAG3"]],
        "min_score": 2,
    },
    {
        "id": "proliferative_t",
        "label": "Proliferative T cell",
        "markers": "STMN1, MKI67, CDK1, TOP2A, UBE2C, PCLAF, CENPF, TYMS, HMGB2, TRAC, CD3D",
        "expected_terms": [["prolifer", "cycling", "cell cycle"], ["MKI67", "TOP2A", "STMN1"], ["T cell", "T-cell"]],
        "min_score": 2,
    },
    {
        "id": "helper_checkpoint_cd4",
        "label": "Checkpoint-high helper T cell",
        "markers": "CXCL13, TNFRSF4, TNFRSF18, BATF, TIGIT, SOX4, TNFRSF25, CTLA4, RORA, XCL1, TNFSF8, PPIA, STAT5A, TOX, PDCD1",
        "expected_terms": [["helper", "CD4", "Tfh"], ["checkpoint", "activated", "PDCD1", "TIGIT", "CTLA4"], ["CXCL13", "TNFRSF4", "BATF"]],
        "min_score": 2,
    },
    {
        "id": "cytotoxic_cd8",
        "label": "Cytotoxic CD8 T cell",
        "markers": "CCL5, GZMA, GZMK, NKG7, GNLY, GZMH, GZMM, PRF1, CXCR3, GZMB, CCL4, IFNG",
        "expected_terms": [["cytotoxic", "effector"], ["CD8"], ["GZMK", "GZMB", "PRF1", "GNLY", "NKG7"]],
        "min_score": 2,
    },
    {
        "id": "pre_exhausted_isg_cd8",
        "label": "Pre-exhausted / IFN-stimulated CD8 T cell",
        "markers": "ISG15, IFI44L, IFI6, IFIT3, IFIT1, IFI44, IFI35, IRF7, IFIT2, LAG3, IFITM1, IFI16, IFI27, IFNG, GZMB, GZMK, PRF1, HAVCR2, IFIH1, GZMA, IRF9, CXCL13, GZMH, IFIT5, PDCD1",
        "expected_terms": [["pre-exhausted", "precursor", "early exhausted", "exhausted"], ["interferon", "IFN", "ISG"], ["CD8", "cytotoxic"]],
        "min_score": 2,
    },
    {
        "id": "naive_t",
        "label": "Naive T cell",
        "markers": "IL7R, TCF7, CCR7, LEF1, SELL, LTB, MAL, NOSIP, PIK3IP1, TRAC, CD3D",
        "expected_terms": [["naive"], ["TCF7", "CCR7", "LEF1", "SELL"], ["T cell", "T-cell"]],
        "min_score": 2,
    },
    {
        "id": "effector_memory_cd8",
        "label": "Effector memory CD8 T cell",
        "markers": "GZMM, IFITM1, GZMK, IFNG, CCL5, CXCR3, CST7, NKG7, GZMA, TRAC, CD8A",
        "expected_terms": [["effector memory", "effector-memory", "memory"], ["CD8"], ["GZMK", "GZMM", "CCL5"]],
        "min_score": 2,
    },
    {
        "id": "naive_memory_cd4",
        "label": "Naive-memory CD4 T cell",
        "markers": "IL7R, TCF7, CCL5, IFITM1, CCR7, LEF1, CD40LG, ANXA1, TRAC, CD3D",
        "expected_terms": [["naive", "memory"], ["CD4", "helper"], ["IL7R", "TCF7", "CCR7"]],
        "min_score": 2,
    },
    {
        "id": "th17",
        "label": "Th17 cell",
        "markers": "IL17A, IL17F, BATF, IL2RA, DUSP4, IL21R, CTLA4, IRF4, CCL20, IL26, BATF3, IL4R, STAT3, RORA",
        "expected_terms": [["Th17", "T helper 17"], ["IL17A", "IL17F", "IL26"], ["CD4", "helper"]],
        "min_score": 2,
    },
    {
        "id": "transitional_memory_cd4",
        "label": "Transitional memory CD4 T cell",
        "markers": "CXCL13, TNFRSF4, TIGIT, IL6ST, PASK, KLRB1, CD40LG, TOX, LEF1, ICOS, CD28, TOX2, CCR7, CD247, RORA, PDCD1",
        "expected_terms": [["transitional", "memory"], ["CD4", "helper"], ["CXCL13", "ICOS", "CD40LG", "PDCD1"]],
        "min_score": 2,
    },
    {
        "id": "recently_activated_cd4",
        "label": "Recently activated CD4 T cell",
        "markers": "CCL4, IFITM1, CD69, PRF1, BCL3, IL7R, IFITM3, TCF7, CD81, CXCR4, GZMK, GZMM, IFITM2",
        "expected_terms": [["activated", "recently activated"], ["CD4", "helper"], ["CD69", "BCL3", "IL7R", "TCF7"]],
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


def marker_dataframe() -> pd.DataFrame:
    return pd.DataFrame({
        "cluster": [case["id"] for case in CASES],
        "markers": [case["markers"] for case in CASES],
    })


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
            "reference_max_content_length": 18000,
            "reference_max_context_length": 18000,
        })

    runCASSIA_subclusters(
        marker=marker_dataframe(),
        major_cluster_info="human pan-cancer tumor-infiltrating T cell",
        output_name=output_base,
        model=MODEL,
        temperature=0,
        provider=PROVIDER,
        n_genes=30,
        tissue="pan-cancer tumor",
        species="human",
        use_reference=use_reference,
        **kwargs,
    )
    usage = get_llm_usage_summary(reset=True)
    result = pd.read_csv(f"{output_base}.csv")
    result["Result ID"] = result["Result ID"].astype(str)
    return result.set_index("Result ID"), usage


def score_existing(outdir: Path) -> int:
    baseline_path = outdir / "hard_t_cell_baseline.csv"
    reference_path = outdir / "hard_t_cell_reference.csv"
    if not baseline_path.exists() or not reference_path.exists():
        print(f"ERROR: missing benchmark outputs in {outdir}")
        return 2

    baseline = pd.read_csv(baseline_path)
    reference = pd.read_csv(reference_path)
    baseline["Result ID"] = baseline["Result ID"].astype(str)
    reference["Result ID"] = reference["Result ID"].astype(str)
    baseline = baseline.set_index("Result ID")
    reference = reference.set_index("Result ID")

    rows = []
    baseline_total = 0
    reference_total = 0
    reference_passed = 0
    improved_or_equal = 0
    for case in CASES:
        base_missing = case["id"] not in baseline.index
        ref_missing = case["id"] not in reference.index
        base_text = "[MISSING OUTPUT ROW]" if base_missing else combine_output_text(baseline.loc[case["id"]])
        ref_text = "[MISSING OUTPUT ROW]" if ref_missing else combine_output_text(reference.loc[case["id"]])
        base_score = score_text(base_text, case["expected_terms"])
        ref_score = score_text(ref_text, case["expected_terms"])
        baseline_total += base_score
        reference_total += ref_score
        reference_passed += int(ref_score >= case["min_score"])
        improved_or_equal += int(ref_score >= base_score)
        rows.append({
            "case_id": case["id"],
            "expected_label": case["label"],
            "expected_terms": ";".join("/".join(group) for group in case["expected_terms"]),
            "min_score": case["min_score"],
            "baseline_score": base_score,
            "baseline_pass": base_score >= case["min_score"],
            "baseline_missing": base_missing,
            "baseline_text": base_text,
            "reference_score": ref_score,
            "reference_pass": ref_score >= case["min_score"],
            "reference_missing": ref_missing,
            "reference_text": ref_text,
            "improved_or_equal": ref_score >= base_score,
        })

    scores = pd.DataFrame(rows)
    scores.to_csv(outdir / "scores.csv", index=False)

    print(f"RESULT_DIR={outdir}")
    print(f"TOTAL baseline={baseline_total} reference={reference_total}")
    print(f"REFERENCE_PASSED={reference_passed}/{len(CASES)}")
    print(f"REFERENCE_IMPROVED_OR_EQUAL={improved_or_equal}/{len(CASES)}")
    print("\nSCORES")
    print(scores[["case_id", "expected_label", "baseline_score", "reference_score", "improved_or_equal"]].to_string(index=False))
    print("\nOUTPUT_TEXT")
    for row in rows:
        print(f"\n[{row['case_id']}] expected: {row['expected_label']}")
        print(f"baseline({row['baseline_score']}): {row['baseline_text']}")
        print(f"reference({row['reference_score']}): {row['reference_text']}")

    return 0 if reference_passed >= 9 and reference_total >= baseline_total else 1


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

    run_id = time.strftime("%Y%m%d_%H%M%S") + "_hard_tcell"
    outdir = RESULTS_DIR / run_id
    outdir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([
        {
            "case_id": case["id"],
            "expected_label": case["label"],
            "expected_terms": ";".join("/".join(group) for group in case["expected_terms"]),
            "min_score": case["min_score"],
            "markers": case["markers"],
        }
        for case in CASES
    ]).to_csv(outdir / "cases.csv", index=False)

    baseline, baseline_usage = run_mode(outdir, "hard_t_cell_baseline", use_reference=False)
    reference, reference_usage = run_mode(outdir, "hard_t_cell_reference", use_reference=True)

    baseline.reset_index().to_csv(outdir / "hard_t_cell_baseline.csv", index=False)
    reference.reset_index().to_csv(outdir / "hard_t_cell_reference.csv", index=False)
    usage = pd.DataFrame([
        {"mode": "hard_t_cell_baseline", **baseline_usage},
        {"mode": "hard_t_cell_reference", **reference_usage},
    ]).drop(columns=["by_model"], errors="ignore")
    usage.to_csv(outdir / "usage.csv", index=False)

    print("\nUSAGE")
    print(usage.to_string(index=False))
    return score_existing(outdir)


if __name__ == "__main__":
    raise SystemExit(main())
