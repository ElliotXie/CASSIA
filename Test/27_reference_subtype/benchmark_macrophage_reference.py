"""
Optional live benchmark for macrophage reference-assisted subclustering.

This script compares macrophage subtype calls with and without reference mode
using OpenRouter.

Evaluation standard:
- Each case has expected literature subtype term groups.
- A row score is the count of expected term groups found in main_cell_type,
  sub_cell_type, or reason. Each group can include synonyms, such as
  ["IFNG", "IFN-γ", "IFN-gamma"].
- A reference call passes a case when it reaches the case's min_score.
- The benchmark passes when reference mode reaches the pass threshold and does
  not underperform baseline on total score.

It intentionally skips when OPENROUTER_API_KEY is not set so normal test runs
do not require live API calls.
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import setup_cassia_imports

setup_cassia_imports()


MODEL = "moonshotai/kimi-k2.6"
PROVIDER = "openrouter"

CASES = [
    {
        "id": "spp1_areg_tam",
        "markers": "SPP1, AREG, EREG, CXCL8, CXCL3, IL1B, TIMP1, LYZ, C1QA, TYROBP",
        "expected_terms": [["SPP1"], ["AREG"], ["TAM", "tumor-associated macrophage"]],
        "min_score": 2,
    },
    {
        "id": "ifng_tam",
        "markers": "CXCL9, CXCL10, GBP1, GBP5, STAT1, WARS1, MMP9, LYZ, C1QA, FCER1G",
        "expected_terms": [["IFNG", "IFN-γ", "IFN-gamma"], ["CXCL9"], ["interferon", "IFN"]],
        "min_score": 2,
    },
    {
        "id": "resident_folr2",
        "markers": "FOLR2, SELENOP, SLC40A1, STAB1, C1QA, C1QB, C1QC, CD163, MRC1, LYZ",
        "expected_terms": [["FOLR2"], ["resident"], ["macrophage", "TAM"]],
        "min_score": 2,
    },
    {
        "id": "apoe_trem2_lipid_tam",
        "markers": "TREM2, APOE, APOC1, GPNMB, LIPA, CTSD, LGMN, PLA2G7, ACP5, LYZ, C1QA",
        "expected_terms": [["APOE"], ["TREM2"], ["lipid", "LAM"]],
        "min_score": 2,
    },
    {
        "id": "c1qc_phagocytic",
        "markers": "C1QA, C1QB, C1QC, APOE, APOC1, CD74, HLA-DRA, HLA-DPA1, PLD4, GPR34, LYZ",
        "expected_terms": [["C1QC", "C1Q"], ["phagocytic", "phagocytosis"], ["macrophage", "TAM"]],
        "min_score": 2,
    },
    {
        "id": "isg15_type_i_ifn",
        "markers": "ISG15, IFIT1, IFIT2, IFIT3, IFITM1, IFITM3, MX1, OAS1, STAT1, LYZ, C1QA",
        "expected_terms": [["ISG15", "ISG"], ["interferon", "IFN"], ["macrophage", "TAM"]],
        "min_score": 2,
    },
    {
        "id": "il1b_tnf_inflammatory",
        "markers": "IL1B, TNF, CXCL8, CXCL1, CXCL2, CCL3, CCL4, NFKBIA, LYZ, FCER1G, C1QA",
        "expected_terms": [["IL1B"], ["TNF"], ["inflammatory", "inflammation"]],
        "min_score": 2,
    },
    {
        "id": "nlrp3_inflammasome",
        "markers": "NLRP3, IL1B, CASP1, CASP4, CXCL8, TNF, NFKBIA, LYZ, C1QA, TYROBP",
        "expected_terms": [["NLRP3"], ["inflammasome"], ["macrophage", "TAM"]],
        "min_score": 2,
    },
    {
        "id": "metallothionein_macrophage",
        "markers": "MT1G, MT1X, MT2A, MT1E, MT1H, MT1F, MIF, SPP1, LDHA, LGALS1, LYZ",
        "expected_terms": [["metallothionein", "MetalloMac"], ["stress", "metal"], ["macrophage", "TAM"]],
        "min_score": 2,
    },
    {
        "id": "heme_iron_macrophage",
        "markers": "HMOX1, SLC40A1, HAMP, CD163, CCL18, LGMN, FTL, FTH1, LYZ, C1QA, APOE",
        "expected_terms": [["heme", "HMOX1"], ["iron"], ["macrophage", "TAM"]],
        "min_score": 2,
    },
    {
        "id": "ecm_remodeling_tam",
        "markers": "COL1A1, COL1A2, COL3A1, SPARC, COL6A1, COL6A2, MMP2, MMP14, LYZ, C1QA, TYROBP",
        "expected_terms": [["ECM", "matrix"], ["remodeling", "remodelling"], ["TAM", "macrophage"]],
        "min_score": 2,
    },
    {
        "id": "heat_shock_stress",
        "markers": "HSPA6, HSPA1A, HSPA1B, DNAJB1, HSPB1, BAG3, FOS, JUN, LYZ, C1QA",
        "expected_terms": [["heat", "HSP"], ["stress"], ["macrophage", "TAM"]],
        "min_score": 2,
    },
]


def _score_row(row, expected_terms):
    text = " ".join(str(row.get(col, "")) for col in ["main_cell_type", "sub_cell_type", "reason"])
    text_lower = text.lower()
    score = 0
    for term_group in expected_terms:
        synonyms = term_group if isinstance(term_group, (list, tuple)) else [term_group]
        if any(str(term).lower() in text_lower for term in synonyms):
            score += 1
    return score, text


def _print_usage_summary(name, summary):
    print(
        f"{name} usage: requests={summary['requests']} "
        f"prompt_tokens={summary['prompt_tokens']} "
        f"completion_tokens={summary['completion_tokens']} "
        f"reasoning_tokens={summary['reasoning_tokens']} "
        f"total_tokens={summary['total_tokens']} "
        f"cost={summary['cost']:.8f}"
    )


def main():
    if not os.environ.get("OPENROUTER_API_KEY"):
        print("SKIP: OPENROUTER_API_KEY is not set.")
        return 0

    from CASSIA import (
        runCASSIA_subclusters,
        reset_llm_usage_log,
        get_llm_usage_summary,
    )

    marker_df = pd.DataFrame({
        "cluster": [case["id"] for case in CASES],
        "markers": [case["markers"] for case in CASES],
    })

    with tempfile.TemporaryDirectory(prefix="cassia_ref_benchmark_") as tmpdir:
        baseline_output = str(Path(tmpdir) / "macrophage_subclusters_baseline")
        reference_output = str(Path(tmpdir) / "macrophage_subclusters_reference")

        reset_llm_usage_log()
        runCASSIA_subclusters(
            marker=marker_df,
            major_cluster_info="human tumor macrophage",
            output_name=baseline_output,
            model=MODEL,
            temperature=0,
            provider=PROVIDER,
            n_genes=20,
            tissue="tumor",
            species="human",
            use_reference=False,
        )
        baseline_usage = get_llm_usage_summary(reset=True)

        runCASSIA_subclusters(
            marker=marker_df,
            major_cluster_info="human tumor macrophage",
            output_name=reference_output,
            model=MODEL,
            temperature=0,
            provider=PROVIDER,
            n_genes=20,
            tissue="tumor",
            species="human",
            use_reference=True,
            reference_model=MODEL,
            reference_cell_type_hint="macrophage",
        )
        reference_usage = get_llm_usage_summary(reset=True)

        baseline = pd.read_csv(f"{baseline_output}.csv").set_index("Result ID")
        reference = pd.read_csv(f"{reference_output}.csv").set_index("Result ID")

        rows = []
        improved_or_equal = 0
        reference_passed = 0
        baseline_total = 0
        reference_total = 0
        for case in CASES:
            base_score, base_text = _score_row(baseline.loc[case["id"]], case["expected_terms"])
            ref_score, ref_text = _score_row(reference.loc[case["id"]], case["expected_terms"])
            baseline_total += base_score
            reference_total += ref_score
            if ref_score >= base_score:
                improved_or_equal += 1
            if ref_score >= case.get("min_score", len(case["expected_terms"])):
                reference_passed += 1
            print(f"{case['id']}: baseline_score={base_score} reference_score={ref_score}")
            print(f"  baseline: {base_text}")
            print(f"  reference: {ref_text}")
            rows.append({
                "case_id": case["id"],
                "expected_terms": ";".join(
                    "/".join(group) if isinstance(group, (list, tuple)) else str(group)
                    for group in case["expected_terms"]
                ),
                "min_score": case.get("min_score", len(case["expected_terms"])),
                "baseline_score": base_score,
                "reference_score": ref_score,
                "reference_passed": ref_score >= case.get("min_score", len(case["expected_terms"])),
                "improved_or_equal": ref_score >= base_score,
            })

        score_path = Path(tmpdir) / "macrophage_reference_benchmark_scores.csv"
        usage_path = Path(tmpdir) / "macrophage_reference_benchmark_usage.csv"
        pd.DataFrame(rows).to_csv(score_path, index=False)
        pd.DataFrame([
            {"run": "baseline", **baseline_usage},
            {"run": "reference", **reference_usage},
        ]).drop(columns=["by_model"], errors="ignore").to_csv(usage_path, index=False)

        print(f"Scores written to {score_path}")
        print(f"Usage summary written to {usage_path}")
        _print_usage_summary("Baseline", baseline_usage)
        _print_usage_summary("Reference", reference_usage)

    pass_threshold = int(len(CASES) * 0.8)
    print(f"Reference mode improved or matched baseline in {improved_or_equal}/{len(CASES)} cases.")
    print(f"Reference mode passed expected-term threshold in {reference_passed}/{len(CASES)} cases.")
    print(f"Total score: baseline={baseline_total} reference={reference_total}")
    return 0 if reference_passed >= pass_threshold and reference_total >= baseline_total else 1


if __name__ == "__main__":
    sys.exit(main())
