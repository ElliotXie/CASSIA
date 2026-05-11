"""
Macrophage subtype reference benchmark.

Compares:
  1. CASSIA subclustering baseline without reference
  2. CASSIA subclustering with agentic subtype reference
  3. Direct one-shot Kimi annotation without CASSIA pipeline/reference

The direct Kimi mode deliberately avoids CASSIA's annotation prompt and
reference agent. It is a raw model baseline for the same marker panel.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = Path(__file__).resolve().parent / "results"

sys.path.insert(0, str(ROOT / "Test" / "shared" / "python"))
from test_utils import setup_cassia_imports  # noqa: E402

setup_cassia_imports()


MODEL = "moonshotai/kimi-k2.6"
PROVIDER = "openrouter"

DEFAULT_CASES = [
    {
        "id": "spp1_areg_tam",
        "markers": "SPP1, AREG, EREG, CXCL8, CXCL3, IL1B, TIMP1, LYZ, C1QA, TYROBP",
        "expected_terms": [["SPP1"], ["AREG"], ["TAM", "tumor-associated macrophage"]],
        "min_score": 2,
    },
    {
        "id": "ifng_tam",
        "markers": "CXCL9, CXCL10, GBP1, GBP5, STAT1, WARS1, MMP9, LYZ, C1QA, FCER1G",
        "expected_terms": [["IFNG", "IFN-gamma", "IFN-γ"], ["CXCL9"], ["interferon", "IFN"]],
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


def parse_expected_terms(value: Any) -> List[List[str]]:
    if isinstance(value, list):
        return [
            [str(term).strip() for term in group if str(term).strip()]
            if isinstance(group, list)
            else [str(group).strip()]
            for group in value
        ]
    return [
        [term.strip() for term in group.split("/") if term.strip()]
        for group in str(value or "").split(";")
        if group.strip()
    ]


def load_cases(cases_csv: Optional[str] = None) -> List[Dict[str, Any]]:
    if cases_csv is None:
        raw_cases = DEFAULT_CASES
    else:
        with open(cases_csv, newline="", encoding="utf-8") as handle:
            raw_cases = list(csv.DictReader(handle))

    cases = []
    for raw_case in raw_cases:
        case_id = raw_case.get("id") or raw_case.get("case_id")
        if not case_id:
            raise ValueError(f"Benchmark case is missing id/case_id: {raw_case}")
        markers = raw_case.get("markers")
        if not markers:
            raise ValueError(f"Benchmark case {case_id} is missing markers")
        cases.append({
            **raw_case,
            "id": str(case_id),
            "markers": str(markers),
            "expected_terms": parse_expected_terms(raw_case.get("expected_terms")),
            "min_score": int(raw_case.get("min_score", 1)),
        })
    return cases


def marker_dataframe(cases: List[Dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame({
        "cluster": [case["id"] for case in cases],
        "markers": [case["markers"] for case in cases],
    })


def score_row(row: pd.Series, expected_terms: List[List[str]]) -> Tuple[int, str]:
    text = " ".join(str(row.get(col, "")) for col in ["main_cell_type", "sub_cell_type", "reason"])
    text_lower = text.lower()
    score = 0
    for term_group in expected_terms:
        synonyms = term_group if isinstance(term_group, (list, tuple)) else [term_group]
        if any(str(term).lower() in text_lower for term in synonyms):
            score += 1
    return score, text


def paper_label_hit(text: str, paper_subtype: Any) -> bool:
    subtype = str(paper_subtype or "").strip()
    if not subtype:
        return False
    label = subtype.split("_", 1)[1] if "_" in subtype else subtype
    text_lower = str(text).lower()
    return subtype.lower() in text_lower or label.lower() in text_lower


def usage_summary_from_openrouter(response_json: Dict[str, Any]) -> Dict[str, Any]:
    usage = response_json.get("usage") or {}
    completion_details = usage.get("completion_tokens_details") or {}
    prompt_details = usage.get("prompt_tokens_details") or {}
    return {
        "requests": 1,
        "prompt_tokens": usage.get("prompt_tokens") or 0,
        "completion_tokens": usage.get("completion_tokens") or 0,
        "reasoning_tokens": completion_details.get("reasoning_tokens") or usage.get("reasoning_tokens") or 0,
        "cached_tokens": prompt_details.get("cached_tokens") or 0,
        "total_tokens": usage.get("total_tokens") or 0,
        "cost": usage.get("cost") or 0.0,
        "response_id": response_json.get("id"),
        "model": response_json.get("model"),
    }


def run_cassia_mode(
    outdir: Path,
    mode_name: str,
    use_reference: bool,
    model: str,
    provider: str,
    cases: List[Dict[str, Any]],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    from CASSIA import get_llm_usage_summary, reset_llm_usage_log, runCASSIA_subclusters

    output_base = str(outdir / mode_name)
    reset_llm_usage_log()
    kwargs = {}
    if use_reference:
        kwargs.update({
            "reference_model": model,
            "reference_cell_type_hint": "macrophage",
        })

    runCASSIA_subclusters(
        marker=marker_dataframe(cases),
        major_cluster_info="human tumor macrophage",
        output_name=output_base,
        model=model,
        temperature=0,
        provider=provider,
        n_genes=20,
        tissue="tumor",
        species="human",
        use_reference=use_reference,
        **kwargs,
    )
    usage = get_llm_usage_summary(reset=True)
    result = pd.read_csv(f"{output_base}.csv")
    return result, usage


def direct_kimi_prompt(cases: List[Dict[str, Any]]) -> str:
    clusters = "\n".join(
        f"- Cluster {case['id']}: {case['markers']}"
        for case in cases
    )
    first_case_id = cases[0]["id"] if cases else "cluster_1"
    return f"""You are a single-cell RNA-seq macrophage subtype annotation expert.

Annotate each subcluster independently from the marker genes. Do not use CASSIA
or any external reference tool. Use concise literature-style macrophage subtype
labels where possible.

Parent cluster: human tumor macrophage
Tissue: tumor
Species: human

Subcluster markers:
{clusters}

Return JSON only with this schema:
{{
  "results": [
    {{
      "cluster_id": "{first_case_id}",
      "main_cell_type": "Macrophage",
      "sub_cell_type": "specific subtype label",
      "reason": "brief marker-based reason"
    }}
  ]
}}"""


def parse_direct_json(text: str) -> List[Dict[str, str]]:
    match = re.search(r"\{[\s\S]*\}", text or "")
    if not match:
        raise RuntimeError(f"Direct Kimi returned non-JSON output: {text[:300]}")
    parsed = json.loads(match.group())
    results = parsed.get("results", [])
    if not isinstance(results, list):
        raise RuntimeError("Direct Kimi JSON did not contain a results list")
    return results


def run_direct_kimi(outdir: Path, model: str, cases: List[Dict[str, Any]]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        raise RuntimeError("OPENROUTER_API_KEY is required for direct Kimi benchmark")

    payload = {
        "model": model,
        "messages": [{"role": "user", "content": direct_kimi_prompt(cases)}],
        "temperature": 0,
        "max_tokens": 5000,
        "reasoning": {"effort": "none", "exclude": True},
    }
    response = requests.post(
        "https://openrouter.ai/api/v1/chat/completions",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        data=json.dumps(payload),
        timeout=180,
    )
    response.raise_for_status()
    response_json = response.json()
    message = response_json["choices"][0]["message"]
    content = message.get("content") or message.get("reasoning") or ""
    rows = []
    marker_map = {case["id"]: case["markers"] for case in cases}
    for result in parse_direct_json(content):
        cluster_id = str(result.get("cluster_id", ""))
        rows.append({
            "Result ID": cluster_id,
            "main_cell_type": result.get("main_cell_type", "Unknown"),
            "sub_cell_type": result.get("sub_cell_type", "Unknown"),
            "key_markers": marker_map.get(cluster_id, ""),
            "reason": result.get("reason", ""),
        })

    df = pd.DataFrame(rows, columns=["Result ID", "main_cell_type", "sub_cell_type", "key_markers", "reason"])
    df.to_csv(outdir / "direct_kimi.csv", index=False)
    (outdir / "direct_kimi_raw.json").write_text(json.dumps(response_json, indent=2), encoding="utf-8")
    return df, usage_summary_from_openrouter(response_json)


def expected_terms_to_string(expected_terms: List[List[str]]) -> str:
    return ";".join(
        "/".join(group) if isinstance(group, (list, tuple)) else str(group)
        for group in expected_terms
    )


def evaluate_modes(
    results: Dict[str, pd.DataFrame],
    outdir: Path,
    cases: List[Dict[str, Any]],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    indexed = {name: df.set_index("Result ID") for name, df in results.items()}
    rows = []
    totals = {name: 0 for name in results}
    passed = {name: 0 for name in results}
    label_hits = {name: 0 for name in results}

    for case in cases:
        row = {
            "case_id": case["id"],
            "paper_subtype": case.get("paper_subtype", ""),
            "expected_terms": expected_terms_to_string(case["expected_terms"]),
            "min_score": case["min_score"],
        }
        for mode_name, df in indexed.items():
            if case["id"] in df.index:
                score, text = score_row(df.loc[case["id"]], case["expected_terms"])
            else:
                score, text = 0, "MISSING_RESULT"
            row[f"{mode_name}_score"] = score
            row[f"{mode_name}_pass"] = score >= case["min_score"]
            label_hit = paper_label_hit(text, case.get("paper_subtype"))
            row[f"{mode_name}_paper_label_hit"] = label_hit
            row[f"{mode_name}_text"] = text
            totals[mode_name] += score
            if score >= case["min_score"]:
                passed[mode_name] += 1
            if label_hit:
                label_hits[mode_name] += 1
        rows.append(row)

    scores = pd.DataFrame(rows)
    scores.to_csv(outdir / "scores.csv", index=False)
    summary = {
        "num_cases": len(cases),
        "totals": totals,
        "passed": passed,
        "pass_rate": {name: passed[name] / len(cases) for name in results},
        "paper_label_hits": label_hits,
        "paper_label_hit_rate": {name: label_hits[name] / len(cases) for name in results},
    }
    return scores, summary


def write_usage(outdir: Path, usage_by_mode: Dict[str, Dict[str, Any]]) -> None:
    rows = []
    for mode, usage in usage_by_mode.items():
        row = {"mode": mode, **usage}
        row.pop("by_model", None)
        rows.append(row)
    pd.DataFrame(rows).to_csv(outdir / "usage.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run macrophage subtype benchmark")
    parser.add_argument("--model", default=MODEL)
    parser.add_argument("--provider", default=PROVIDER)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument(
        "--cases-csv",
        default=None,
        help="Optional benchmark case CSV with case_id, markers, expected_terms, and min_score columns.",
    )
    args = parser.parse_args()

    if not os.environ.get("OPENROUTER_API_KEY"):
        print("SKIP: OPENROUTER_API_KEY is not set.")
        return 0

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.output_dir) if args.output_dir else RESULTS_DIR / timestamp
    outdir.mkdir(parents=True, exist_ok=True)
    cases = load_cases(args.cases_csv)
    marker_dataframe(cases).to_csv(outdir / "inputs.csv", index=False)
    pd.DataFrame(cases).to_csv(outdir / "cases.csv", index=False)

    print(f"Writing benchmark outputs to {outdir}")

    results = {}
    usage_by_mode = {}

    results["cassia_baseline"], usage_by_mode["cassia_baseline"] = run_cassia_mode(
        outdir=outdir,
        mode_name="cassia_baseline",
        use_reference=False,
        model=args.model,
        provider=args.provider,
        cases=cases,
    )
    results["cassia_reference"], usage_by_mode["cassia_reference"] = run_cassia_mode(
        outdir=outdir,
        mode_name="cassia_reference",
        use_reference=True,
        model=args.model,
        provider=args.provider,
        cases=cases,
    )
    results["direct_kimi"], usage_by_mode["direct_kimi"] = run_direct_kimi(outdir, args.model, cases)

    scores, summary = evaluate_modes(results, outdir, cases)
    write_usage(outdir, usage_by_mode)

    summary.update({
        "model": args.model,
        "provider": args.provider,
        "cases_csv": args.cases_csv,
        "output_dir": str(outdir),
        "usage": usage_by_mode,
    })
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(scores[[
        "case_id",
        "cassia_baseline_score",
        "cassia_reference_score",
        "direct_kimi_score",
        "cassia_baseline_paper_label_hit",
        "cassia_reference_paper_label_hit",
        "direct_kimi_paper_label_hit",
    ]].to_string(index=False))
    print(json.dumps({
        "totals": summary["totals"],
        "passed": summary["passed"],
        "pass_rate": summary["pass_rate"],
        "paper_label_hits": summary["paper_label_hits"],
        "paper_label_hit_rate": summary["paper_label_hit_rate"],
        "usage": {
            mode: {
                "requests": usage.get("requests"),
                "total_tokens": usage.get("total_tokens"),
                "cost": usage.get("cost"),
            }
            for mode, usage in usage_by_mode.items()
        },
    }, indent=2))

    if summary["passed"]["cassia_reference"] < summary["passed"]["cassia_baseline"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
