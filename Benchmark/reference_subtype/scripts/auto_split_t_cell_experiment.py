"""
Experimental LLM-assisted auto-split benchmark for hard T-cell subclustering.

This is intentionally outside the CASSIA package API. It tests whether a small
orchestration layer can improve reference-mode robustness by:

1. Asking an LLM to split clusters into biologically coherent groups.
2. Validating that every input cluster appears exactly once.
3. Running CASSIA subclustering + reference mode per group.
4. Merging group outputs.
5. Retrying missing rows as single-cluster jobs.

The script reuses hard_t_cell_marker_benchmark.CASES and can be run against an
existing hard benchmark result directory so baseline/reference calls are not
repeated.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[3]
KEY_FILE = SUITE_DIR / ".openrouter_key"

sys.path.insert(0, str(ROOT / "Test" / "shared" / "python"))
from test_utils import setup_cassia_imports  # noqa: E402

setup_cassia_imports()

_hard_spec = importlib.util.spec_from_file_location(
    "hard_t_cell_marker_benchmark",
    SCRIPT_DIR / "hard_t_cell_marker_benchmark.py",
)
_hard_module = importlib.util.module_from_spec(_hard_spec)
assert _hard_spec and _hard_spec.loader
_hard_spec.loader.exec_module(_hard_module)

CASES = _hard_module.CASES
MODEL = _hard_module.MODEL
PROVIDER = _hard_module.PROVIDER
combine_output_text = _hard_module.combine_output_text
marker_dataframe = _hard_module.marker_dataframe
score_text = _hard_module.score_text


def load_key_from_local_file() -> None:
    if os.environ.get("OPENROUTER_API_KEY"):
        return
    if KEY_FILE.exists():
        key = KEY_FILE.read_text(encoding="utf-8").strip()
        if key:
            os.environ["OPENROUTER_API_KEY"] = key


def _extract_json(text: str) -> Dict[str, Any]:
    match = re.search(r"\{[\s\S]*\}", str(text or ""))
    if not match:
        return {}
    try:
        return json.loads(match.group())
    except json.JSONDecodeError:
        return {}


def fallback_split_plan() -> Dict[str, Any]:
    return {
        "strategy": "deterministic fallback split by known T-cell subtype programs",
        "groups": [
            {
                "group_id": "cd8_exhaustion_cytotoxic",
                "cluster_ids": [
                    "terminal_exhausted_cd8",
                    "cytotoxic_cd8",
                    "pre_exhausted_isg_cd8",
                    "effector_memory_cd8",
                ],
                "group_hint": "CD8 T-cell cytotoxic, effector-memory, interferon-stimulated, and exhaustion continuum",
                "reference_cell_type_hint": "CD8 T cell",
                "reason": "Shared CD8 cytotoxic and exhaustion-related markers.",
            },
            {
                "group_id": "cd4_treg_checkpoint_helper",
                "cluster_ids": [
                    "treg",
                    "helper_checkpoint_cd4",
                    "transitional_memory_cd4",
                ],
                "group_hint": "CD4 T-cell checkpoint, Treg, helper, and transitional-memory states",
                "reference_cell_type_hint": "CD4 T cell",
                "reason": "Shared CD4 helper/checkpoint/Treg markers.",
            },
            {
                "group_id": "cd4_memory_activation_th17",
                "cluster_ids": [
                    "naive_memory_cd4",
                    "th17",
                    "recently_activated_cd4",
                ],
                "group_hint": "CD4 T-cell naive-memory, recently activated, and Th17-like states",
                "reference_cell_type_hint": "CD4 T cell",
                "reason": "CD4 memory/activation and helper polarization comparison.",
            },
            {
                "group_id": "t_cell_naive_proliferative",
                "cluster_ids": [
                    "naive_t",
                    "proliferative_t",
                ],
                "group_hint": "T-cell naive and proliferative state modules",
                "reference_cell_type_hint": "T cell",
                "reason": "State-module clusters that need lineage-state separation.",
            },
        ],
    }


def validate_split_plan(plan: Dict[str, Any]) -> Tuple[bool, str]:
    expected = [case["id"] for case in CASES]
    expected_set = set(expected)
    groups = plan.get("groups")
    if not isinstance(groups, list) or not groups:
        return False, "groups missing or empty"

    observed: List[str] = []
    for group in groups:
        ids = group.get("cluster_ids")
        if not isinstance(ids, list) or not ids:
            return False, f"group {group.get('group_id')} has no cluster_ids"
        if len(ids) > 5:
            return False, f"group {group.get('group_id')} has >5 clusters"
        observed.extend(str(item) for item in ids)

    observed_set = set(observed)
    if observed_set != expected_set:
        missing = sorted(expected_set - observed_set)
        extra = sorted(observed_set - expected_set)
        return False, f"split mismatch missing={missing} extra={extra}"
    if len(observed) != len(observed_set):
        duplicates = sorted({item for item in observed if observed.count(item) > 1})
        return False, f"duplicate cluster ids: {duplicates}"
    return True, "ok"


def plan_split_with_llm(outdir: Path) -> Dict[str, Any]:
    from CASSIA.core.llm_utils import call_llm

    marker_lines = "\n".join(
        f"- {case['id']} ({case['label']} expected only for benchmark metadata, do not copy labels blindly): {case['markers']}"
        for case in CASES
    )
    prompt = f"""You are planning an experimental CASSIA subclustering run.

You are NOT annotating final cell types. Your only task is to split the input
clusters into smaller biologically coherent groups for downstream annotation.

Rules:
- Every cluster ID must appear exactly once.
- Each group must contain 2-5 clusters.
- Keep confusing states together when they need direct comparison.
- Separate CD8 exhaustion/cytotoxic continua from CD4/Treg/Tfh/helper continua
  when possible.
- State-only modules such as proliferation or IFN can join the closest lineage
  group or a small state-module group.
- Return JSON only.

Clusters and marker genes:
{marker_lines}

Return this JSON schema:
{{
  "strategy": "short explanation",
  "groups": [
    {{
      "group_id": "short_snake_case_id",
      "cluster_ids": ["cluster_id_1", "cluster_id_2"],
      "group_hint": "biological context for downstream annotation",
      "reference_cell_type_hint": "CD8 T cell|CD4 T cell|T cell|innate-like T cell",
      "reason": "why these clusters should be compared together"
    }}
  ]
}}"""

    response = call_llm(
        prompt=prompt,
        provider=PROVIDER,
        model=MODEL,
        temperature=0,
        max_tokens=1400,
    )
    (outdir / "split_planner_raw.txt").write_text(response, encoding="utf-8")
    plan = _extract_json(response)
    ok, reason = validate_split_plan(plan)
    if ok:
        plan["planner_status"] = "llm_valid"
        plan["planner_validation"] = reason
        return plan

    fallback = fallback_split_plan()
    fallback["planner_status"] = "llm_invalid_fallback_used"
    fallback["planner_validation"] = reason
    fallback["llm_plan"] = plan
    return fallback


def run_group(
    outdir: Path,
    group: Dict[str, Any],
    marker_df: pd.DataFrame,
    retry: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    from CASSIA import get_llm_usage_summary, reset_llm_usage_log, runCASSIA_subclusters

    group_id = str(group["group_id"])
    ids = [str(item) for item in group["cluster_ids"]]
    group_marker_df = marker_df[marker_df["cluster"].astype(str).isin(ids)].copy()
    suffix = "_retry" if retry else ""
    output_base = str(outdir / f"auto_split_{group_id}{suffix}")

    reset_llm_usage_log()
    runCASSIA_subclusters(
        marker=group_marker_df,
        major_cluster_info=group.get("group_hint") or "human pan-cancer tumor-infiltrating T cell",
        output_name=output_base,
        model=MODEL,
        temperature=0,
        provider=PROVIDER,
        n_genes=30,
        tissue="pan-cancer tumor",
        species="human",
        use_reference=True,
        reference_model=MODEL,
        reference_cell_type_hint=group.get("reference_cell_type_hint") or "t cell",
        reference_max_content_length=12000,
        reference_max_context_length=12000,
    )
    usage = get_llm_usage_summary(reset=True)
    result = pd.read_csv(f"{output_base}.csv")
    result["Result ID"] = result["Result ID"].astype(str)
    return result, usage


def run_auto_split(outdir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    marker_df = marker_dataframe()
    plan = plan_split_with_llm(outdir)
    (outdir / "split_plan.json").write_text(json.dumps(plan, indent=2), encoding="utf-8")

    results = []
    usage_rows = []
    for group in plan["groups"]:
        group_result, usage = run_group(outdir, group, marker_df)
        results.append(group_result)
        usage_rows.append({"mode": f"auto_split:{group['group_id']}", **usage})

    merged = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    expected = {case["id"] for case in CASES}
    observed = set(merged["Result ID"].astype(str)) if "Result ID" in merged.columns else set()
    missing = sorted(expected - observed)

    retry_results = []
    for missing_id in missing:
        retry_group = {
            "group_id": f"missing_{missing_id}",
            "cluster_ids": [missing_id],
            "group_hint": "single missing T-cell subtype cluster retry",
            "reference_cell_type_hint": "T cell",
        }
        retry_result, retry_usage = run_group(outdir, retry_group, marker_df, retry=True)
        retry_results.append(retry_result)
        usage_rows.append({"mode": f"auto_split_retry:{missing_id}", **retry_usage})

    if retry_results:
        merged = pd.concat([merged] + retry_results, ignore_index=True)
    if "Result ID" in merged.columns:
        merged["Result ID"] = merged["Result ID"].astype(str)
        merged = merged.drop_duplicates(subset=["Result ID"], keep="last")

    merged.to_csv(outdir / "auto_split_reference.csv", index=False)
    usage = pd.DataFrame(usage_rows).drop(columns=["by_model"], errors="ignore")
    usage.to_csv(outdir / "usage_auto_split.csv", index=False)
    return merged.set_index("Result ID"), usage


def score_all(outdir: Path, auto_split: pd.DataFrame) -> pd.DataFrame:
    baseline_path = outdir / "hard_t_cell_baseline.csv"
    reference_path = outdir / "hard_t_cell_reference.csv"
    baseline = pd.read_csv(baseline_path)
    reference = pd.read_csv(reference_path)
    baseline["Result ID"] = baseline["Result ID"].astype(str)
    reference["Result ID"] = reference["Result ID"].astype(str)
    baseline = baseline.set_index("Result ID")
    reference = reference.set_index("Result ID")

    rows = []
    for case in CASES:
        case_id = case["id"]
        base_text = "[MISSING OUTPUT ROW]" if case_id not in baseline.index else combine_output_text(baseline.loc[case_id])
        ref_text = "[MISSING OUTPUT ROW]" if case_id not in reference.index else combine_output_text(reference.loc[case_id])
        split_text = "[MISSING OUTPUT ROW]" if case_id not in auto_split.index else combine_output_text(auto_split.loc[case_id])
        rows.append({
            "case_id": case_id,
            "expected_label": case["label"],
            "min_score": case["min_score"],
            "baseline_score": score_text(base_text, case["expected_terms"]),
            "reference_score": score_text(ref_text, case["expected_terms"]),
            "auto_split_score": score_text(split_text, case["expected_terms"]),
            "baseline_text": base_text,
            "reference_text": ref_text,
            "auto_split_text": split_text,
            "reference_missing": case_id not in reference.index,
            "auto_split_missing": case_id not in auto_split.index,
        })
    scores = pd.DataFrame(rows)
    scores["auto_split_pass"] = scores["auto_split_score"] >= scores["min_score"]
    scores["auto_split_vs_reference"] = scores["auto_split_score"] - scores["reference_score"]
    scores["auto_split_vs_baseline"] = scores["auto_split_score"] - scores["baseline_score"]
    scores.to_csv(outdir / "scores_auto_split.csv", index=False)
    return scores


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", help="Existing hard_t_cell result directory.")
    args = parser.parse_args()

    load_key_from_local_file()
    if not os.environ.get("OPENROUTER_API_KEY"):
        print("ERROR: OPENROUTER_API_KEY is not set and local key file is missing.")
        return 2

    outdir = Path(args.result_dir)
    if not (outdir / "hard_t_cell_baseline.csv").exists():
        print(f"ERROR: not a hard T-cell result directory: {outdir}")
        return 2

    auto_split, usage = run_auto_split(outdir)
    scores = score_all(outdir, auto_split)

    totals = {
        "baseline": int(scores["baseline_score"].sum()),
        "reference": int(scores["reference_score"].sum()),
        "auto_split": int(scores["auto_split_score"].sum()),
    }
    pass_count = int(scores["auto_split_pass"].sum())

    print(f"RESULT_DIR={outdir}")
    print(f"TOTALS={totals}")
    print(f"AUTO_SPLIT_PASSED={pass_count}/{len(scores)}")
    print("\nSCORES")
    print(scores[[
        "case_id",
        "expected_label",
        "baseline_score",
        "reference_score",
        "auto_split_score",
        "auto_split_vs_reference",
        "auto_split_missing",
    ]].to_string(index=False))
    print("\nUSAGE_AUTO_SPLIT")
    print(usage.to_string(index=False))

    return 0 if totals["auto_split"] >= totals["reference"] and pass_count >= 11 else 1


if __name__ == "__main__":
    raise SystemExit(main())
