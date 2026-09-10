"""
Run the core CASSIA benchmark against runCASSIA_batch.

The benchmark manifest contains multiple tissues, so this runner groups cases by
dataset/tissue/species and calls runCASSIA_batch once per group. The resulting
summary CSVs are merged back into one run directory.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = SUITE_DIR / "results"
DEFAULT_CASES = SUITE_DIR / "cases" / "core_50.csv"
DEFAULT_MODEL = "openai/gpt-5.6-terra"
DEFAULT_PROVIDER = "openrouter"

sys.path.insert(0, str(ROOT / "CASSIA_python"))

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_core_batch_results import score_run_dir  # noqa: E402


def slugify(value: object) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", str(value or "").lower()).strip("_")
    return slug or "group"


def provider_env_var(provider: str) -> Optional[str]:
    provider = provider.lower()
    if provider == "openrouter":
        return "OPENROUTER_API_KEY"
    if provider == "openai":
        return "OPENAI_API_KEY"
    if provider == "anthropic":
        return "ANTHROPIC_API_KEY"
    if provider.startswith("http"):
        return "CUSTOMIZED_API_KEY"
    return None


def load_local_key_if_present(provider: str) -> None:
    env_var = provider_env_var(provider)
    if not env_var or os.environ.get(env_var):
        return

    key_file = SUITE_DIR / f".{provider.lower()}_key"
    if key_file.exists():
        key = key_file.read_text(encoding="utf-8").strip()
        if key:
            os.environ[env_var] = key


def require_api_key(provider: str, allow_free_api: bool) -> None:
    load_local_key_if_present(provider)
    env_var = provider_env_var(provider)
    if allow_free_api or provider.lower().startswith("http"):
        return
    if env_var and os.environ.get(env_var):
        return
    raise RuntimeError(
        f"{env_var or 'provider API key'} is not set. Set it in the environment "
        f"or create {SUITE_DIR / ('.' + provider.lower() + '_key')}. "
        "Use --allow-free-api only for smoke tests, not benchmark runs."
    )


def load_cases(path: Path, difficulties: Optional[Iterable[str]], limit: Optional[int]) -> pd.DataFrame:
    cases = pd.read_csv(path)
    required = {
        "case_id",
        "dataset",
        "tissue",
        "species",
        "expected_cell_type",
        "difficulty",
        "marker_list",
    }
    missing = sorted(required - set(cases.columns))
    if missing:
        raise ValueError(f"Cases file is missing columns: {', '.join(missing)}")

    cases["case_id"] = cases["case_id"].astype(str)
    if difficulties:
        allowed = {value.lower() for value in difficulties}
        cases = cases[cases["difficulty"].astype(str).str.lower().isin(allowed)]
    if limit is not None:
        cases = cases.head(limit)
    if cases.empty:
        raise ValueError("No benchmark cases selected")
    return cases.reset_index(drop=True)


def marker_dataframe(group: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "cluster": group["case_id"].astype(str),
        "markers": group["marker_list"].astype(str),
    })


def group_id_for(group_key: Dict[str, str]) -> str:
    return "_".join(slugify(group_key[key]) for key in ["dataset", "tissue", "species"])


def flatten_usage(usage: Dict[str, object]) -> Dict[str, object]:
    flat = dict(usage or {})
    if "by_model" in flat:
        flat["by_model"] = json.dumps(flat["by_model"], sort_keys=True)
    return flat


def write_run_config(outdir: Path, args: argparse.Namespace, cases: pd.DataFrame) -> None:
    config = {
        "model": args.model,
        "provider": args.provider,
        "temperature": args.temperature,
        "max_workers": args.max_workers,
        "validator_involvement": args.validator_involvement,
        "reasoning": args.reasoning,
        "use_reference": args.use_reference,
        "reference_model": args.reference_model,
        "cases_csv": str(args.cases_csv),
        "num_cases": int(len(cases)),
        "difficulties": args.difficulty,
    }
    (outdir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")


def run_group(
    outdir: Path,
    group_id: str,
    group: pd.DataFrame,
    args: argparse.Namespace,
) -> Dict[str, object]:
    from CASSIA import get_llm_usage_summary, reset_llm_usage_log, runCASSIA_batch

    output_base = outdir / group_id / "cassia_batch"
    output_base.parent.mkdir(parents=True, exist_ok=True)
    tissue = str(group["tissue"].iloc[0])
    species = str(group["species"].iloc[0])

    reset_llm_usage_log()
    kwargs = {
        "marker": marker_dataframe(group),
        "output_name": str(output_base),
        "n_genes": args.n_genes,
        "model": args.model,
        "temperature": args.temperature,
        "tissue": tissue,
        "species": species,
        "additional_info": args.additional_info,
        "celltype_column": "cluster",
        "gene_column_name": "markers",
        "max_workers": args.max_workers,
        "provider": args.provider,
        "max_retries": args.max_retries,
        "validator_involvement": args.validator_involvement,
        "reasoning": args.reasoning,
        "use_reference": args.use_reference,
        "reference_model": args.reference_model,
        "reference_cell_type_hint": args.reference_cell_type_hint,
        "validate_api_key_before_start": not args.no_validate_api_key,
        "auto_convert_ids": not args.no_auto_convert_ids,
    }
    runCASSIA_batch(**kwargs)

    summary_path = Path(f"{output_base}_summary.csv")
    conversations_path = Path(f"{output_base}_conversations.json")
    report_path = Path(f"{output_base}_report.html")
    if not summary_path.exists():
        raise RuntimeError(f"Expected summary CSV was not created: {summary_path}")

    usage = flatten_usage(get_llm_usage_summary(reset=True))
    return {
        "group_id": group_id,
        "summary_path": str(summary_path),
        "conversations_path": str(conversations_path),
        "report_path": str(report_path),
        "usage": usage,
    }


def combine_group_outputs(group_results: List[Dict[str, object]], cases: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    case_meta = cases.drop(columns=["marker_list"]).copy()
    frames = []
    for result in group_results:
        group_id = str(result["group_id"])
        summary_path = Path(str(result["summary_path"]))
        frame = pd.read_csv(summary_path)
        frame.insert(0, "run_group", group_id)
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.merge(
        case_meta,
        left_on="Cluster ID",
        right_on="case_id",
        how="left",
        suffixes=("", "_case"),
    )
    combined.to_csv(outdir / "combined_predictions.csv", index=False)
    return combined


def print_dry_run(cases: pd.DataFrame) -> None:
    print(f"Selected {len(cases)} cases")
    for (dataset, tissue, species), group in cases.groupby(["dataset", "tissue", "species"], sort=True):
        print(f"- {dataset} / {tissue} / {species}: {len(group)} cases")
        for row in group[["case_id", "expected_cell_type", "difficulty"]].itertuples(index=False):
            print(f"  - {row.case_id}: {row.expected_cell_type} ({row.difficulty})")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the core runCASSIA_batch benchmark")
    parser.add_argument("--cases-csv", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--provider", default=DEFAULT_PROVIDER)
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--max-workers", type=int, default=5)
    parser.add_argument("--max-retries", type=int, default=1)
    parser.add_argument("--n-genes", type=int, default=50)
    parser.add_argument("--validator-involvement", default="v1")
    parser.add_argument("--reasoning", default=None)
    parser.add_argument("--additional-info", default=None)
    parser.add_argument("--use-reference", action="store_true")
    parser.add_argument("--reference-model", default=None)
    parser.add_argument("--reference-cell-type-hint", default=None)
    parser.add_argument("--no-validate-api-key", action="store_true")
    parser.add_argument("--no-auto-convert-ids", action="store_true")
    parser.add_argument("--allow-free-api", action="store_true")
    parser.add_argument("--difficulty", action="append", choices=["easy", "medium", "hard"], default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cases = load_cases(args.cases_csv, args.difficulty, args.limit)
    if args.dry_run:
        print_dry_run(cases)
        return 0

    require_api_key(args.provider, args.allow_free_api)

    run_id = args.run_id or time.strftime("%Y%m%d_%H%M%S")
    outdir = args.output_dir or RESULTS_DIR / run_id
    outdir.mkdir(parents=True, exist_ok=True)

    cases.to_csv(outdir / "cases.csv", index=False)
    cases[["case_id", "marker_list"]].rename(
        columns={"case_id": "cluster", "marker_list": "markers"}
    ).to_csv(outdir / "inputs.csv", index=False)
    write_run_config(outdir, args, cases)

    print(f"Writing benchmark outputs to {outdir}")
    group_results: List[Dict[str, object]] = []
    usage_rows = []

    for (dataset, tissue, species), group in cases.groupby(["dataset", "tissue", "species"], sort=True):
        group_key = {"dataset": dataset, "tissue": tissue, "species": species}
        group_id = group_id_for(group_key)
        print(f"\nRunning {group_id}: {len(group)} cases")
        result = run_group(outdir, group_id, group.reset_index(drop=True), args)
        group_results.append(result)
        usage_rows.append({
            "group_id": group_id,
            "dataset": dataset,
            "tissue": tissue,
            "species": species,
            "cases": int(len(group)),
            **result["usage"],
        })

    combine_group_outputs(group_results, cases, outdir)
    pd.DataFrame(usage_rows).to_csv(outdir / "usage.csv", index=False)
    (outdir / "group_outputs.json").write_text(
        json.dumps(group_results, indent=2),
        encoding="utf-8",
    )

    score_run_dir(outdir)
    print(f"\nCompleted run: {outdir}")
    print(f"Scores: {outdir / 'scores.csv'}")
    print(f"Summary: {outdir / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
