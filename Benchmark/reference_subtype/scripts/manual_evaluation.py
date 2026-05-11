"""
Fixed manual-evaluation rubric for reference subtype benchmarks.

The benchmark's automatic expected-term score is useful triage, but the final
quality call should be a structured human review. This module makes that review
reproducible by validating one row per case and mode, computing totals, and
writing summary artifacts.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


RUBRIC: List[Tuple[str, int, str]] = [
    (
        "program_accuracy",
        4,
        "Biological subtype/program match: 4 exact, 3 correct broad program with minor miss, "
        "2 partial/hybrid undercalled, 1 weak, 0 wrong.",
    ),
    (
        "subtype_specificity",
        2,
        "Subtype resolution: 2 subtype-level label, 1 broad macrophage/TAM state, 0 generic/wrong.",
    ),
    (
        "marker_evidence",
        2,
        "Marker reasoning: 2 uses defining markers and conflicts correctly, 1 partial, 0 unsupported.",
    ),
    (
        "ambiguity_handling",
        1,
        "Handles hybrid/alternative states and avoids overclaiming when needed.",
    ),
    (
        "traceability",
        1,
        "Output is traceable to a consensus program or paper alias. Exact paper alias is helpful but not required.",
    ),
]

MAX_SCORE = sum(max_score for _name, max_score, _description in RUBRIC)
PASS_MIN_TOTAL = 7
PASS_MIN_PROGRAM_ACCURACY = 3
EXCELLENT_MIN_TOTAL = 9


def metric_names() -> List[str]:
    return [name for name, _max_score, _description in RUBRIC]


def rubric_dict() -> Dict[str, Dict[str, object]]:
    return {
        name: {"max_score": max_score, "description": description}
        for name, max_score, description in RUBRIC
    }


def load_manual_scores(scores_csv: str | Path) -> pd.DataFrame:
    path = Path(scores_csv)
    df = pd.read_csv(path)
    required = {"case_id", "mode", "notes", *metric_names()}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Manual score file is missing columns: {', '.join(missing)}")

    for metric, max_score, _description in RUBRIC:
        df[metric] = pd.to_numeric(df[metric], errors="raise")
        invalid = df[(df[metric] < 0) | (df[metric] > max_score)]
        if not invalid.empty:
            bad_cases = ", ".join(invalid["case_id"].astype(str).head(5))
            raise ValueError(
                f"Metric {metric} must be between 0 and {max_score}; invalid cases: {bad_cases}"
            )

    duplicated = df.duplicated(["case_id", "mode"], keep=False)
    if duplicated.any():
        duplicates = df.loc[duplicated, ["case_id", "mode"]].drop_duplicates()
        duplicate_text = ", ".join(
            f"{row.case_id}/{row.mode}" for row in duplicates.itertuples()
        )
        raise ValueError(f"Duplicate manual score rows: {duplicate_text}")

    return df


def summarize_manual_scores(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    scored = df.copy()
    scored["manual_total"] = scored[metric_names()].sum(axis=1)
    scored["manual_pass"] = (
        (scored["manual_total"] >= PASS_MIN_TOTAL)
        & (scored["program_accuracy"] >= PASS_MIN_PROGRAM_ACCURACY)
    )
    scored["manual_excellent"] = scored["manual_total"] >= EXCELLENT_MIN_TOTAL

    rows = []
    for mode, group in scored.groupby("mode", sort=True):
        rows.append({
            "mode": mode,
            "cases": int(len(group)),
            "total_score": float(group["manual_total"].sum()),
            "max_score": int(len(group) * MAX_SCORE),
            "mean_score": float(group["manual_total"].mean()),
            "pass_count": int(group["manual_pass"].sum()),
            "pass_rate": float(group["manual_pass"].mean()),
            "excellent_count": int(group["manual_excellent"].sum()),
            "excellent_rate": float(group["manual_excellent"].mean()),
            **{
                f"mean_{metric}": float(group[metric].mean())
                for metric in metric_names()
            },
        })
    summary_df = pd.DataFrame(rows)

    summary = {
        "rubric": rubric_dict(),
        "max_score_per_case": MAX_SCORE,
        "pass_rule": {
            "manual_total_min": PASS_MIN_TOTAL,
            "program_accuracy_min": PASS_MIN_PROGRAM_ACCURACY,
        },
        "excellent_rule": {"manual_total_min": EXCELLENT_MIN_TOTAL},
        "by_mode": summary_df.to_dict(orient="records"),
    }
    return scored, summary_df, summary


def write_manual_summary(scores_csv: str | Path, output_dir: str | Path | None = None) -> Dict[str, Path]:
    input_path = Path(scores_csv)
    outdir = Path(output_dir) if output_dir else input_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_manual_scores(input_path)
    scored, summary_df, summary = summarize_manual_scores(df)

    scored_path = outdir / "manual_scores_scored.csv"
    summary_csv_path = outdir / "manual_score_summary.csv"
    summary_json_path = outdir / "manual_score_summary.json"

    scored.to_csv(scored_path, index=False)
    summary_df.to_csv(summary_csv_path, index=False)
    summary_json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return {
        "scored": scored_path,
        "summary_csv": summary_csv_path,
        "summary_json": summary_json_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize manual subtype benchmark scores")
    parser.add_argument("scores_csv", help="Manual score CSV with one row per case and mode")
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    outputs = write_manual_summary(args.scores_csv, args.output_dir)
    for label, path in outputs.items():
        print(f"{label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
