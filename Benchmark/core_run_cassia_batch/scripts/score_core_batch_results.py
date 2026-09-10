"""
Deterministic triage scoring for the core runCASSIA_batch benchmark.

This is intentionally lightweight. It checks whether the prediction text
contains expected broad/specific terms, then leaves final biological judgement
to manual review when a run is close or ambiguous.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CASES = SUITE_DIR / "cases" / "core_50.csv"

PREDICTION_TEXT_COLUMNS = [
    "Predicted General Cell Type",
    "Predicted Detailed Cell Type",
    "Possible Mixed Cell Types",
]


def parse_expected_terms(value: object) -> List[List[str]]:
    groups: List[List[str]] = []
    for group in str(value or "").split(";"):
        terms = [term.strip() for term in group.split("/") if term.strip()]
        if terms:
            groups.append(terms)
    return groups


def normalize_for_search(value: object) -> str:
    text = str(value or "").lower()
    text = text.replace("-", " ")
    text = re.sub(r"\s+", " ", text)
    return text


def compact(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def prediction_text(row: pd.Series) -> str:
    return " | ".join(str(row.get(column, "")) for column in PREDICTION_TEXT_COLUMNS)


def group_hit(text: str, term_group: Iterable[str]) -> bool:
    normalized_text = normalize_for_search(text)
    compact_text = compact(text)
    for term in term_group:
        normalized_term = normalize_for_search(term)
        if normalized_term and normalized_term in normalized_text:
            return True
        compact_term = compact(term)
        if compact_term and compact_term in compact_text:
            return True
    return False


def score_text(text: str, expected_terms: List[List[str]]) -> Tuple[int, List[str], List[str]]:
    hit_groups: List[str] = []
    missed_groups: List[str] = []
    for group in expected_terms:
        rendered = "/".join(group)
        if group_hit(text, group):
            hit_groups.append(rendered)
        else:
            missed_groups.append(rendered)
    return len(hit_groups), hit_groups, missed_groups


def load_cases(cases_csv: str | Path) -> pd.DataFrame:
    cases = pd.read_csv(cases_csv)
    required = {
        "case_id",
        "dataset",
        "tissue",
        "species",
        "expected_cell_type",
        "difficulty",
        "case_type",
        "expected_terms",
        "min_score",
    }
    missing = sorted(required - set(cases.columns))
    if missing:
        raise ValueError(f"Cases file is missing columns: {', '.join(missing)}")
    cases["case_id"] = cases["case_id"].astype(str)
    cases["min_score"] = pd.to_numeric(cases["min_score"], errors="raise").astype(int)
    return cases


def load_predictions(predictions_csv: str | Path) -> pd.DataFrame:
    predictions = pd.read_csv(predictions_csv)
    required = {"Cluster ID", *PREDICTION_TEXT_COLUMNS}
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise ValueError(f"Predictions file is missing columns: {', '.join(missing)}")
    predictions["Cluster ID"] = predictions["Cluster ID"].astype(str)
    return predictions


def summarize(scores: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, object]]:
    rows = []
    group_columns = ["dataset", "difficulty", "case_type"]
    for column in group_columns:
        for value, group in scores.groupby(column, sort=True):
            rows.append({
                "group_by": column,
                "group": value,
                "cases": int(len(group)),
                "pass_count": int(group["term_pass"].sum()),
                "pass_rate": float(group["term_pass"].mean()),
                "mean_term_score": float(group["term_score"].mean()),
                "missing_count": int(group["missing_prediction"].sum()),
            })

    summary_df = pd.DataFrame(rows)
    summary = {
        "num_cases": int(len(scores)),
        "pass_count": int(scores["term_pass"].sum()),
        "pass_rate": float(scores["term_pass"].mean()),
        "mean_term_score": float(scores["term_score"].mean()),
        "missing_count": int(scores["missing_prediction"].sum()),
        "by_group": summary_df.to_dict(orient="records"),
    }
    return summary_df, summary


def score_predictions(cases: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    indexed_predictions = predictions.drop_duplicates("Cluster ID").set_index("Cluster ID")
    rows = []

    for case in cases.to_dict(orient="records"):
        case_id = str(case["case_id"])
        if case_id in indexed_predictions.index:
            pred = indexed_predictions.loc[case_id]
            text = prediction_text(pred)
            missing = False
            output_row = pred.to_dict()
        else:
            text = ""
            missing = True
            output_row = {column: "" for column in predictions.columns}
            output_row["Cluster ID"] = case_id

        expected_terms = parse_expected_terms(case["expected_terms"])
        term_score, hit_groups, missed_groups = score_text(text, expected_terms)
        min_score = int(case["min_score"])

        rows.append({
            **case,
            **output_row,
            "prediction_text": text,
            "term_score": term_score,
            "term_pass": (term_score >= min_score) and not missing,
            "hit_groups": ";".join(hit_groups),
            "missed_groups": ";".join(missed_groups),
            "missing_prediction": missing,
        })

    return pd.DataFrame(rows)


def score_run_dir(
    run_dir: str | Path,
    cases_csv: str | Path | None = None,
    predictions_csv: str | Path | None = None,
) -> Dict[str, Path]:
    run_path = Path(run_dir)
    cases_path = Path(cases_csv) if cases_csv else run_path / "cases.csv"
    predictions_path = Path(predictions_csv) if predictions_csv else run_path / "combined_predictions.csv"

    cases = load_cases(cases_path)
    predictions = load_predictions(predictions_path)
    scores = score_predictions(cases, predictions)
    summary_df, summary = summarize(scores)

    scores_path = run_path / "scores.csv"
    summary_csv_path = run_path / "score_summary.csv"
    summary_json_path = run_path / "summary.json"

    scores.to_csv(scores_path, index=False)
    summary_df.to_csv(summary_csv_path, index=False)
    summary_json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    return {
        "scores": scores_path,
        "summary_csv": summary_csv_path,
        "summary_json": summary_json_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Score a core runCASSIA_batch benchmark run")
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--cases-csv", type=Path, default=None)
    parser.add_argument("--predictions-csv", type=Path, default=None)
    args = parser.parse_args()

    outputs = score_run_dir(args.run_dir, args.cases_csv, args.predictions_csv)
    for label, path in outputs.items():
        print(f"{label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
