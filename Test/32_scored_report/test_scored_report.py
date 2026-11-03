"""Regression coverage for issue #25: scored reports retain annotation history."""

import json
import sys
from pathlib import Path

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "CASSIA_python"))

from CASSIA.evaluation import scoring


def test_scored_report_restores_companion_json(tmp_path, monkeypatch):
    summary = tmp_path / "cassia_results_summary.csv"
    scored = tmp_path / "cassia_annotation_scored.csv"
    pd.DataFrame([
        {
            "Cluster ID": 0,
            "Predicted General Cell Type": "T cell",
            "Marker List": "CD3D, CD3E",
            "Tissue": "Large Intestine",
            "Species": "Human",
        },
        {
            "Cluster ID": "cluster 2",
            "Predicted General Cell Type": "B cell",
            "Marker List": "MS4A1, CD79A",
            "Tissue": "Large Intestine",
            "Species": "Human",
        },
    ]).to_csv(summary, index=False)
    conversations = {
        "0": {
            "annotations": ["T-cell annotation details"],
            "validations": ["VALIDATION PASSED: T-cell markers agree"],
            "formatting": '{"main_cell_type": "T cell"}',
        },
        "cluster 2": {
            "annotations": ["B-cell annotation details"],
            "validations": ["VALIDATION PASSED: B-cell markers agree"],
            "formatting": '{"main_cell_type": "B cell"}',
        },
    }
    (tmp_path / "cassia_results_conversations.json").write_text(
        json.dumps(conversations), encoding="utf-8"
    )

    seen_histories = []

    def fake_score(_origin, _markers, history, **_kwargs):
        seen_histories.append(history)
        return 91, "Scoring explanation for this cluster"

    monkeypatch.setattr(scoring, "score_single_analysis", fake_score)
    scoring.runCASSIA_score_batch(
        str(summary), str(scored), max_workers=1, model="test-model"
    )

    report = scored.with_name("cassia_annotation_scored_report.html").read_text(
        encoding="utf-8"
    )
    assert len(seen_histories) == 2
    assert any("T-cell annotation details" in item for item in seen_histories)
    assert any("B-cell annotation details" in item for item in seen_histories)
    assert "T-cell annotation details" in report
    assert "B-cell annotation details" in report
    assert "T-cell markers agree" in report
    assert "B-cell markers agree" in report
    assert "Scoring explanation for this cluster" in report
    assert "No annotation data available." not in report


def test_already_scored_input_still_generates_report(tmp_path):
    scored_input = tmp_path / "results_summary.csv"
    output = tmp_path / "scored.csv"
    pd.DataFrame([{
        "Cluster ID": "cluster 1",
        "Predicted General Cell Type": "T cell",
        "Marker List": "CD3D",
        "Tissue": "Large Intestine",
        "Species": "Human",
        "Score": 92,
        "Scoring_Reasoning": "Existing scoring explanation",
    }]).to_csv(scored_input, index=False)
    (tmp_path / "results_conversations.json").write_text(
        json.dumps({"cluster 1": {"annotations": ["Existing annotation detail"]}}),
        encoding="utf-8",
    )

    scoring.runCASSIA_score_batch(str(scored_input), str(output), model="test-model")

    report = output.with_name("scored_report.html").read_text(encoding="utf-8")
    assert "Existing annotation detail" in report
    assert "Existing scoring explanation" in report
