"""
CASSIA Test 28: CLI Core
========================
Fast tests for the CLI run-folder and parser logic. These do not call external
agent CLIs or LLM APIs.

Usage:
    python test_cli_core.py
"""

import csv
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest


sys.path.insert(0, str(Path(__file__).parent.parent.parent / "CASSIA_python"))

from CASSIA.cli import main
from CASSIA.cli.backends import AGENT_BACKENDS, AgentCLIBackend, render_agent_argv
from CASSIA.cli.boost import (
    build_boost_followup_prompt,
    build_boost_prompt,
    extract_check_genes,
    extract_candidate_set,
    get_cluster_top_markers,
    parse_gene_args,
    query_marker_genes,
)
from CASSIA.cli.judge import build_judge_prompt, normalize_judgments
from CASSIA.cli.agent_validation import (
    build_revision_prompt,
    build_validator_prompt,
    run_validated_annotation,
)
from CASSIA.cli.runner import (
    MarkerCluster,
    build_annotation_prompt,
    build_annotation_prompt_v1,
    extract_json_object,
    load_marker_clusters,
    run_api_annotation,
)
from CASSIA.cli.result_schema import (
    ANNOTATION_SCHEMA_VERSION,
    normalize_annotation_payload,
    normalize_fused_boost_payload,
)


@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


def _capture_help(argv):
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        try:
            main(argv)
        except SystemExit as exc:
            assert exc.code == 0
    return buffer.getvalue()


def _capture_main(argv):
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        code = main(argv)
    return code, buffer.getvalue()


def test_cli_help_includes_workflow_examples():
    top_help = _capture_help(["--help"])
    validate_help = _capture_help(["validate", "--help"])
    examples_help = _capture_help(["examples", "--help"])
    guide_help = _capture_help(["guide", "--help"])
    annotate_help = _capture_help(["annotate", "--help"])
    consensus_help = _capture_help(["consensus", "--help"])
    judge_help = _capture_help(["judge", "--help"])
    boost_auto_help = _capture_help(["boost", "auto", "--help"])
    boost_run_help = _capture_help(["boost", "run", "--help"])
    subcluster_help = _capture_help(["subcluster", "run", "--help"])

    assert "Start an annotation:" in top_help
    assert "Other common workflows:" in top_help
    assert "cassia validate markers.csv" in top_help
    assert "cassia examples --out cassia_example" in top_help
    assert "Validate marker CSV structure" in validate_help
    assert "Create marker CSVs" in examples_help
    assert "coding agents" in guide_help
    assert "cassia annotate -i markers.csv --backend codex-cli" in top_help
    assert "Examples:" in annotate_help
    assert "cassia consensus --inputs runs/codex/summary.csv" in consensus_help
    assert "cassia boost auto --run runs/brain_codex" in boost_auto_help
    assert "--mode {review,fused}" in boost_run_help
    assert "--fused-prompt-version {v2-compact,v2,v3}" in boost_run_help
    assert "Default: unlimited" in boost_run_help
    assert "--model MODEL" in boost_run_help
    assert "cassia subcluster run --markers cd8_subcluster_markers.csv" in subcluster_help
    assert "--mode {one-shot,validated,fused-boost}" in annotate_help
    assert "Annotation modes:" in top_help
    assert "Agent subscription backends:" in annotate_help
    assert "opencode" in annotate_help
    assert "results.json" in annotate_help
    assert "help" in top_help
    assert "cassia judge" in top_help
    assert "stable-judge-v1.1" in judge_help
    assert "--judge-reasoning-effort" in judge_help


def test_cli_help_command_supports_nested_topics():
    code, empty_help = _capture_main([])
    assert code == 0
    assert "Start an annotation:" in empty_help

    code, top_help = _capture_main(["help"])
    assert code == 0
    assert "Annotation modes:" in top_help

    code, annotate_help = _capture_main(["help", "annotate"])
    assert code == 0
    assert "usage: cassia annotate" in annotate_help
    assert "fused-boost" in annotate_help

    code, nested_help = _capture_main(["help", "boost", "run"])
    assert code == 0
    assert "usage: cassia boost run" in nested_help
    assert "--fused-prompt-version" in nested_help


def test_api_annotation_limit_is_applied_before_batch_dispatch(tmp_dir, monkeypatch):
    marker_path = tmp_dir / "markers.csv"
    marker_path.write_text(
        "cluster,markers\n"
        'first,"CD3D,CD3E,TRBC1"\n'
        'second,"MS4A1,CD79A,CD37"\n',
        encoding="utf-8",
    )
    captured = {}

    def fake_run_cassia_batch(**kwargs):
        captured.update(kwargs)

    import CASSIA.engine.tools_function as tools_function

    monkeypatch.setattr(tools_function, "runCASSIA_batch", fake_run_cassia_batch)
    args = SimpleNamespace(
        input=str(marker_path),
        out=str(tmp_dir / "run"),
        backend="openrouter",
        model="deepseek/deepseek-v4-flash-0731",
        temperature=0,
        tissue="blood",
        species="human",
        additional_info=None,
        celltype_column="cluster",
        gene_column="markers",
        n_genes=50,
        max_workers=1,
        max_retries=1,
        ranking_method="avg_log2FC",
        ascending=None,
        validator_involvement="v0",
        reasoning=None,
        use_reference=False,
        reference_model=None,
        reference_cell_type_hint=None,
        skip_api_key_validation=True,
        quiet=True,
        limit=1,
    )

    assert run_api_annotation(args) == 0
    assert list(captured["marker"]["cluster"]) == ["first"]
    assert captured["celltype_column"] == "cluster"
    assert captured["gene_column_name"] == "markers"


def test_cli_agent_guide_can_print_and_save(tmp_dir):
    code, output = _capture_main(["guide"])
    assert code == 0
    assert "# CASSIA CLI Guide for Coding Agents" in output
    assert "cassia annotate" in output
    assert "cassia.annotation.v1" in output

    out_path = tmp_dir / "CASSIA_AGENT_GUIDE.md"
    code, output = _capture_main(["guide", "--out", str(out_path)])
    assert code == 0
    assert out_path.exists()
    assert "Operating rules" in out_path.read_text(encoding="utf-8")


def test_canonical_result_schema_preserves_workflow_fields():
    standard = normalize_annotation_payload({
        "main_cell_type": "T cell",
        "sub_cell_types": "CD8 T cell, activated T cell",
        "possible_mixed_cell_types": None,
    }, cluster_id="3", markers=["CD3D", "TRAC"])
    assert standard["schema_version"] == ANNOTATION_SCHEMA_VERSION
    assert standard["cluster_id"] == "3"
    assert standard["sub_cell_types"] == ["CD8 T cell", "activated T cell"]

    fused = normalize_fused_boost_payload({
        "final_cell_type": "T cell",
        "final_sub_cell_type": "CD8 T cell",
        "ranked_sub_cell_types": ["activated T cell", "CD8 T cell"],
        "possible_mixed_cell_types": "NK cell",
        "changed_from_original": "false",
        "supporting_markers": "CD3D, TRAC",
    })
    assert fused["schema_version"] == ANNOTATION_SCHEMA_VERSION
    assert fused["annotation_mode"] == "fused_boost"
    assert fused["main_cell_type"] == fused["final_cell_type"] == "T cell"
    assert fused["sub_cell_types"] == ["CD8 T cell", "activated T cell"]
    assert fused["possible_mixed_cell_types"] == ["NK cell"]
    assert fused["changed_from_original"] is False


def test_cli_examples_generates_runnable_project(tmp_dir):
    out_dir = tmp_dir / "cassia_example"
    code, output = _capture_main([
        "examples",
        "--out",
        str(out_dir),
        "--backend",
        "shell",
    ])
    assert code == 0
    assert "Created CASSIA example project" in output
    for relative_path in [
        "README.md",
        "markers.csv",
        "raw_markers.csv",
        "subcluster_markers.csv",
        "toy_agent.py",
        "run.sh",
        "run_offline.sh",
        "consensus_inputs/codex/summary.csv",
        "consensus_inputs/claude/summary.csv",
    ]:
        assert (out_dir / relative_path).exists()

    code, validate_output = _capture_main(["validate", str(out_dir / "markers.csv")])
    assert code == 0
    assert "Status: OK" in validate_output

    consensus_out = out_dir / "runs" / "consensus_test.csv"
    code = main([
        "consensus",
        "--inputs",
        str(out_dir / "consensus_inputs" / "codex"),
        str(out_dir / "consensus_inputs" / "claude"),
        "--out",
        str(consensus_out),
        "--no-html",
    ])
    assert code == 0
    assert consensus_out.exists()

    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).parent.parent.parent / "CASSIA_python")
    env["CASSIA_CMD"] = f"{sys.executable} -m CASSIA.cli"
    env["PYTHON_BIN"] = sys.executable
    completed = subprocess.run(
        ["bash", "run_offline.sh"],
        cwd=str(out_dir),
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout
    assert (out_dir / "runs" / "offline_annotation" / "summary.csv").exists()
    assert (out_dir / "runs" / "offline_subcluster" / "subcluster_results.csv").exists()
    assert (out_dir / "runs" / "example_consensus.csv").exists()


def test_cli_validate_preformatted_marker_list(tmp_dir):
    marker_path = tmp_dir / "markers.csv"
    marker_path.write_text(
        "cluster,markers\n"
        '0,"CD3D, CD3E, TRAC, IL7R, LTB"\n'
        '1,"MS4A1, CD79A, CD74, CD79B, BANK1"\n',
        encoding="utf-8",
    )
    code, output = _capture_main([
        "validate",
        str(marker_path),
        "--backend",
        "codex-cli",
        "--tissue",
        "blood",
    ])
    assert code == 0
    assert "Status: OK" in output
    assert "Detected format: preformatted_marker_list" in output
    assert "Clusters: 2" in output
    assert "cassia annotate" in output


def test_cli_validate_long_marker_table_json(tmp_dir):
    marker_path = tmp_dir / "raw_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "0,CD3D,2.5,0.91,0.04,1e-20\n"
        "0,TRAC,2.1,0.88,0.05,1e-18\n"
        "1,MS4A1,3.2,0.88,0.03,1e-30\n",
        encoding="utf-8",
    )
    code, output = _capture_main([
        "validate",
        str(marker_path),
        "--json",
        "--n-genes",
        "2",
    ])
    payload = json.loads(output)
    assert code == 0
    assert payload["status"] == "warning"
    assert payload["format"] == "seurat_long"
    assert payload["detected_columns"]["cluster_column"] == "cluster"
    assert payload["summary"]["cluster_count"] == 2
    assert payload["warnings"]


def test_cli_validate_missing_ranking_columns_errors(tmp_dir):
    marker_path = tmp_dir / "bad_long_markers.csv"
    marker_path.write_text(
        "cluster,gene,p_val_adj\n"
        "0,CD3D,1e-20\n"
        "1,MS4A1,1e-30\n",
        encoding="utf-8",
    )
    code, output = _capture_main(["validate", str(marker_path), "--celltype-column", "cluster", "--gene-column", "gene"])
    assert code == 1
    assert "Status: ERROR" in output
    assert "Missing: avg_log2FC" in output


def test_extract_json_object():
    payload = extract_json_object(
        "agent log\n```json\n"
        '{"main_cell_type": "T cell", "sub_cell_types": ["CD8 T cell"]}'
        "\n```\n"
    )
    assert payload["main_cell_type"] == "T cell"
    assert payload["sub_cell_types"] == ["CD8 T cell"]


def test_load_marker_clusters_preformatted(tmp_dir):
    marker_path = tmp_dir / "markers.csv"
    marker_path.write_text(
        "cluster,markers\n"
        '0,"CD3D, CD3E, TRAC"\n',
        encoding="utf-8",
    )
    clusters = load_marker_clusters(marker_path, n_genes=2)
    assert len(clusters) == 1
    assert clusters[0].cluster_id == "0"
    assert clusters[0].markers == ["CD3D", "CD3E"]


def test_load_marker_clusters_long_table_with_gene_column_override(tmp_dir):
    marker_path = tmp_dir / "raw_markers_for_loader.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "0,CD3D,2.5,0.91,0.04,1e-20\n"
        "0,TRAC,2.1,0.88,0.05,1e-18\n"
        "1,MS4A1,3.2,0.88,0.03,1e-30\n",
        encoding="utf-8",
    )
    clusters = load_marker_clusters(
        marker_path,
        n_genes=2,
        celltype_column="cluster",
        gene_column="gene",
    )
    by_cluster = {cluster.cluster_id: cluster.markers for cluster in clusters}
    assert by_cluster["0"] == ["CD3D", "TRAC"]
    assert by_cluster["1"] == ["MS4A1"]


def test_cli_dry_run(tmp_dir):
    marker_path = tmp_dir / "markers.csv"
    run_dir = tmp_dir / "run"
    marker_path.write_text(
        "cluster,markers\n"
        '0,"CD3D, CD3E, TRAC"\n',
        encoding="utf-8",
    )
    code = main([
        "annotate",
        "--input",
        str(marker_path),
        "--out",
        str(run_dir),
        "--backend",
        "claude-cli",
        "--dry-run",
        "--limit",
        "1",
    ])
    assert code == 0
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "dry-run"
    assert len(list((run_dir / "prompts").glob("*.md"))) == 1


def test_codex_backend_skips_git_repo_check(tmp_dir):
    argv = render_agent_argv(
        AGENT_BACKENDS["codex-cli"],
        "hello",
        tmp_dir / "prompt.md",
        {"agent_output_file": str(tmp_dir / "out.txt")},
    )
    assert "--skip-git-repo-check" in argv


def test_opencode_backend_uses_noninteractive_pure_run_and_model(tmp_dir, monkeypatch):
    captured = {}

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        return SimpleNamespace(
            returncode=0,
            stdout='OpenCode\n{"main_cell_type":"T cell"}',
            stderr="",
        )

    monkeypatch.setattr(subprocess, "run", fake_run)
    backend = AgentCLIBackend(
        "opencode",
        model="opencode/ling-3.0-flash-fin-free",
    )
    output = backend.run(
        "annotate this cluster",
        tmp_dir / "prompt.md",
        tmp_dir,
        {"agent_output_file": ""},
    )

    assert output.endswith('{"main_cell_type":"T cell"}')
    assert captured["argv"][1:3] == ["run", "--pure"]
    assert "annotate this cluster" in captured["argv"]
    assert captured["argv"][-2:] == ["--model", "opencode/ling-3.0-flash-fin-free"]


def test_annotation_prompt_uses_cassia_analysis_structure():
    prompt = build_annotation_prompt(
        MarkerCluster(cluster_id="0", markers=["CD3D", "CD3E", "TRAC"]),
        tissue="large intestine",
        species="human",
    )
    assert "professional computational biologist" in prompt
    assert "key functional or pathway markers" in prompt
    assert "key cell type markers" in prompt
    assert "Cross-check" in prompt
    assert "top three most probable sub-cell types" in prompt
    assert '"main_cell_type"' in prompt
    assert "FINAL ANNOTATION COMPLETED" not in prompt


def test_v1_prompt_and_validator_use_original_cassia_text():
    cluster = MarkerCluster(cluster_id="0", markers=["CD3D", "CD3E", "TRAC"])
    prompt = build_annotation_prompt_v1(
        cluster,
        tissue="blood",
        species="human",
    )
    assert "rewarded $10000" in prompt
    assert "FINAL ANNOTATION COMPLETED" in prompt
    assert "Your task is to annotate a single-cell human dataset from blood tissue" in prompt
    assert '"main_cell_type"' in prompt
    assert '"confidence": <integer 0-10>' in prompt

    validator_prompt = build_validator_prompt(
        annotation_response="T cell supported by CD3D and TRAC",
        marker_list=cluster.markers,
        tissue="blood",
        involvement="v1",
    )
    assert "You are an expert biologist specializing in single-cell analysis" in validator_prompt
    assert "Please validate the following annotation result" in validator_prompt
    assert "Marker List: CD3D, CD3E, TRAC" in validator_prompt


def test_validated_annotation_reuses_saved_initial_response():
    calls = []

    def fake_call(stage, attempt, prompt):
        calls.append((stage, attempt, prompt))
        if stage == "validation" and attempt == 1:
            return "VALIDATION FAILED: B-cell markers are absent."
        if stage == "annotation" and attempt == 2:
            assert "Previous annotation attempt failed validation" in prompt
            return '{"main_cell_type":"T cell","sub_cell_types":["CD3 T cell"]}'
        return "VALIDATION PASSED"

    run = run_validated_annotation(
        initial_prompt="original saved prompt",
        initial_response='{\"main_cell_type\":\"B cell\",\"sub_cell_types\":[\"B cell\"]}',
        marker_list=["CD3D", "CD3E", "TRAC"],
        tissue="blood",
        call_agent=fake_call,
    )
    assert run.validation_passed is True
    assert run.validation_attempts == 2
    assert run.history[0]["reused"] is True
    assert [(stage, attempt) for stage, attempt, _ in calls] == [
        ("validation", 1),
        ("annotation", 2),
        ("validation", 2),
    ]
    assert '"T cell"' in run.final_response


def test_self_reflect_v2_validator_and_revision_are_adversarial():
    validator_prompt = build_validator_prompt(
        annotation_response="Final call: B cell because CD3D is an immune marker.",
        marker_list=["CD3D", "CD3E", "TRAC", "IL7R"],
        tissue="blood",
        additional_info="human PBMC",
        involvement="self-reflect-v2",
        species="human",
    )
    assert "try to\nfalsify" in validator_prompt
    assert "1. CD3D" in validator_prompt
    assert "Tissue: blood" in validator_prompt
    assert "Species: human" in validator_prompt
    assert "exact primary" in validator_prompt

    revision_prompt = build_revision_prompt(
        original_prompt="Annotate these markers.",
        previous_response="Final call: B cell",
        validation_feedback="VALIDATION FAILED; T-cell program is stronger.",
        involvement="self-reflect-v2",
        marker_list=["CD3D", "CD3E", "TRAC"],
        tissue="blood",
    )
    assert "fresh annotation" in revision_prompt
    assert "do not merely defend" in revision_prompt
    assert "Ranked markers: CD3D, CD3E, TRAC" in revision_prompt


def test_cli_shell_backend(tmp_dir):
    marker_path = tmp_dir / "shell_markers.csv"
    run_dir = tmp_dir / "shell_run"
    fake_agent = tmp_dir / "fake_agent.py"
    marker_path.write_text(
        "cluster,markers\n"
        '0,"CD3D, CD3E, TRAC"\n',
        encoding="utf-8",
    )
    fake_agent.write_text(
        "import json\n"
        "print(json.dumps({\n"
        "    'main_cell_type': 'T cell',\n"
        "    'sub_cell_types': ['CD3 T cell'],\n"
        "    'possible_mixed_cell_types': [],\n"
        "    'confidence': 'high',\n"
        "    'evidence': 'CD3D, CD3E, and TRAC support a T cell annotation.'\n"
        "}));\n",
        encoding="utf-8",
    )
    code = main([
        "annotate",
        "--input",
        str(marker_path),
        "--out",
        str(run_dir),
        "--backend",
        "shell",
        "--command-template",
        f"{sys.executable} {fake_agent}",
    ])
    assert code == 0
    results = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
    assert results["0"]["analysis_result"]["main_cell_type"] == "T cell"
    assert results["0"]["analysis_result"]["schema_version"] == ANNOTATION_SCHEMA_VERSION
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["result_schema_version"] == ANNOTATION_SCHEMA_VERSION
    assert (run_dir / "summary.csv").exists()
    assert (run_dir / "report.md").exists()
    assert (run_dir / "report.html").exists()
    assert "CASSIA Agent CLI Annotation Report" in (run_dir / "report.html").read_text(encoding="utf-8")
    assert main(["report", str(run_dir)]) == 0


def test_cli_validated_shell_backend(tmp_dir):
    marker_path = tmp_dir / "validated_markers.csv"
    run_dir = tmp_dir / "validated_run"
    fake_agent = tmp_dir / "fake_validated_agent.py"
    marker_path.write_text(
        "cluster,markers\n"
        '0,"CD3D, CD3E, TRAC"\n',
        encoding="utf-8",
    )
    fake_agent.write_text(
        "import json, pathlib, sys\n"
        "prompt = pathlib.Path(sys.argv[1]).read_text(encoding='utf-8')\n"
        "if 'Please validate the following annotation result' in prompt:\n"
        "    print('VALIDATION PASSED' if '\"main_cell_type\": \"T cell\"' in prompt else "
        "'VALIDATION FAILED: marker evidence supports a T cell, not a B cell.')\n"
        "elif 'Previous annotation attempt failed validation' in prompt:\n"
        "    print(json.dumps({'main_cell_type': 'T cell', 'sub_cell_types': ['CD3 T cell'], "
        "'possible_mixed_cell_types': [], 'confidence': 'high', "
        "'evidence': 'CD3D, CD3E, and TRAC support T cells.'}))\n"
        "else:\n"
        "    print(json.dumps({'main_cell_type': 'B cell', 'sub_cell_types': ['B cell'], "
        "'possible_mixed_cell_types': [], 'confidence': 'low', "
        "'evidence': 'Initial deliberately incorrect result.'}))\n",
        encoding="utf-8",
    )
    code = main([
        "annotate",
        "--input",
        str(marker_path),
        "--out",
        str(run_dir),
        "--backend",
        "shell",
        "--command-template",
        f"{sys.executable} {fake_agent} {{prompt_file}}",
        "--workflow",
        "validated",
        "--prompt-version",
        "v1",
    ])
    assert code == 0
    results = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
    result = results["0"]
    assert result["analysis_result"]["main_cell_type"] == "T cell"
    assert result["validation_passed"] is True
    assert result["validation_attempts"] == 2
    history_path = Path(result["validation_history_file"])
    assert history_path.exists()
    history = json.loads(history_path.read_text(encoding="utf-8"))
    assert [item["stage"] for item in history["history"]] == [
        "annotation", "validation", "annotation", "validation"
    ]
    html_report = (run_dir / "report.html").read_text(encoding="utf-8")
    assert "CASSIA Agent CLI Annotation Report" in html_report
    assert "Validation PASSED" in html_report


def test_boost_query_marker_genes(tmp_dir):
    marker_path = tmp_dir / "raw_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "0,CD3D,2.5,0.91,0.04,1e-20\n"
        "0,MS4A1,-0.4,0.02,0.60,0.02\n"
        "1,MS4A1,3.2,0.88,0.03,1e-30\n",
        encoding="utf-8",
    )
    result = query_marker_genes(marker_path, ["CD3D", "MS4A1", "NKG7"], cluster="0")
    assert result[result["gene"] == "CD3D"]["avg_log2FC"].iloc[0] == 2.5
    assert result[result["gene"] == "MS4A1"]["cluster"].iloc[0] == "0"
    assert result[result["gene"] == "NKG7"]["status"].iloc[0] == "not_found"


def test_boost_gene_queries_are_unlimited_by_default():
    genes = [f"GENE{index}" for index in range(40)]
    text = f"<check_genes>{','.join(genes)}</check_genes>"
    assert extract_check_genes(text) == genes
    assert extract_check_genes(text, max_genes=20) == genes[:20]


def test_boost_gene_query_parser_ignores_explanatory_tag_mentions():
    text = (
        "Re-evaluation after the previous `<check_genes>` panel.\n"
        "The evidence supports another targeted request.\n"
        "<check_genes>NCR1,NCR3,TRGC1,TRGC2</check_genes>"
    )
    assert extract_check_genes(text) == ["NCR1", "NCR3", "TRGC1", "TRGC2"]


def test_fused_boost_prompt_has_no_prior_annotation_or_default_gene_cap():
    prompt = build_boost_prompt(
        cluster="0",
        major_cluster_info="cells in human lung",
        top_markers=["EPCAM", "KRT19", "SCGB1A1"],
        mode="fused",
    )
    assert "No prior annotation is available or trusted" in prompt
    assert "Analyze functional markers, cell-type markers" in prompt
    assert "original CASSIA functional-marker" not in prompt
    assert "Given the following cluster markers" not in prompt
    assert "There is no numerical gene cap" in prompt
    assert "Original CASSIA annotation:" not in prompt
    assert "Do not output final JSON before receiving marker-query results" in prompt


def test_fused_boost_legacy_v2_remains_available_for_reproduction():
    prompt = build_boost_prompt(
        cluster="0",
        major_cluster_info="cells in human lung",
        top_markers=["EPCAM", "KRT19", "SCGB1A1"],
        mode="fused",
        prompt_variant="v2",
    )
    assert "original CASSIA functional-marker" in prompt
    assert "CASSIA ACTIVE-EVIDENCE EXTENSION" in prompt


def test_fused_boost_v3_is_hierarchy_first_and_answer_agnostic():
    prompt = build_boost_prompt(
        cluster="test cluster",
        major_cluster_info="cells in human tissue",
        top_markers=["GENEA", "GENEB", "GENEC"],
        mode="fused",
        prompt_variant="v3",
    )
    assert "HIERARCHICAL CALIBRATION" in prompt
    assert "broad lineage/general cell type" in prompt
    assert "strongest alternative" in prompt
    assert "There is no numerical gene cap" in prompt
    assert "expected_cell_type" not in prompt
    assert "ground_truth" not in prompt


def test_fused_boost_v14_keeps_v2_initial_prompt_and_adds_final_gate():
    kwargs = {
        "cluster": "test cluster",
        "major_cluster_info": "cells in human tissue",
        "top_markers": ["GENEA", "GENEB", "GENEC"],
        "mode": "fused",
    }
    assert build_boost_prompt(**kwargs, prompt_variant="v14") == build_boost_prompt(
        **kwargs, prompt_variant="v2"
    )
    followup = build_boost_followup_prompt(
        transcript="prior transcript",
        query_text="returned statistics",
        is_final_round=True,
        prompt_variant="v14",
    )
    assert "Final sufficiency gate" in followup
    assert "coherent reciprocal program" in followup
    assert "outside the supplied dataset context" in followup
    assert "expected_cell_type" not in followup
    assert "ground_truth" not in followup


def test_candidate_boost_builds_auditable_three_candidate_tournament():
    prompt = build_boost_prompt(
        cluster="test cluster",
        major_cluster_info="cells in human blood",
        top_markers=["CD3D", "NKG7", "TRAC"],
        mode="fused",
        prompt_variant="candidate",
        candidate_count=3,
    )
    assert "propose exactly 3" in prompt
    assert "<candidate_set>" in prompt
    assert "candidate_audit" in prompt
    assert "There is no numerical gene cap" in prompt
    assert "ground_truth" not in prompt

    candidates = extract_candidate_set(
        '<candidate_set>{"candidates":['
        '{"label":"T cell","broad_lineage":"lymphoid","positive_markers":["CD3D"]},'
        '{"label":"NK cell","exclusion_markers":["TRAC"]},'
        '{"label":"NKT cell"}'
        ']}</candidate_set>\n<check_genes>CD3D,NKG7,TRAC</check_genes>'
    )
    assert [item["label"] for item in candidates] == ["T cell", "NK cell", "NKT cell"]
    followup = build_boost_followup_prompt(
        transcript="candidate slate",
        query_text="returned statistics",
        is_final_round=False,
        prompt_variant="candidate",
    )
    assert "explicit tournament" in followup
    assert "do not simply defend the initial rank 1" in followup


def test_experimental_fused_prompts_are_answer_agnostic_and_mechanistically_distinct():
    common = {
        "cluster": "test cluster",
        "major_cluster_info": "cells in human lung",
        "top_markers": ["EPCAM", "KRT8", "SCGB1A1"],
        "mode": "fused",
    }
    open_world = build_boost_prompt(**common, prompt_variant="open_world")
    program_first = build_boost_prompt(**common, prompt_variant="program_first")
    assert "OPEN-WORLD probe" in open_world
    assert "mandatory null-slate reconstruction" in open_world
    assert "phenotype-blind" in program_first
    assert "do not emit or rank any cell-type" in program_first
    assert "<program_plan>" in program_first
    for prompt in (open_world, program_first):
        assert "<check_genes>" in prompt
        assert "final_sub_cell_type" in prompt
        assert "expected_cell_type" not in prompt
        assert "ground_truth" not in prompt

    open_followup = build_boost_followup_prompt(
        "prior transcript", "returned statistics", False, "open_world"
    )
    program_followup = build_boost_followup_prompt(
        "prior transcript", "returned statistics", False, "program_first"
    )
    assert "without privileging the initial shortlist" in open_followup
    assert "map coherent stable programs to lineage first" in program_followup


def test_branch_search_fuses_breadth_and_depth_without_answer_leakage():
    prompt = build_boost_prompt(
        cluster="test cluster",
        major_cluster_info="cells in human lung",
        top_markers=["EPCAM", "KRT8", "SCGB1A1"],
        mode="fused",
        prompt_variant="branch_search",
    )
    assert "exactly three genuinely distinct" in prompt
    assert "BREADTH ROUND" in prompt
    assert "DEPTH ROUND" in prompt
    assert "OPEN BRANCH RULE" in prompt
    assert "Complete both the breadth and depth evidence rounds" in prompt
    assert "There is no numerical gene cap" in prompt
    assert "expected_cell_type" not in prompt
    assert "ground_truth" not in prompt

    followup = build_boost_followup_prompt(
        "prior transcript", "returned statistics", False, "branch_search"
    )
    assert "mutable branch ledger" in followup
    assert "distinct depth panel" in followup
    assert "both breadth and depth marker-query rounds" in followup


def test_signature_tool_prompts_preserve_v2_and_require_answer_blind_requests():
    common = {
        "cluster": "test cluster",
        "major_cluster_info": "cells in human lung",
        "top_markers": ["EPCAM", "KRT8", "SCGB1A1"],
        "mode": "fused",
    }
    baseline = build_boost_prompt(**common, prompt_variant="v2")
    gsea = build_boost_prompt(**common, prompt_variant="gsea_tool")
    ucell = build_boost_prompt(**common, prompt_variant="ucell_tool")

    assert gsea.startswith(baseline)
    assert ucell.startswith(baseline)
    assert "<gsea_request>" in gsea
    assert "weighted preranked gsea" in gsea.lower()
    assert "<ucell_request>" in ucell
    assert "per-cell Mann-Whitney rank score" in ucell
    for prompt in (gsea, ucell):
        assert "Propose 2-4 genuinely competing, answer-agnostic signatures" in prompt
        assert "Never infer the hidden target label" in prompt
        assert "expected_cell_type" not in prompt
        assert "ground_truth" not in prompt

    gsea_followup = build_boost_followup_prompt(
        "prior transcript", "returned statistics", False, "gsea_tool"
    )
    ucell_followup = build_boost_followup_prompt(
        "prior transcript", "returned statistics", False, "ucell_tool"
    )
    assert "Interpret GSEA at the program level" in gsea_followup
    assert "Interpret UCell distributions" in ucell_followup


def test_codex_backend_pins_reasoning_effort_without_editing_user_config(tmp_path, monkeypatch):
    captured = {}

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        return SimpleNamespace(returncode=0, stdout='{"ok":true}', stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    backend = AgentCLIBackend(
        "codex-cli",
        model="gpt-5.6-luna",
        reasoning_effort="max",
    )
    output = backend.run(
        "judge this",
        tmp_path / "prompt.md",
        tmp_path,
        {"agent_output_file": ""},
    )
    assert output == '{"ok":true}'
    assert "--model" in captured["argv"]
    assert "gpt-5.6-luna" in captured["argv"]
    assert "--config" in captured["argv"]
    assert 'model_reasoning_effort="max"' in captured["argv"]


def test_stable_judge_scores_only_top1_and_keeps_lower_ranks_diagnostic():
    expected = [{
        "case_id": "case-1",
        "mode": "pred_A",
        "source": "",
        "expected_label": "CD8 effector memory T cell",
        "expected_terms": "",
        "marker_genes": "CD3D, CD8A, GZMK",
        "prediction": {
            "missing_output": False,
            "main_cell_type": "T cell",
            "top1_subtype": "naive CD8 T cell",
            "top2_subtype": "CD8 effector memory T cell",
            "top3_subtype": "NK cell",
        },
    }]
    payload = {
        "protocol_version": "stable-judge-v1.1",
        "judgments": [{
            "case_id": "case-1",
            "mode": "pred_A",
            "top1_verdict": "wrong",
            "matched_rank": 2,
            "lineage_score": 2,
            "subtype_score": 0,
            "state_score": 0,
            "marker_support_score": 0,
            "major_error": "wrong_subtype",
            "rationale": "The exact truth appears only at rank 2.",
        }],
    }
    scored = normalize_judgments(payload, expected).iloc[0]
    assert scored["top1_verdict"] == "wrong"
    assert bool(scored["lower_rank_only_recovery"]) is True
    assert bool(scored["judge_pass"]) is False
    assert "top1_subtype is the ONLY ranked subtype" in build_judge_prompt(expected)


def test_fused_schema_preserves_named_displaced_candidate_without_failing():
    normalized = normalize_fused_boost_payload({
        "final_cell_type": "Endocrine cell",
        "final_sub_cell_type": "Epsilon cell",
        "changed_from_original": "Alpha cell",
    })
    assert normalized["changed_from_original"] is None
    assert normalized["changed_from_candidate"] == "Alpha cell"


def test_boost_top_markers_filter_low_pct_artifacts(tmp_dir):
    marker_path = tmp_dir / "raw_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "0,RARE_ARTIFACT,5.0,0.001,0.0,1e-6\n"
        "0,CD3D,2.5,0.91,0.04,1e-20\n"
        "0,TRAC,2.1,0.88,0.05,1e-18\n",
        encoding="utf-8",
    )
    genes, top_rows = get_cluster_top_markers(marker_path, cluster="0", n_genes=2)
    assert genes == ["CD3D", "TRAC"]
    assert "RARE_ARTIFACT" not in top_rows["gene"].tolist()


def test_cli_boost_query(tmp_dir):
    marker_path = tmp_dir / "raw_markers.csv"
    out_path = tmp_dir / "query.json"
    marker_path.write_text(
        "group,names,logfoldchanges,pvals_adj\n"
        "0,CD3D,2.5,1e-20\n"
        "0,TRAC,2.1,1e-10\n",
        encoding="utf-8",
    )
    code = main([
        "boost",
        "query",
        "--markers",
        str(marker_path),
        "--cluster",
        "0",
        "--genes",
        "CD3D, TRAC",
        "--format",
        "json",
        "--out",
        str(out_path),
    ])
    assert code == 0
    assert '"gene": "CD3D"' in out_path.read_text(encoding="utf-8")
    assert parse_gene_args(["CD3D, TRAC", "CD3D"]) == ["CD3D", "TRAC"]


def test_cli_boost_run_shell_backend(tmp_dir):
    run_dir = tmp_dir / "boost_source_run"
    run_dir.mkdir()
    (run_dir / "results.json").write_text(
        json.dumps({
            "0": {
                "analysis_result": {
                    "main_cell_type": "T cell",
                    "sub_cell_types": ["CD3 T cell"],
                    "confidence": "medium",
                }
            }
        }),
        encoding="utf-8",
    )
    marker_path = tmp_dir / "raw_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "0,CD3D,2.5,0.91,0.04,1e-20\n"
        "0,TRAC,2.1,0.88,0.05,1e-18\n"
        "0,MS4A1,-0.4,0.02,0.60,0.02\n",
        encoding="utf-8",
    )
    fake_agent = tmp_dir / "fake_boost_agent.py"
    fake_agent.write_text(
        "import json, pathlib, sys\n"
        "prompt = pathlib.Path(sys.argv[1]).read_text()\n"
        "if 'Latest CASSIA marker query results' in prompt:\n"
        "    print(json.dumps({\n"
        "        'final_cell_type': 'T cell',\n"
        "        'final_sub_cell_type': 'CD3 T cell',\n"
        "        'confidence': 'high',\n"
        "        'changed_from_original': False,\n"
        "        'checked_genes': ['CD3D', 'TRAC', 'MS4A1'],\n"
        "        'supporting_markers': ['CD3D', 'TRAC'],\n"
        "        'refuting_markers': ['MS4A1'],\n"
        "        'alternatives': [],\n"
        "        'evidence': 'CD3D and TRAC support T cells; MS4A1 argues against B cells.'\n"
        "    }))\n"
        "else:\n"
        "    print('Evaluate T cell versus B cell. <check_genes>CD3D,TRAC,MS4A1</check_genes>\\n' + json.dumps({\n"
        "        'final_cell_type': 'B cell',\n"
        "        'final_sub_cell_type': 'Premature answer that must not bypass the query',\n"
        "        'confidence': 'low'\n"
        "    }))\n",
        encoding="utf-8",
    )
    code = main([
        "boost",
        "run",
        "--run",
        str(run_dir),
        "--markers",
        str(marker_path),
        "--cluster",
        "0",
        "--backend",
        "shell",
        "--command-template",
        f"{sys.executable} {fake_agent} {{prompt_file}}",
        "--iterations",
        "3",
    ])
    assert code == 0
    boost_dir = run_dir / "boost" / "0"
    final_result = json.loads((boost_dir / "final.json").read_text(encoding="utf-8"))
    manifest = json.loads((boost_dir / "boost_manifest.json").read_text(encoding="utf-8"))
    assert final_result["final_cell_type"] == "T cell"
    assert final_result["schema_version"] == ANNOTATION_SCHEMA_VERSION
    assert final_result["main_cell_type"] == "T cell"
    assert final_result["cluster_id"] == "0"
    assert "CD3D" in manifest["checked_genes"]
    assert (boost_dir / "queries" / "round_001.csv").exists()
    assert (boost_dir / "transcript.md").exists()
    assert (boost_dir / "summary.html").exists()
    assert (boost_dir / "summary_tags.txt").exists()
    assert manifest["html_report"].endswith("summary.html")
    assert manifest["result_schema_version"] == ANNOTATION_SCHEMA_VERSION
    assert manifest["outputs"]["final_json"].endswith("final.json")
    assert "CASSIA Cell Type Annotation Summary" in (boost_dir / "summary.html").read_text(encoding="utf-8")
    report_code = main(["report", str(boost_dir)])
    assert report_code == 0
    assert "CASSIA Cell Type Annotation Summary" in (boost_dir / "summary.html").read_text(encoding="utf-8")


def test_cli_boost_fused_dry_run_does_not_require_prior_annotation(tmp_dir):
    run_dir = tmp_dir / "fused_boost_run"
    run_dir.mkdir()
    marker_path = tmp_dir / "fused_boost_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "0,EPCAM,3.0,0.90,0.05,1e-20\n"
        "0,KRT19,2.5,0.80,0.10,1e-18\n"
        "0,SCGB1A1,2.1,0.70,0.03,1e-15\n",
        encoding="utf-8",
    )
    code = main([
        "boost",
        "run",
        "--run",
        str(run_dir),
        "--markers",
        str(marker_path),
        "--cluster",
        "0",
        "--mode",
        "fused",
        "--dry-run",
    ])
    assert code == 0
    boost_dir = run_dir / "boost" / "0"
    prompt = (boost_dir / "prompts" / "round_001.md").read_text(encoding="utf-8")
    manifest = json.loads((boost_dir / "boost_manifest.json").read_text(encoding="utf-8"))
    assert "No prior annotation is available or trusted" in prompt
    assert "There is no numerical gene cap" in prompt
    assert manifest["mode"] == "fused"


def test_cli_annotate_fused_boost_mode_dry_run(tmp_dir):
    marker_path = tmp_dir / "annotate_fused_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "0,EPCAM,3.0,0.90,0.05,1e-20\n"
        "0,KRT19,2.5,0.80,0.10,1e-18\n",
        encoding="utf-8",
    )
    output_dir = tmp_dir / "annotate_fused_run"
    code = main([
        "annotate",
        "--input",
        str(marker_path),
        "--cluster",
        "0",
        "--mode",
        "fused-boost",
        "--backend",
        "codex-cli",
        "--dry-run",
        "--out",
        str(output_dir),
    ])
    assert code == 0
    manifest = json.loads((output_dir / "boost_manifest.json").read_text(encoding="utf-8"))
    assert manifest["mode"] == "fused"
    assert manifest["cluster"] == "0"
    assert (output_dir / "prompts" / "round_001.md").exists()
    assert manifest["annotation_source"] is None


def test_cli_boost_auto_plan_selects_uncertain_clusters(tmp_dir):
    run_dir = tmp_dir / "boost_auto_plan_run"
    run_dir.mkdir()
    (run_dir / "results.json").write_text(
        json.dumps({
            "0": {
                "analysis_result": {
                    "main_cell_type": "T cell",
                    "sub_cell_types": ["CD3 T cell"],
                    "confidence": "high",
                    "possible_mixed_cell_types": [],
                    "evidence": "CD3D, CD3E, and TRAC strongly support T cells.",
                }
            },
            "1": {
                "analysis_result": {
                    "main_cell_type": "Macrophage",
                    "sub_cell_types": ["TAM-like macrophage"],
                    "confidence": "low",
                    "possible_mixed_cell_types": [],
                    "evidence": "Ambiguous macrophage versus dendritic marker support.",
                }
            },
            "2": {
                "analysis_result": {
                    "main_cell_type": "B cell",
                    "sub_cell_types": ["activated B cell"],
                    "confidence": "medium",
                    "possible_mixed_cell_types": ["plasma cell"],
                    "evidence": "Mixed immunoglobulin and MS4A1 signal.",
                }
            },
        }),
        encoding="utf-8",
    )
    marker_path = tmp_dir / "boost_auto_plan_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "1,LST1,2.8,0.90,0.05,1e-20\n"
        "1,IL1B,2.2,0.70,0.04,1e-12\n"
        "2,MS4A1,3.1,0.88,0.03,1e-30\n"
        "2,MZB1,2.4,0.64,0.02,1e-18\n",
        encoding="utf-8",
    )
    code = main([
        "boost",
        "auto",
        "--run",
        str(run_dir),
        "--markers",
        str(marker_path),
        "--backend",
        "shell",
        "--plan-only",
        "--max-clusters",
        "2",
    ])
    assert code == 0
    auto_dir = run_dir / "boost" / "_auto"
    manifest = json.loads((auto_dir / "boost_auto_manifest.json").read_text(encoding="utf-8"))
    selected_clusters = {candidate["cluster"] for candidate in manifest["selected_candidates"]}
    assert manifest["status"] == "plan-only"
    assert selected_clusters == {"1", "2"}
    assert (auto_dir / "auto_summary.csv").exists()
    assert (auto_dir / "auto_report.html").exists()


def test_cli_boost_auto_shell_backend_runs_selected_cluster(tmp_dir):
    run_dir = tmp_dir / "boost_auto_run"
    run_dir.mkdir()
    (run_dir / "results.json").write_text(
        json.dumps({
            "1": {
                "analysis_result": {
                    "main_cell_type": "Macrophage",
                    "sub_cell_types": ["TAM-like macrophage"],
                    "confidence": "low",
                    "possible_mixed_cell_types": [],
                    "evidence": "Ambiguous myeloid markers require boost review.",
                }
            }
        }),
        encoding="utf-8",
    )
    marker_path = tmp_dir / "boost_auto_markers.csv"
    marker_path.write_text(
        "cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj\n"
        "1,LST1,2.8,0.90,0.05,1e-20\n"
        "1,IL1B,2.2,0.70,0.04,1e-12\n"
        "1,MS4A7,2.0,0.62,0.04,1e-10\n",
        encoding="utf-8",
    )
    fake_agent = tmp_dir / "fake_boost_auto_agent.py"
    fake_agent.write_text(
        "import json\n"
        "print(json.dumps({\n"
        "    'final_cell_type': 'Macrophage',\n"
        "    'final_sub_cell_type': 'Inflammatory macrophage',\n"
        "    'confidence': 'high',\n"
        "    'changed_from_original': False,\n"
        "    'checked_genes': ['LST1', 'IL1B', 'MS4A7'],\n"
        "    'supporting_markers': ['LST1', 'IL1B', 'MS4A7'],\n"
        "    'refuting_markers': [],\n"
        "    'alternatives': [],\n"
        "    'evidence': 'LST1, IL1B, and MS4A7 support an inflammatory macrophage annotation.'\n"
        "}));\n",
        encoding="utf-8",
    )
    code = main([
        "boost",
        "auto",
        "--run",
        str(run_dir),
        "--markers",
        str(marker_path),
        "--backend",
        "shell",
        "--command-template",
        f"{sys.executable} {fake_agent} {{prompt_file}}",
        "--max-clusters",
        "1",
        "--iterations",
        "2",
    ])
    assert code == 0
    auto_dir = run_dir / "boost" / "_auto"
    cluster_dir = run_dir / "boost" / "1"
    manifest = json.loads((auto_dir / "boost_auto_manifest.json").read_text(encoding="utf-8"))
    final_result = json.loads((cluster_dir / "final.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "completed"
    assert manifest["results"][0]["status"] == "completed"
    assert final_result["final_sub_cell_type"] == "Inflammatory macrophage"
    assert (auto_dir / "auto_summary.csv").exists()
    assert (auto_dir / "auto_report.html").exists()
    assert (cluster_dir / "summary.html").exists()


def test_cli_subcluster_dry_run(tmp_dir):
    marker_path = tmp_dir / "subcluster_markers.csv"
    out_dir = tmp_dir / "subcluster_dry_run"
    marker_path.write_text(
        "cluster,markers\n"
        'exhausted,"HAVCR2, TIGIT, LAG3, CXCL13"\n'
        'memory,"IL7R, CCR7, TCF7, LEF1"\n',
        encoding="utf-8",
    )
    code = main([
        "subcluster",
        "run",
        "--markers",
        str(marker_path),
        "--major-cluster-info",
        "CD8 T cell",
        "--backend",
        "shell",
        "--dry-run",
        "--out",
        str(out_dir),
    ])
    assert code == 0
    manifest = json.loads((out_dir / "subcluster_manifest.json").read_text(encoding="utf-8"))
    prompt = (out_dir / "prompts" / "subcluster_prompt.md").read_text(encoding="utf-8")
    assert manifest["status"] == "dry-run"
    assert "Parent cluster context: CD8 T cell" in prompt
    assert "Subcluster exhausted" in prompt
    assert not (out_dir / "subcluster_results.csv").exists()


def test_cli_subcluster_run_shell_backend(tmp_dir):
    marker_path = tmp_dir / "subcluster_markers_run.csv"
    out_dir = tmp_dir / "subcluster_run"
    fake_agent = tmp_dir / "fake_subcluster_agent.py"
    marker_path.write_text(
        "cluster,markers\n"
        'exhausted,"HAVCR2, TIGIT, LAG3, CXCL13"\n'
        'memory,"IL7R, CCR7, TCF7, LEF1"\n',
        encoding="utf-8",
    )
    fake_agent.write_text(
        "import json, pathlib, sys\n"
        "prompt = pathlib.Path(sys.argv[1]).read_text()\n"
        "assert 'Subcluster exhausted' in prompt\n"
        "assert 'Subcluster memory' in prompt\n"
        "print(json.dumps({\n"
        "    'parent_cluster': 'CD8 T cell',\n"
        "    'subclusters': [\n"
        "        {\n"
        "            'cluster_id': 'exhausted',\n"
        "            'main_cell_type': 'Exhausted CD8 T cell',\n"
        "            'sub_cell_type': 'CXCL13-positive exhausted CD8 T cell',\n"
        "            'key_markers': ['HAVCR2', 'TIGIT', 'LAG3', 'CXCL13'],\n"
        "            'reason': 'HAVCR2, TIGIT, LAG3, and CXCL13 support an exhausted state.'\n"
        "        },\n"
        "        {\n"
        "            'cluster_id': 'memory',\n"
        "            'main_cell_type': 'Memory CD8 T cell',\n"
        "            'sub_cell_type': 'Central memory CD8 T cell',\n"
        "            'key_markers': ['IL7R', 'CCR7', 'TCF7', 'LEF1'],\n"
        "            'reason': 'IL7R, CCR7, TCF7, and LEF1 support central memory identity.'\n"
        "        }\n"
        "    ]\n"
        "}));\n",
        encoding="utf-8",
    )
    code = main([
        "subcluster",
        "run",
        "--markers",
        str(marker_path),
        "--major-cluster-info",
        "CD8 T cell",
        "--backend",
        "shell",
        "--command-template",
        f"{sys.executable} {fake_agent} {{prompt_file}}",
        "--out",
        str(out_dir),
    ])
    assert code == 0
    manifest = json.loads((out_dir / "subcluster_manifest.json").read_text(encoding="utf-8"))
    with (out_dir / "subcluster_results.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert manifest["status"] == "completed"
    assert rows[0]["Result ID"] == "exhausted"
    assert rows[0]["sub_cell_type"] == "CXCL13-positive exhausted CD8 T cell"
    assert rows[1]["Result ID"] == "memory"
    assert (out_dir / "raw" / "subcluster_response.txt").exists()
    assert "Subclustering Annotation Report" in (out_dir / "subcluster_report.html").read_text(encoding="utf-8")


def test_cli_consensus_summary_csvs(tmp_dir):
    run_a = tmp_dir / "consensus_run_a"
    run_b = tmp_dir / "consensus_run_b"
    run_c = tmp_dir / "consensus_run_c"
    for run_dir in (run_a, run_b, run_c):
        run_dir.mkdir()
    headers = [
        "Cluster ID",
        "Predicted General Cell Type",
        "Predicted Detailed Cell Type",
        "Confidence",
        "Evidence",
    ]
    rows_by_run = {
        run_a: [
            ["0", "T cell", "CD8 T cell", "high", "CD3D and CD8A support T cells."],
            ["1", "Macrophage", "Inflammatory macrophage", "medium", "LST1 and IL1B support macrophages."],
        ],
        run_b: [
            ["0", "T cell", "CD8 T cell", "high", "TRAC and CD8B support T cells."],
            ["1", "Dendritic cell", "cDC", "low", "FCER1A suggested dendritic identity."],
        ],
        run_c: [
            ["0", "NK cell", "NK cell", "medium", "NKG7 suggested NK identity."],
            ["1", "Macrophage", "TAM-like macrophage", "high", "C1QA and APOE support macrophages."],
        ],
    }
    for run_dir, rows in rows_by_run.items():
        with (run_dir / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(headers)
            writer.writerows(rows)

    out_path = tmp_dir / "consensus.csv"
    code = main([
        "consensus",
        "--inputs",
        str(run_a),
        str(run_b / "summary.csv"),
        str(run_c),
        "--out",
        str(out_path),
    ])
    assert code == 0
    with out_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    by_cluster = {row["cluster_id"]: row for row in rows}
    assert by_cluster["0"]["consensus_main_cell_type"] == "T cell"
    assert by_cluster["0"]["consensus_sub_cell_type"] == "CD8 T cell"
    assert by_cluster["0"]["status"] == "consensus"
    assert by_cluster["1"]["consensus_main_cell_type"] == "Macrophage"
    assert by_cluster["1"]["status"] == "partial_consensus"
    assert (tmp_dir / "consensus.html").exists()
    assert "CASSIA Consensus Report" in (tmp_dir / "consensus.html").read_text(encoding="utf-8")


def test_cli_consensus_subcluster_glob(tmp_dir):
    base_dir = tmp_dir / "subcluster_consensus"
    run_a = base_dir / "sub_a"
    run_b = base_dir / "sub_b"
    run_a.mkdir(parents=True)
    run_b.mkdir(parents=True)
    headers = ["Result ID", "main_cell_type", "sub_cell_type", "key_markers", "reason"]
    with (run_a / "subcluster_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerow(["exhausted", "Exhausted CD8 T cell", "CXCL13-positive exhausted CD8 T cell", "HAVCR2, CXCL13", ""])
        writer.writerow(["memory", "Memory CD8 T cell", "Central memory CD8 T cell", "IL7R, CCR7", ""])
    with (run_b / "subcluster_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerow(["exhausted", "Exhausted CD8 T cell", "CXCL13-positive exhausted CD8 T cell", "TIGIT, CXCL13", ""])
        writer.writerow(["memory", "Naive CD8 T cell", "Naive-like CD8 T cell", "TCF7, LEF1", ""])

    out_path = tmp_dir / "subcluster_consensus.csv"
    code = main([
        "consensus",
        "--glob",
        str(base_dir / "*"),
        "--out",
        str(out_path),
        "--no-html",
    ])
    assert code == 0
    with out_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    by_cluster = {row["cluster_id"]: row for row in rows}
    assert by_cluster["exhausted"]["status"] == "consensus"
    assert by_cluster["exhausted"]["consensus_main_cell_type"] == "Exhausted CD8 T cell"
    assert by_cluster["memory"]["status"] == "conflict"
    assert not (tmp_dir / "subcluster_consensus.html").exists()


def test_openrouter_deepseek_reasoning_and_timeout(monkeypatch):
    import CASSIA.core.llm_utils as llm_utils

    captured = {}

    class FakeResponse:
        def raise_for_status(self):
            return None

        def json(self):
            return {
                "id": "test-response",
                "model": "deepseek/deepseek-v4-flash-0731",
                "choices": [{"message": {"content": "OK"}, "finish_reason": "stop"}],
                "usage": {"prompt_tokens": 1, "completion_tokens": 1, "total_tokens": 2},
            }

    def fake_post(url, headers, data, timeout):
        captured.update({"url": url, "data": json.loads(data), "timeout": timeout})
        return FakeResponse()

    monkeypatch.setattr(llm_utils.requests, "post", fake_post)
    result = llm_utils.call_llm(
        prompt="Reply: OK",
        provider="openrouter",
        model="deepseek/deepseek-v4-flash-0731",
        api_key="test-key",
        temperature=0,
        max_tokens=32,
        reasoning={"effort": "high"},
    )

    assert result == "OK"
    assert captured["data"]["reasoning"] == {"effort": "high"}
    assert captured["timeout"] == 600


def test_direct_deepseek_uses_standard_key_and_native_thinking(monkeypatch):
    import openai
    import CASSIA.core.llm_utils as llm_utils

    captured = {}

    class FakeCompletions:
        def create(self, **kwargs):
            captured["request"] = kwargs
            return SimpleNamespace(
                id="test-response",
                model="deepseek-v4-flash",
                usage=None,
                choices=[SimpleNamespace(message=SimpleNamespace(content="OK"))],
            )

    class FakeClient:
        def __init__(self):
            self.chat = SimpleNamespace(completions=FakeCompletions())

    def fake_openai(**kwargs):
        captured["client"] = kwargs
        return FakeClient()

    monkeypatch.setenv("DEEPSEEK_API_KEY", "direct-test-key")
    monkeypatch.delenv("CUSTOMIZED_API_KEY", raising=False)
    monkeypatch.setattr(openai, "OpenAI", fake_openai)

    result = llm_utils.call_llm(
        prompt="Reply: OK",
        provider="https://api.deepseek.com",
        model="deepseek-v4-flash",
        temperature=0,
        max_tokens=32,
        reasoning={"effort": "low"},
    )

    assert result == "OK"
    assert captured["client"] == {
        "api_key": "direct-test-key",
        "base_url": "https://api.deepseek.com",
    }
    assert captured["request"]["extra_body"] == {
        "thinking": {"type": "enabled"},
        "reasoning_effort": "low",
    }


def run_all_tests():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        test_cli_help_includes_workflow_examples()
        test_cli_help_command_supports_nested_topics()
        test_cli_agent_guide_can_print_and_save(tmp_dir)
        test_canonical_result_schema_preserves_workflow_fields()
        test_cli_examples_generates_runnable_project(tmp_dir)
        test_cli_validate_preformatted_marker_list(tmp_dir)
        test_cli_validate_long_marker_table_json(tmp_dir)
        test_cli_validate_missing_ranking_columns_errors(tmp_dir)
        test_extract_json_object()
        test_load_marker_clusters_preformatted(tmp_dir)
        test_load_marker_clusters_long_table_with_gene_column_override(tmp_dir)
        test_cli_dry_run(tmp_dir)
        test_codex_backend_skips_git_repo_check(tmp_dir)
        test_annotation_prompt_uses_cassia_analysis_structure()
        test_v1_prompt_and_validator_use_original_cassia_text()
        test_validated_annotation_reuses_saved_initial_response()
        test_cli_shell_backend(tmp_dir)
        test_cli_validated_shell_backend(tmp_dir)
        test_boost_query_marker_genes(tmp_dir)
        test_boost_gene_queries_are_unlimited_by_default()
        test_boost_gene_query_parser_ignores_explanatory_tag_mentions()
        test_fused_boost_prompt_has_no_prior_annotation_or_default_gene_cap()
        test_fused_boost_legacy_v2_remains_available_for_reproduction()
        test_fused_boost_v3_is_hierarchy_first_and_answer_agnostic()
        test_fused_boost_v14_keeps_v2_initial_prompt_and_adds_final_gate()
        test_candidate_boost_builds_auditable_three_candidate_tournament()
        test_experimental_fused_prompts_are_answer_agnostic_and_mechanistically_distinct()
        test_branch_search_fuses_breadth_and_depth_without_answer_leakage()
        test_signature_tool_prompts_preserve_v2_and_require_answer_blind_requests()
        test_stable_judge_scores_only_top1_and_keeps_lower_ranks_diagnostic()
        test_fused_schema_preserves_named_displaced_candidate_without_failing()
        test_boost_top_markers_filter_low_pct_artifacts(tmp_dir)
        test_cli_boost_query(tmp_dir)
        test_cli_boost_run_shell_backend(tmp_dir)
        test_cli_boost_fused_dry_run_does_not_require_prior_annotation(tmp_dir)
        test_cli_annotate_fused_boost_mode_dry_run(tmp_dir)
        test_cli_boost_auto_plan_selects_uncertain_clusters(tmp_dir)
        test_cli_boost_auto_shell_backend_runs_selected_cluster(tmp_dir)
        test_cli_subcluster_dry_run(tmp_dir)
        test_cli_subcluster_run_shell_backend(tmp_dir)
        test_cli_consensus_summary_csvs(tmp_dir)
        test_cli_consensus_subcluster_glob(tmp_dir)
    print("CLI core tests passed")
    return True


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
