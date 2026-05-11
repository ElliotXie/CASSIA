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


sys.path.insert(0, str(Path(__file__).parent.parent.parent / "CASSIA_python"))

from CASSIA.cli import main
from CASSIA.cli.backends import AGENT_BACKENDS, render_agent_argv
from CASSIA.cli.boost import get_cluster_top_markers, parse_gene_args, query_marker_genes
from CASSIA.cli.runner import MarkerCluster, build_annotation_prompt, extract_json_object, load_marker_clusters


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
    annotate_help = _capture_help(["annotate", "--help"])
    consensus_help = _capture_help(["consensus", "--help"])
    boost_auto_help = _capture_help(["boost", "auto", "--help"])
    subcluster_help = _capture_help(["subcluster", "run", "--help"])

    assert "Common workflows:" in top_help
    assert "cassia validate markers.csv" in top_help
    assert "cassia examples --out cassia_example" in top_help
    assert "Validate marker CSV structure" in validate_help
    assert "Create marker CSVs" in examples_help
    assert "cassia annotate -i markers.csv --backend codex-cli" in top_help
    assert "Examples:" in annotate_help
    assert "cassia consensus --inputs runs/codex/summary.csv" in consensus_help
    assert "cassia boost auto --run runs/brain_codex" in boost_auto_help
    assert "cassia subcluster run --markers cd8_subcluster_markers.csv" in subcluster_help


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
    assert (run_dir / "summary.csv").exists()
    assert (run_dir / "report.md").exists()


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
        "    print('Evaluate T cell versus B cell. <check_genes>CD3D,TRAC,MS4A1</check_genes>')\n",
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
    assert "CD3D" in manifest["checked_genes"]
    assert (boost_dir / "queries" / "round_001.csv").exists()
    assert (boost_dir / "transcript.md").exists()
    assert (boost_dir / "summary.html").exists()
    assert (boost_dir / "summary_tags.txt").exists()
    assert manifest["html_report"].endswith("summary.html")
    assert "CASSIA Cell Type Annotation Summary" in (boost_dir / "summary.html").read_text(encoding="utf-8")


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


def run_all_tests():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        test_cli_help_includes_workflow_examples()
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
        test_cli_shell_backend(tmp_dir)
        test_boost_query_marker_genes(tmp_dir)
        test_boost_top_markers_filter_low_pct_artifacts(tmp_dir)
        test_cli_boost_query(tmp_dir)
        test_cli_boost_run_shell_backend(tmp_dir)
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
