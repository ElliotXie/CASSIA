"""
CASSIA Test 28: CLI Core
========================
Fast tests for the CLI run-folder and parser logic. These do not call external
agent CLIs or LLM APIs.

Usage:
    python test_cli_core.py
"""

import json
import sys
import tempfile
from pathlib import Path


sys.path.insert(0, str(Path(__file__).parent.parent.parent / "CASSIA_python"))

from CASSIA.cli import main
from CASSIA.cli.backends import AGENT_BACKENDS, render_agent_argv
from CASSIA.cli.boost import get_cluster_top_markers, parse_gene_args, query_marker_genes
from CASSIA.cli.runner import MarkerCluster, build_annotation_prompt, extract_json_object, load_marker_clusters


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


def run_all_tests():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        test_extract_json_object()
        test_load_marker_clusters_preformatted(tmp_dir)
        test_cli_dry_run(tmp_dir)
        test_codex_backend_skips_git_repo_check(tmp_dir)
        test_annotation_prompt_uses_cassia_analysis_structure()
        test_cli_shell_backend(tmp_dir)
        test_boost_query_marker_genes(tmp_dir)
        test_boost_top_markers_filter_low_pct_artifacts(tmp_dir)
        test_cli_boost_query(tmp_dir)
        test_cli_boost_run_shell_backend(tmp_dir)
    print("CLI core tests passed")
    return True


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
