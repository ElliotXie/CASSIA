"""Generate runnable example projects for the CASSIA CLI."""

from __future__ import annotations

import csv
import textwrap
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence


MARKER_ROWS = [
    {
        "cluster": "T_cell",
        "markers": "CD3D, CD3E, TRAC, IL7R, LTB, CCR7, CD27, TCF7, LEF1, SELL",
    },
    {
        "cluster": "B_cell",
        "markers": "MS4A1, CD79A, CD79B, CD74, BANK1, CD37, HLA-DRA, CD52, LTB, TCL1A",
    },
    {
        "cluster": "Macrophage",
        "markers": "LYZ, LST1, C1QA, C1QB, APOE, MS4A7, TYROBP, FCER1G, AIF1, CST3",
    },
]

RAW_MARKER_ROWS = [
    ("T_cell", "CD3D", 3.4, 0.94, 0.06, "1e-40"),
    ("T_cell", "CD3E", 3.1, 0.92, 0.05, "1e-38"),
    ("T_cell", "TRAC", 2.9, 0.90, 0.08, "1e-35"),
    ("T_cell", "IL7R", 2.2, 0.70, 0.12, "1e-18"),
    ("T_cell", "NKG7", 0.4, 0.18, 0.10, "0.02"),
    ("B_cell", "MS4A1", 3.2, 0.88, 0.03, "1e-36"),
    ("B_cell", "CD79A", 3.0, 0.85, 0.04, "1e-34"),
    ("B_cell", "CD79B", 2.8, 0.81, 0.05, "1e-30"),
    ("B_cell", "BANK1", 2.0, 0.62, 0.03, "1e-16"),
    ("Macrophage", "LYZ", 3.5, 0.95, 0.20, "1e-45"),
    ("Macrophage", "LST1", 3.1, 0.90, 0.10, "1e-35"),
    ("Macrophage", "C1QA", 2.7, 0.78, 0.02, "1e-30"),
    ("Macrophage", "APOE", 2.4, 0.70, 0.05, "1e-22"),
]

SUBCLUSTER_ROWS = [
    {
        "cluster": "exhausted",
        "markers": "HAVCR2, TIGIT, LAG3, CXCL13, PDCD1, TOX, CTLA4, ENTPD1",
    },
    {
        "cluster": "memory",
        "markers": "IL7R, CCR7, TCF7, LEF1, SELL, CD27, LTB, MAL",
    },
    {
        "cluster": "effector",
        "markers": "GZMB, PRF1, NKG7, CCL5, IFNG, GNLY, CST7, GZMA",
    },
]

SUMMARY_HEADERS = [
    "Cluster ID",
    "Predicted General Cell Type",
    "Predicted Detailed Cell Type",
    "Possible Mixed Cell Types",
    "Marker Number",
    "Marker List",
    "Backend",
    "Confidence",
    "Evidence",
]

SUMMARY_ROWS_BY_SOURCE = {
    "codex": [
        [
            "T_cell",
            "T cell",
            "Naive/central memory T cell",
            "",
            "10",
            MARKER_ROWS[0]["markers"],
            "codex-cli",
            "high",
            "CD3D, CD3E, TRAC, IL7R, CCR7, and TCF7 support a T cell memory-like identity.",
        ],
        [
            "B_cell",
            "B cell",
            "MS4A1-positive naive B cell",
            "",
            "10",
            MARKER_ROWS[1]["markers"],
            "codex-cli",
            "high",
            "MS4A1, CD79A, CD79B, and BANK1 support B cell identity.",
        ],
        [
            "Macrophage",
            "Macrophage",
            "C1Q/APOE macrophage",
            "",
            "10",
            MARKER_ROWS[2]["markers"],
            "codex-cli",
            "medium",
            "LYZ, LST1, C1QA, C1QB, APOE, and TYROBP support macrophage identity.",
        ],
    ],
    "claude": [
        [
            "T_cell",
            "T cell",
            "Central memory T cell",
            "",
            "10",
            MARKER_ROWS[0]["markers"],
            "claude-cli",
            "high",
            "CD3D, CD3E, TRAC, CCR7, IL7R, and TCF7 support a central memory T cell label.",
        ],
        [
            "B_cell",
            "B cell",
            "Naive B cell",
            "",
            "10",
            MARKER_ROWS[1]["markers"],
            "claude-cli",
            "high",
            "MS4A1, CD79A, CD79B, CD74, and BANK1 support a B cell annotation.",
        ],
        [
            "Macrophage",
            "Macrophage",
            "APOE/C1Q macrophage",
            "",
            "10",
            MARKER_ROWS[2]["markers"],
            "claude-cli",
            "medium",
            "LYZ, LST1, C1QA, C1QB, APOE, and MS4A7 support macrophage identity.",
        ],
    ],
}


def _write_csv(path: Path, headers: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def _write_dict_csv(path: Path, headers: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _toy_agent_source() -> str:
    return textwrap.dedent(
        """\
        #!/usr/bin/env python3
        import json
        import pathlib
        import sys


        def print_json(payload):
            print(json.dumps(payload, ensure_ascii=False))


        prompt = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")

        if "Subcluster exhausted" in prompt or '"subclusters"' in prompt:
            print_json({
                "parent_cluster": "CD8 T cell",
                "subclusters": [
                    {
                        "cluster_id": "exhausted",
                        "main_cell_type": "Exhausted CD8 T cell",
                        "sub_cell_type": "CXCL13-positive exhausted CD8 T cell",
                        "key_markers": ["HAVCR2", "TIGIT", "LAG3", "CXCL13"],
                        "reason": "HAVCR2, TIGIT, LAG3, and CXCL13 support an exhausted T cell state."
                    },
                    {
                        "cluster_id": "memory",
                        "main_cell_type": "Memory CD8 T cell",
                        "sub_cell_type": "Central memory CD8 T cell",
                        "key_markers": ["IL7R", "CCR7", "TCF7", "LEF1"],
                        "reason": "IL7R, CCR7, TCF7, and LEF1 support central memory identity."
                    },
                    {
                        "cluster_id": "effector",
                        "main_cell_type": "Effector CD8 T cell",
                        "sub_cell_type": "Cytotoxic effector CD8 T cell",
                        "key_markers": ["GZMB", "PRF1", "NKG7", "CCL5"],
                        "reason": "GZMB, PRF1, NKG7, and CCL5 support cytotoxic effector identity."
                    }
                ]
            })
        elif "Cluster ID: B_cell" in prompt or "MS4A1" in prompt:
            print_json({
                "main_cell_type": "B cell",
                "sub_cell_types": ["Naive B cell", "MS4A1-positive B cell", "antigen-presenting B cell"],
                "possible_mixed_cell_types": [],
                "key_functional_markers": [{"genes": ["CD74", "HLA-DRA"], "interpretation": "antigen presentation"}],
                "key_cell_type_markers": [{"genes": ["MS4A1", "CD79A", "CD79B"], "supports": "B cell", "interpretation": "canonical B cell receptor markers"}],
                "confidence": "high",
                "evidence": "MS4A1, CD79A, CD79B, CD74, and BANK1 support a B cell annotation."
            })
        elif "Cluster ID: Macrophage" in prompt or "C1QA" in prompt:
            print_json({
                "main_cell_type": "Macrophage",
                "sub_cell_types": ["C1Q/APOE macrophage", "tissue macrophage", "monocyte-derived macrophage"],
                "possible_mixed_cell_types": [],
                "key_functional_markers": [{"genes": ["C1QA", "C1QB", "APOE"], "interpretation": "complement and lipid handling"}],
                "key_cell_type_markers": [{"genes": ["LYZ", "LST1", "MS4A7"], "supports": "macrophage", "interpretation": "myeloid/macrophage markers"}],
                "confidence": "high",
                "evidence": "LYZ, LST1, C1QA, C1QB, APOE, and MS4A7 support macrophage identity."
            })
        else:
            print_json({
                "main_cell_type": "T cell",
                "sub_cell_types": ["Central memory T cell", "Naive T cell", "CD4/CD8 unresolved T cell"],
                "possible_mixed_cell_types": [],
                "key_functional_markers": [{"genes": ["IL7R", "CCR7", "TCF7"], "interpretation": "memory/naive T cell program"}],
                "key_cell_type_markers": [{"genes": ["CD3D", "CD3E", "TRAC"], "supports": "T cell", "interpretation": "canonical T cell receptor complex markers"}],
                "confidence": "high",
                "evidence": "CD3D, CD3E, TRAC, IL7R, CCR7, and TCF7 support a T cell annotation."
            })
        """
    )


def _readme_text(backend: str) -> str:
    return textwrap.dedent(
        f"""\
        # CASSIA CLI Example Project

        This folder was generated by `cassia examples`.

        Files:
        - `markers.csv`: preformatted cluster marker lists for `cassia annotate`
        - `raw_markers.csv`: long differential-expression marker table for `cassia validate` and `cassia boost query`
        - `subcluster_markers.csv`: marker lists for `cassia subcluster run`
        - `consensus_inputs/*/summary.csv`: small summary files for `cassia consensus`
        - `toy_agent.py`: deterministic offline shell backend for smoke tests
        - `run_offline.sh`: no-API smoke test using `toy_agent.py`
        - `run.sh`: same workflow using the selected backend

        Quick start:

        ```bash
        bash run_offline.sh
        ```

        Real agent run:

        ```bash
        bash run.sh {backend}
        ```

        If `cassia` is not on PATH, run scripts with:

        ```bash
        CASSIA_CMD="python -m CASSIA.cli" bash run_offline.sh
        ```
        """
    )


def _run_script(backend: str) -> str:
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env bash
        set -euo pipefail

        ROOT="$(pwd)"
        PYTHON_BIN="${{PYTHON_BIN:-python}}"
        CASSIA_CMD="${{CASSIA_CMD:-$PYTHON_BIN -m CASSIA.cli}}"
        BACKEND="${{1:-{backend}}}"
        TOY_AGENT="$ROOT/toy_agent.py"
        AGENT_ARGS=(--backend "$BACKEND")
        if [[ "$BACKEND" == "shell" ]]; then
          AGENT_ARGS=(--backend shell --command-template "$PYTHON_BIN $TOY_AGENT {{prompt_file}}")
        fi

        echo "== Validate marker inputs =="
        $CASSIA_CMD validate "$ROOT/markers.csv"
        $CASSIA_CMD validate "$ROOT/raw_markers.csv" --celltype-column cluster --gene-column gene --n-genes 5

        echo "== Annotate clusters with $BACKEND =="
        $CASSIA_CMD annotate \\
          --input "$ROOT/markers.csv" \\
          "${{AGENT_ARGS[@]}}" \\
          --tissue blood \\
          --species human \\
          --out "$ROOT/runs/agent_annotation"

        echo "== Query local marker evidence =="
        $CASSIA_CMD boost query \\
          --markers "$ROOT/raw_markers.csv" \\
          --cluster T_cell \\
          --genes CD3D,TRAC,NKG7

        echo "== Plan annotation boost candidates =="
        $CASSIA_CMD boost auto \\
          --run "$ROOT/runs/agent_annotation" \\
          --markers "$ROOT/raw_markers.csv" \\
          --backend "$BACKEND" \\
          --max-clusters 2 \\
          --plan-only

        echo "== Dry-run subcluster prompt =="
        $CASSIA_CMD subcluster run \\
          --markers "$ROOT/subcluster_markers.csv" \\
          --major-cluster-info "CD8 T cell in blood" \\
          "${{AGENT_ARGS[@]}}" \\
          --dry-run \\
          --out "$ROOT/runs/subcluster_dry_run"

        echo "== Build consensus from example summaries =="
        $CASSIA_CMD consensus \\
          --inputs "$ROOT/consensus_inputs/codex" "$ROOT/consensus_inputs/claude" \\
          --out "$ROOT/runs/example_consensus.csv"
        """
    )


def _offline_script() -> str:
    return textwrap.dedent(
        """\
        #!/usr/bin/env bash
        set -euo pipefail

        ROOT="$(pwd)"
        PYTHON_BIN="${PYTHON_BIN:-python}"
        CASSIA_CMD="${CASSIA_CMD:-$PYTHON_BIN -m CASSIA.cli}"
        TOY_AGENT="$ROOT/toy_agent.py"
        SHELL_AGENT="$PYTHON_BIN $TOY_AGENT {prompt_file}"

        echo "== Validate marker inputs =="
        $CASSIA_CMD validate "$ROOT/markers.csv"
        $CASSIA_CMD validate "$ROOT/raw_markers.csv" --celltype-column cluster --gene-column gene --n-genes 5

        echo "== Offline annotation smoke test =="
        $CASSIA_CMD annotate \\
          --input "$ROOT/markers.csv" \\
          --backend shell \\
          --command-template "$SHELL_AGENT" \\
          --tissue blood \\
          --species human \\
          --out "$ROOT/runs/offline_annotation"

        echo "== Query local marker evidence =="
        $CASSIA_CMD boost query \\
          --markers "$ROOT/raw_markers.csv" \\
          --cluster T_cell \\
          --genes CD3D,TRAC,NKG7

        echo "== Offline boost planning =="
        $CASSIA_CMD boost auto \\
          --run "$ROOT/runs/offline_annotation" \\
          --markers "$ROOT/raw_markers.csv" \\
          --backend shell \\
          --plan-only \\
          --max-clusters 2

        echo "== Offline subcluster annotation =="
        $CASSIA_CMD subcluster run \\
          --markers "$ROOT/subcluster_markers.csv" \\
          --major-cluster-info "CD8 T cell in blood" \\
          --backend shell \\
          --command-template "$SHELL_AGENT" \\
          --out "$ROOT/runs/offline_subcluster"

        echo "== Consensus from example summaries =="
        $CASSIA_CMD consensus \\
          --inputs "$ROOT/consensus_inputs/codex" "$ROOT/consensus_inputs/claude" \\
          --out "$ROOT/runs/example_consensus.csv"
        """
    )


def create_example_project(out_dir: Path, backend: str, force: bool = False) -> List[Path]:
    """Create a runnable CASSIA CLI example project."""
    out_dir = out_dir.expanduser()
    if out_dir.exists() and any(out_dir.iterdir()) and not force:
        raise FileExistsError(f"Output directory already exists and is not empty: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    written: List[Path] = []
    markers_path = out_dir / "markers.csv"
    _write_dict_csv(markers_path, ["cluster", "markers"], MARKER_ROWS)
    written.append(markers_path)

    raw_path = out_dir / "raw_markers.csv"
    _write_csv(raw_path, ["cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"], RAW_MARKER_ROWS)
    written.append(raw_path)

    subcluster_path = out_dir / "subcluster_markers.csv"
    _write_dict_csv(subcluster_path, ["cluster", "markers"], SUBCLUSTER_ROWS)
    written.append(subcluster_path)

    for source, rows in SUMMARY_ROWS_BY_SOURCE.items():
        summary_path = out_dir / "consensus_inputs" / source / "summary.csv"
        _write_csv(summary_path, SUMMARY_HEADERS, rows)
        written.append(summary_path)

    files = {
        "README.md": _readme_text(backend),
        "toy_agent.py": _toy_agent_source(),
        "run.sh": _run_script(backend),
        "run_offline.sh": _offline_script(),
    }
    for filename, content in files.items():
        path = out_dir / filename
        path.write_text(content, encoding="utf-8")
        written.append(path)
        if filename.endswith(".sh") or filename == "toy_agent.py":
            path.chmod(0o755)

    (out_dir / "runs").mkdir(exist_ok=True)
    return written


def run_examples(args: Any) -> int:
    """Generate an example project from CLI arguments."""
    out_dir = Path(args.out)
    try:
        written = create_example_project(out_dir=out_dir, backend=args.backend, force=args.force)
    except FileExistsError as exc:
        print(str(exc))
        print("Use --force to write example files into this directory.")
        return 1

    print(f"Created CASSIA example project: {out_dir}")
    print(f"Wrote {len(written)} file(s).")
    print("")
    print("Try:")
    print(f"  cd {out_dir}")
    print("  bash run_offline.sh")
    print(f"  bash run.sh {args.backend}")
    return 0
