"""Command-line interface for CASSIA."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

from CASSIA import __version__

from .backends import is_agent_backend, is_api_backend, list_backends
from .boost import parse_gene_args, query_marker_genes, run_boost, run_boost_auto, write_query_output
from .consensus import run_consensus
from .runner import dispatch_annotation, generate_markdown_report, resume_run
from .subcluster import run_subcluster
from .validate import run_validate


HELP_FORMATTER = argparse.RawDescriptionHelpFormatter

TOP_LEVEL_EPILOG = """Common workflows:
  cassia doctor
  cassia backends list
  cassia validate markers.csv

  cassia annotate -i markers.csv --backend codex-cli --tissue brain --species human --out runs/brain_codex
  cassia boost auto --run runs/brain_codex --markers raw_markers.csv --backend codex-cli
  cassia subcluster run --markers cd8_markers.csv --major-cluster-info "CD8 T cell in tumor" --backend codex-cli
  cassia consensus --inputs runs/codex/summary.csv runs/claude/summary.csv --out runs/consensus.csv

Run "cassia COMMAND --help" for command-specific options and examples.
"""

VALIDATE_EPILOG = """Examples:
  cassia validate markers.csv
  cassia validate raw_findallmarkers.csv --celltype-column cluster --gene-column gene --ranking-method avg_log2FC
  cassia validate markers.csv --json
"""

ANNOTATE_EPILOG = """Examples:
  cassia annotate -i markers.csv --backend codex-cli --tissue brain --species human --out runs/brain_codex

  cassia annotate -i findallmarkers.csv --celltype-column cluster --gene-column gene \\
    --ranking-method avg_log2FC --n-genes 50 --backend claude-cli --out runs/claude

  cassia annotate -i markers.csv --backend shell \\
    --command-template 'python my_agent.py {prompt_file}' --out runs/custom_agent
"""

BOOST_EPILOG = """Examples:
  cassia boost query --markers raw_findallmarkers.csv --cluster 3 --genes CD3D,CD3E,TRAC
  cassia boost run --run runs/brain_codex --markers raw_findallmarkers.csv --cluster 3 --backend codex-cli
  cassia boost auto --run runs/brain_codex --markers raw_findallmarkers.csv --backend codex-cli --max-clusters 5
"""

BOOST_QUERY_EPILOG = """Examples:
  cassia boost query --markers raw_findallmarkers.csv --cluster 3 --genes CD3D,CD3E,TRAC
  cassia boost query -m raw_findallmarkers.csv --cluster macrophage LST1 C1QA APOE --format json
"""

BOOST_RUN_EPILOG = """Example:
  cassia boost run --run runs/brain_codex --markers raw_findallmarkers.csv \\
    --cluster 3 --backend codex-cli --iterations 5
"""

BOOST_AUTO_EPILOG = """Examples:
  cassia boost auto --run runs/brain_codex --markers raw_findallmarkers.csv --backend codex-cli --max-clusters 5
  cassia boost auto --run runs/brain_codex --markers raw_findallmarkers.csv --plan-only
"""

SUBCLUSTER_EPILOG = """Example:
  cassia subcluster run --markers cd8_subcluster_markers.csv \\
    --major-cluster-info "CD8 T cell in human tumor" --backend codex-cli --out runs/cd8_subcluster
"""

CONSENSUS_EPILOG = """Examples:
  cassia consensus --inputs runs/codex/summary.csv runs/claude/summary.csv --out runs/consensus.csv
  cassia consensus --inputs runs/codex runs/claude runs/cursor --out runs/consensus.csv
  cassia consensus --glob 'runs/*/summary.csv' --threshold 0.75 --out runs/consensus.csv
"""


def _print_backend_table(backends: Dict[str, Dict[str, object]]) -> None:
    print("Backend        Kind       Available  Details")
    print("-------------  ---------  ---------  -------")
    for name in sorted(backends):
        info = backends[name]
        available = "yes" if info.get("available") else "no"
        if info.get("kind") == "api":
            detail = f"env: {info.get('requires')}"
        elif info.get("executable"):
            detail = f"exe: {info.get('executable')}"
            if info.get("path"):
                detail += f" ({info.get('path')})"
        else:
            detail = str(info.get("description", ""))
        print(f"{name:13}  {str(info.get('kind')):9}  {available:9}  {detail}")


def cmd_backends_list(args: argparse.Namespace) -> int:
    backends = list_backends()
    if args.json:
        print(json.dumps(backends, indent=2))
    else:
        _print_backend_table(backends)
    return 0


def cmd_doctor(args: argparse.Namespace) -> int:
    backends = list_backends()
    checks = {
        "cassia_version": __version__,
        "python": sys.version.split()[0],
        "api_env": {
            name: bool(os.environ.get(str(info.get("requires"))))
            for name, info in backends.items()
            if info.get("kind") == "api"
        },
        "agent_cli": {
            name: {
                "available": bool(info.get("available")),
                "executable": info.get("executable"),
                "path": info.get("path"),
            }
            for name, info in backends.items()
            if info.get("kind") == "agent-cli" and name != "shell"
        },
    }

    if args.json:
        print(json.dumps(checks, indent=2))
    else:
        print(f"CASSIA version: {checks['cassia_version']}")
        print(f"Python: {checks['python']}")
        print("")
        print("API keys:")
        for name, available in checks["api_env"].items():
            print(f"  {name}: {'set' if available else 'missing'}")
        print("")
        print("Agent CLIs:")
        for name, info in checks["agent_cli"].items():
            status = "found" if info["available"] else "missing"
            path = f" ({info['path']})" if info.get("path") else ""
            print(f"  {name}: {status}{path}")

    if args.strict:
        any_backend = any(
            info.get("available")
            for name, info in backends.items()
            if name != "shell"
        )
        return 0 if any_backend else 1
    return 0


def cmd_init(args: argparse.Namespace) -> int:
    config_dir = Path(args.directory) / ".cassia"
    config_dir.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "config.json"
    if config_path.exists() and not args.force:
        print(f"CASSIA config already exists: {config_path}")
        return 0

    config = {
        "default_backend": args.backend,
        "tissue": args.tissue,
        "species": args.species,
        "created_by": "cassia init",
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"Created {config_path}")
    return 0


def cmd_annotate(args: argparse.Namespace) -> int:
    if not is_api_backend(args.backend) and not is_agent_backend(args.backend):
        raise SystemExit(
            f"Unknown backend '{args.backend}'. Run 'cassia backends list' to see supported backends."
        )
    if args.backend == "shell" and not args.command_template:
        raise SystemExit("--command-template is required when --backend shell is used")
    return dispatch_annotation(args)


def cmd_report(args: argparse.Namespace) -> int:
    report_path = generate_markdown_report(Path(args.run_dir))
    print(f"Wrote {report_path}")
    return 0


def cmd_resume(args: argparse.Namespace) -> int:
    return resume_run(Path(args.run_dir))


def cmd_boost_query(args: argparse.Namespace) -> int:
    genes = parse_gene_args((args.genes or []) + (args.gene or []))
    if not genes:
        raise SystemExit("Provide at least one gene with --genes or positional gene arguments")

    result = query_marker_genes(
        marker_path=Path(args.markers),
        genes=genes,
        cluster=args.cluster,
        gene_column=args.gene_column,
        cluster_column=args.cluster_column,
    )
    text = write_query_output(result, output_format=args.format, out=Path(args.out) if args.out else None)
    if args.out:
        print(f"Wrote {args.out}")
    else:
        print(text)
    return 0


def cmd_boost_run(args: argparse.Namespace) -> int:
    if not is_agent_backend(args.backend):
        raise SystemExit("cassia boost run requires an agent CLI backend: claude-cli, codex-cli, cursor-agent, or shell")
    if args.backend == "shell" and not args.command_template:
        raise SystemExit("--command-template is required when --backend shell is used")
    return run_boost(args)


def cmd_boost_auto(args: argparse.Namespace) -> int:
    if not is_agent_backend(args.backend):
        raise SystemExit("cassia boost auto requires an agent CLI backend: claude-cli, codex-cli, cursor-agent, or shell")
    if args.backend == "shell" and not args.command_template and not (args.plan_only or args.dry_run):
        raise SystemExit("--command-template is required when --backend shell is used")
    return run_boost_auto(args)


def cmd_subcluster_run(args: argparse.Namespace) -> int:
    if not is_agent_backend(args.backend):
        raise SystemExit("cassia subcluster run requires an agent CLI backend: claude-cli, codex-cli, cursor-agent, or shell")
    if args.backend == "shell" and not args.command_template and not args.dry_run:
        raise SystemExit("--command-template is required when --backend shell is used")
    return run_subcluster(args)


def cmd_consensus(args: argparse.Namespace) -> int:
    return run_consensus(args)


def cmd_validate(args: argparse.Namespace) -> int:
    return run_validate(args)


def add_common_annotation_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-i", "--input", required=True, help="Input marker CSV file.")
    parser.add_argument("-o", "--out", help="Run output directory.")
    parser.add_argument(
        "--backend",
        default="openrouter",
        help=(
            "Backend/provider to use: openrouter, openai, anthropic, claude-cli, "
            "codex-cli, cursor-agent, shell, or an OpenAI-compatible base URL."
        ),
    )
    parser.add_argument("--model", help="Model name for API-backed CASSIA runs.")
    parser.add_argument("--temperature", type=float, help="Temperature for API-backed CASSIA runs.")
    parser.add_argument("--tissue", default="lung", help="Tissue context.")
    parser.add_argument("--species", default="human", help="Species context.")
    parser.add_argument("--additional-info", help="Extra dataset context for annotation.")
    parser.add_argument("--celltype-column", help="Cluster/cell type column name.")
    parser.add_argument("--gene-column", help="Gene or marker-list column name.")
    parser.add_argument("--n-genes", type=int, default=50, help="Top marker count per cluster.")
    parser.add_argument("--max-workers", type=int, default=10, help="Parallel workers for API backend.")
    parser.add_argument("--max-retries", type=int, default=1, help="Retries for API backend.")
    parser.add_argument(
        "--ranking-method",
        default="avg_log2FC",
        choices=["avg_log2FC", "p_val_adj", "pct_diff", "Score"],
        help="Marker ranking method for differential-expression tables.",
    )
    sort_group = parser.add_mutually_exclusive_group()
    sort_group.add_argument("--ascending", action="store_true", dest="ascending")
    sort_group.add_argument("--descending", action="store_false", dest="ascending")
    parser.set_defaults(ascending=None)
    parser.add_argument("--validator-involvement", default="v1", help="CASSIA validator mode for API backend.")
    parser.add_argument("--reasoning", help="Reasoning effort for API models that support it.")
    parser.add_argument("--use-reference", action="store_true", help="Use CASSIA reference retrieval in API backend.")
    parser.add_argument("--reference-model", help="Model used by reference retrieval.")
    parser.add_argument("--reference-cell-type-hint", help="Parent lineage hint for reference retrieval.")
    parser.add_argument("--skip-api-key-validation", action="store_true", help="Skip API key preflight validation.")
    parser.add_argument(
        "--command-template",
        "--shell-command",
        dest="command_template",
        help=(
            "Shell command template for --backend shell. Placeholders: "
            "{prompt}, {prompt_file}, {input}, {out}, {cluster}."
        ),
    )
    parser.add_argument("--timeout", type=int, default=900, help="Agent CLI timeout per cluster in seconds.")
    parser.add_argument("--dry-run", action="store_true", help="Write prompts and manifest without calling a backend.")
    parser.add_argument("--resume", action="store_true", help="Skip clusters already present in results.json.")
    parser.add_argument("--keep-going", action="store_true", help="Continue local agent runs after a cluster fails.")
    parser.add_argument("--limit", type=int, help="Limit clusters, useful for smoke tests.")
    parser.add_argument("--quiet", action="store_true", help="Reduce API backend progress output.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cassia",
        description=(
            "CASSIA command-line interface for single-cell marker annotation, "
            "agent-native boost review, subcluster annotation, and consensus voting."
        ),
        formatter_class=HELP_FORMATTER,
        epilog=TOP_LEVEL_EPILOG,
    )
    parser.add_argument("--version", action="version", version=f"CASSIA {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Create a local .cassia config file.")
    init_parser.add_argument("--directory", default=".", help="Directory where .cassia/config.json is created.")
    init_parser.add_argument("--backend", default="openrouter", help="Default backend to write into config.")
    init_parser.add_argument("--tissue", default="lung", help="Default tissue context.")
    init_parser.add_argument("--species", default="human", help="Default species context.")
    init_parser.add_argument("--force", action="store_true", help="Overwrite an existing config.")
    init_parser.set_defaults(func=cmd_init)

    validate_parser = subparsers.add_parser(
        "validate",
        help="Validate marker CSV input before running annotation.",
        description="Validate marker CSV structure, inferred columns, ranking columns, and prepared marker counts.",
        formatter_class=HELP_FORMATTER,
        epilog=VALIDATE_EPILOG,
    )
    validate_parser.add_argument("input", help="Input marker CSV file.")
    validate_parser.add_argument("--celltype-column", help="Cluster/cell type column name.")
    validate_parser.add_argument("--gene-column", help="Gene or marker-list column name.")
    validate_parser.add_argument("--n-genes", type=int, default=50, help="Top marker count per cluster.")
    validate_parser.add_argument(
        "--ranking-method",
        default="avg_log2FC",
        choices=["avg_log2FC", "p_val_adj", "pct_diff", "Score"],
        help="Marker ranking method for differential-expression tables.",
    )
    validate_sort_group = validate_parser.add_mutually_exclusive_group()
    validate_sort_group.add_argument("--ascending", action="store_true", dest="ascending")
    validate_sort_group.add_argument("--descending", action="store_false", dest="ascending")
    validate_parser.set_defaults(ascending=None)
    validate_parser.add_argument("--limit", type=int, help="Limit clusters checked after preparation.")
    validate_parser.add_argument("--tissue", help="Tissue context used only for the suggested command.")
    validate_parser.add_argument("--species", default="human", help="Species used only for the suggested command.")
    validate_parser.add_argument("--backend", default="codex-cli", help="Backend used only for the suggested command.")
    validate_parser.add_argument("--json", action="store_true", help="Print machine-readable validation diagnostics.")
    validate_parser.add_argument("--strict", action="store_true", help="Exit non-zero when warnings are present.")
    validate_parser.set_defaults(func=cmd_validate)

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Run CASSIA annotation from marker CSV data.",
        description="Annotate clusters from marker CSV data using an API backend or a local agent CLI.",
        formatter_class=HELP_FORMATTER,
        epilog=ANNOTATE_EPILOG,
    )
    add_common_annotation_args(annotate_parser)
    annotate_parser.set_defaults(func=cmd_annotate)

    backends_parser = subparsers.add_parser("backends", help="Inspect available CASSIA backends.")
    backends_subparsers = backends_parser.add_subparsers(dest="backends_command")
    list_parser = backends_subparsers.add_parser("list", help="List known backends.")
    list_parser.add_argument("--json", action="store_true", help="Print backend metadata as JSON.")
    list_parser.set_defaults(func=cmd_backends_list)
    backends_parser.set_defaults(func=cmd_backends_list, json=False)

    doctor_parser = subparsers.add_parser("doctor", help="Check CASSIA CLI environment readiness.")
    doctor_parser.add_argument("--json", action="store_true", help="Print diagnostics as JSON.")
    doctor_parser.add_argument("--strict", action="store_true", help="Exit non-zero if no backend is available.")
    doctor_parser.set_defaults(func=cmd_doctor)

    report_parser = subparsers.add_parser("report", help="Generate a Markdown report for a run folder.")
    report_parser.add_argument("run_dir", help="Run directory containing results.json.")
    report_parser.set_defaults(func=cmd_report)

    resume_parser = subparsers.add_parser("resume", help="Resume an agent CLI run folder.")
    resume_parser.add_argument("run_dir", help="Run directory containing run_manifest.json.")
    resume_parser.set_defaults(func=cmd_resume)

    boost_parser = subparsers.add_parser(
        "boost",
        help="Annotation boost helper commands.",
        description="Inspect marker evidence and re-check uncertain CASSIA annotations.",
        formatter_class=HELP_FORMATTER,
        epilog=BOOST_EPILOG,
    )
    boost_subparsers = boost_parser.add_subparsers(dest="boost_command", required=True)
    query_parser = boost_subparsers.add_parser(
        "query",
        help="Query marker statistics for genes, optionally within one cluster.",
        description="Query local marker statistics for genes, optionally within one cluster.",
        formatter_class=HELP_FORMATTER,
        epilog=BOOST_QUERY_EPILOG,
    )
    query_parser.add_argument("-m", "--markers", required=True, help="Raw marker table CSV.")
    query_parser.add_argument("-g", "--genes", action="append", help="Comma/space separated genes to query.")
    query_parser.add_argument("gene", nargs="*", help="Additional genes to query.")
    query_parser.add_argument("--cluster", help="Cluster ID/name to filter before querying.")
    query_parser.add_argument("--gene-column", help="Gene column name. Defaults to auto-detection.")
    query_parser.add_argument("--cluster-column", help="Cluster column name. Defaults to auto-detection.")
    query_parser.add_argument("--format", choices=["table", "csv", "json"], default="table", help="Output format.")
    query_parser.add_argument("-o", "--out", help="Optional output file.")
    query_parser.set_defaults(func=cmd_boost_query)

    run_parser = boost_subparsers.add_parser(
        "run",
        help="Run an agent CLI annotation boost loop for one cluster.",
        description="Run the annotation boost loop for one selected cluster using a local agent CLI.",
        formatter_class=HELP_FORMATTER,
        epilog=BOOST_RUN_EPILOG,
    )
    run_parser.add_argument("--run", required=True, help="Existing CASSIA run directory.")
    run_parser.add_argument("-m", "--markers", required=True, help="Raw marker table CSV.")
    run_parser.add_argument("--cluster", required=True, help="Cluster ID/name to boost.")
    run_parser.add_argument("--major-cluster-info", default="single-cell RNA-seq dataset", help="Dataset context for the boost prompt.")
    run_parser.add_argument("--backend", default="codex-cli", help="Agent backend: codex-cli, claude-cli, cursor-agent, or shell.")
    run_parser.add_argument(
        "--command-template",
        "--shell-command",
        dest="command_template",
        help="Shell command template for --backend shell. Placeholders: {prompt}, {prompt_file}, {input}, {out}, {cluster}.",
    )
    run_parser.add_argument("--strategy", choices=["breadth", "depth"], default="breadth", help="Boost search strategy.")
    run_parser.add_argument("--iterations", type=int, default=5, help="Maximum agent rounds before finalization.")
    run_parser.add_argument("--n-genes", type=int, default=50, help="Top raw markers to include in the initial prompt.")
    run_parser.add_argument("--max-genes-per-round", type=int, default=20, help="Maximum genes accepted from each <check_genes> request.")
    run_parser.add_argument("--gene-column", help="Gene column name. Defaults to auto-detection.")
    run_parser.add_argument("--cluster-column", help="Cluster column name. Defaults to auto-detection.")
    run_parser.add_argument(
        "--ranking-method",
        default="avg_log2FC",
        choices=["avg_log2FC", "p_val_adj", "p_val", "pct.1", "pct.2", "Score"],
        help="Marker ranking column for selecting top markers.",
    )
    sort_group = run_parser.add_mutually_exclusive_group()
    sort_group.add_argument("--ascending", action="store_true", dest="ascending")
    sort_group.add_argument("--descending", action="store_false", dest="ascending")
    run_parser.set_defaults(ascending=None)
    run_parser.add_argument("--additional-task", help="Optional extra question, e.g. check if this is malignant.")
    run_parser.add_argument("--timeout", type=int, default=900, help="Agent CLI timeout per round in seconds.")
    run_parser.add_argument("--dry-run", action="store_true", help="Write the initial prompt and manifest without calling an agent.")
    run_parser.add_argument("-o", "--out", help="Optional boost output directory.")
    run_parser.set_defaults(func=cmd_boost_run)

    auto_parser = boost_subparsers.add_parser(
        "auto",
        help="Auto-select uncertain clusters from a run folder and boost them.",
        description="Select low-confidence, mixed, or ambiguous clusters from a run folder and boost them.",
        formatter_class=HELP_FORMATTER,
        epilog=BOOST_AUTO_EPILOG,
    )
    auto_parser.add_argument("--run", required=True, help="Existing CASSIA run directory.")
    auto_parser.add_argument("-m", "--markers", required=True, help="Raw marker table CSV.")
    auto_parser.add_argument("--major-cluster-info", default="single-cell RNA-seq dataset", help="Dataset context for the boost prompt.")
    auto_parser.add_argument("--backend", default="codex-cli", help="Agent backend: codex-cli, claude-cli, cursor-agent, or shell.")
    auto_parser.add_argument(
        "--command-template",
        "--shell-command",
        dest="command_template",
        help="Shell command template for --backend shell. Placeholders: {prompt}, {prompt_file}, {input}, {out}, {cluster}.",
    )
    auto_parser.add_argument("--strategy", choices=["breadth", "depth"], default="breadth", help="Boost search strategy.")
    auto_parser.add_argument("--iterations", type=int, default=5, help="Maximum agent rounds before finalization.")
    auto_parser.add_argument("--n-genes", type=int, default=50, help="Top raw markers to include in each boost prompt.")
    auto_parser.add_argument("--max-genes-per-round", type=int, default=20, help="Maximum genes accepted from each <check_genes> request.")
    auto_parser.add_argument("--gene-column", help="Gene column name. Defaults to auto-detection.")
    auto_parser.add_argument("--cluster-column", help="Cluster column name. Defaults to auto-detection.")
    auto_parser.add_argument(
        "--ranking-method",
        default="avg_log2FC",
        choices=["avg_log2FC", "p_val_adj", "p_val", "pct.1", "pct.2", "Score"],
        help="Marker ranking column for selecting top markers.",
    )
    auto_sort_group = auto_parser.add_mutually_exclusive_group()
    auto_sort_group.add_argument("--ascending", action="store_true", dest="ascending")
    auto_sort_group.add_argument("--descending", action="store_false", dest="ascending")
    auto_parser.set_defaults(ascending=None)
    auto_parser.add_argument("--additional-task", help="Optional extra instruction appended to each auto-selected boost prompt.")
    auto_parser.add_argument("--timeout", type=int, default=900, help="Agent CLI timeout per round in seconds.")
    auto_parser.add_argument("--dry-run", action="store_true", help="Write boost prompts and manifests without calling an agent.")
    auto_parser.add_argument("-o", "--out", help="Optional auto output directory. Defaults to RUN/boost/_auto.")
    auto_parser.add_argument("--max-clusters", type=int, default=5, help="Maximum clusters to boost.")
    auto_parser.add_argument(
        "--confidence",
        action="append",
        choices=["low", "medium", "high", "unknown"],
        help="Confidence level to prioritize. Repeatable. Defaults to low, medium, and unknown.",
    )
    auto_parser.add_argument("--only-low-confidence", action="store_true", help="Only select clusters with low confidence.")
    auto_parser.add_argument("--target-lineage", action="append", help="Comma-separated lineage terms to focus on, e.g. macrophage,T cell.")
    auto_parser.add_argument("--all", action="store_true", help="Allow all clusters to be selected, ordered after risk-scored clusters.")
    auto_parser.add_argument("--plan-only", action="store_true", help="Write candidate plan and aggregate reports without running boost.")
    auto_parser.add_argument("--force", action="store_true", help="Re-run clusters even if final.json already exists.")
    auto_parser.add_argument("--fail-fast", action="store_true", help="Stop after the first failed boost run.")
    auto_parser.set_defaults(func=cmd_boost_auto)

    subcluster_parser = subparsers.add_parser(
        "subcluster",
        help="Agent-native subcluster annotation commands.",
        description="Annotate subclusters inside one parent population using a local agent CLI.",
        formatter_class=HELP_FORMATTER,
        epilog=SUBCLUSTER_EPILOG,
    )
    subcluster_subparsers = subcluster_parser.add_subparsers(dest="subcluster_command", required=True)
    subcluster_run_parser = subcluster_subparsers.add_parser(
        "run",
        help="Annotate parent-cluster subclusters from a marker table using an agent CLI.",
        description="Annotate parent-cluster subclusters from a marker table using an agent CLI.",
        formatter_class=HELP_FORMATTER,
        epilog=SUBCLUSTER_EPILOG,
    )
    subcluster_run_parser.add_argument("-m", "--markers", required=True, help="Subcluster marker table CSV.")
    subcluster_run_parser.add_argument(
        "--major-cluster-info",
        required=True,
        help="Parent cluster context, e.g. 'CD8 T cell in human tumor'.",
    )
    subcluster_run_parser.add_argument("--backend", default="codex-cli", help="Agent backend: codex-cli, claude-cli, cursor-agent, or shell.")
    subcluster_run_parser.add_argument(
        "--command-template",
        "--shell-command",
        dest="command_template",
        help="Shell command template for --backend shell. Placeholders: {prompt}, {prompt_file}, {input}, {out}, {cluster}.",
    )
    subcluster_run_parser.add_argument("-o", "--out", help="Output directory. Defaults to cassia_runs/subcluster_TIMESTAMP.")
    subcluster_run_parser.add_argument("--tissue", help="Tissue context.")
    subcluster_run_parser.add_argument("--species", help="Species context.")
    subcluster_run_parser.add_argument("--additional-context", help="Optional context appended to the subcluster prompt.")
    subcluster_run_parser.add_argument("--n-genes", type=int, default=50, help="Top marker count per subcluster.")
    subcluster_run_parser.add_argument("--cluster-column", help="Subcluster ID column name. Defaults to auto-detection.")
    subcluster_run_parser.add_argument("--gene-column", help="Gene or marker-list column name. Defaults to auto-detection.")
    subcluster_run_parser.add_argument(
        "--ranking-method",
        default="avg_log2FC",
        choices=["avg_log2FC", "p_val_adj", "pct_diff", "Score"],
        help="Marker ranking method for differential-expression tables.",
    )
    subcluster_sort_group = subcluster_run_parser.add_mutually_exclusive_group()
    subcluster_sort_group.add_argument("--ascending", action="store_true", dest="ascending")
    subcluster_sort_group.add_argument("--descending", action="store_false", dest="ascending")
    subcluster_run_parser.set_defaults(ascending=None)
    subcluster_run_parser.add_argument("--limit", type=int, help="Limit subclusters, useful for smoke tests.")
    subcluster_run_parser.add_argument("--timeout", type=int, default=900, help="Agent CLI timeout in seconds.")
    subcluster_run_parser.add_argument("--dry-run", action="store_true", help="Write prompt and manifest without calling an agent.")
    subcluster_run_parser.add_argument("--use-reference", action="store_true", help="Retrieve CASSIA subtype references before prompting.")
    subcluster_run_parser.add_argument("--reference-provider", help="Provider for reference retrieval. Defaults to openrouter.")
    subcluster_run_parser.add_argument("--reference-model", help="Model used by reference retrieval.")
    subcluster_run_parser.add_argument("--reference-cell-type-hint", help="Optional parent-lineage hint for reference retrieval.")
    subcluster_run_parser.add_argument("--reference-depth", choices=["summary", "detailed"], default="detailed", help="Reference extraction depth.")
    subcluster_run_parser.add_argument("--reference-max-content-length", type=int, default=5000, help="Maximum characters per reference result.")
    subcluster_run_parser.add_argument("--reference-max-context-length", type=int, default=12000, help="Maximum total reference context characters.")
    subcluster_run_parser.set_defaults(func=cmd_subcluster_run)

    consensus_parser = subparsers.add_parser(
        "consensus",
        help="Build deterministic consensus from multiple CASSIA summary CSVs.",
        description=(
            "Build deterministic consensus from multiple CASSIA summary, "
            "subcluster, or boost-auto CSV outputs without calling an LLM."
        ),
        formatter_class=HELP_FORMATTER,
        epilog=CONSENSUS_EPILOG,
    )
    consensus_parser.add_argument(
        "--inputs",
        nargs="+",
        help="Input CASSIA summary/subcluster/boost-auto CSV files or run directories.",
    )
    consensus_parser.add_argument(
        "--glob",
        dest="glob_patterns",
        action="append",
        help="Glob pattern for input CSV files or run directories. Repeatable.",
    )
    consensus_parser.add_argument("-o", "--out", default="consensus.csv", help="Output consensus CSV path.")
    consensus_parser.add_argument(
        "--html",
        help="Output HTML report path. Defaults to the CSV output path with .html suffix.",
    )
    consensus_parser.add_argument("--no-html", action="store_true", help="Do not write the HTML consensus report.")
    consensus_parser.add_argument(
        "--threshold",
        type=float,
        default=2 / 3,
        help="Minimum top-vote fraction for consensus status. Defaults to two-thirds.",
    )
    consensus_parser.add_argument("--cluster-column", help="Cluster ID column override.")
    consensus_parser.add_argument("--main-column", help="Broad/main annotation column override.")
    consensus_parser.add_argument("--sub-column", help="Subtype/detail annotation column override.")
    consensus_parser.set_defaults(func=cmd_consensus)

    return parser


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
