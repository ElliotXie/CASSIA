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
from .boost import parse_gene_args, query_marker_genes, run_boost, write_query_output
from .runner import dispatch_annotation, generate_markdown_report, resume_run


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
    parser = argparse.ArgumentParser(prog="cassia", description="CASSIA command-line interface.")
    parser.add_argument("--version", action="version", version=f"CASSIA {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Create a local .cassia config file.")
    init_parser.add_argument("--directory", default=".", help="Directory where .cassia/config.json is created.")
    init_parser.add_argument("--backend", default="openrouter", help="Default backend to write into config.")
    init_parser.add_argument("--tissue", default="lung", help="Default tissue context.")
    init_parser.add_argument("--species", default="human", help="Default species context.")
    init_parser.add_argument("--force", action="store_true", help="Overwrite an existing config.")
    init_parser.set_defaults(func=cmd_init)

    annotate_parser = subparsers.add_parser("annotate", help="Run CASSIA annotation from marker CSV data.")
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

    boost_parser = subparsers.add_parser("boost", help="Annotation boost helper commands.")
    boost_subparsers = boost_parser.add_subparsers(dest="boost_command", required=True)
    query_parser = boost_subparsers.add_parser(
        "query",
        help="Query marker statistics for genes, optionally within one cluster.",
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

    return parser


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
