"""CLI handlers and parser wiring for integrated automation experiments."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .policy import AUTOMATION_STRATEGIES


def _print_json(value: Any) -> None:
    print(json.dumps(value, indent=2, ensure_ascii=False, default=str))


def cmd_auto(args: argparse.Namespace) -> int:
    """Run the bounded clustering-and-annotation controller."""
    from .automation import run_auto_agent

    return run_auto_agent(args)


def cmd_compare(args: argparse.Namespace) -> int:
    """Compare saved automation runs using audited, deterministic metrics."""
    from .experiment import render_automation_comparison, summarize_automation_run

    run_dirs = [Path(value).expanduser().resolve() for value in args.runs]
    if args.json:
        _print_json([summarize_automation_run(path) for path in run_dirs])
        return 0
    report = render_automation_comparison(run_dirs)
    if args.out:
        output = Path(args.out).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report, encoding="utf-8")
        print(f"Wrote {output}")
    else:
        print(report, end="")
    return 0


def add_automation_subparsers(subparsers: argparse._SubParsersAction) -> None:
    """Register the ``agent auto`` and ``agent compare`` command surfaces."""
    auto = subparsers.add_parser(
        "auto",
        help="Run recoverable clustering + annotation with a coding-agent controller.",
        description=(
            "Run an end-to-end, audited clustering-and-annotation experiment. "
            "One-shot R transactions version the Seurat state, enforce topology "
            "budgets, and write checkpoints plus exact cell memberships.\n\n"
            "Examples:\n"
            "  cassia agent auto obj.rds --out runs/fixed --strategy fixed\n"
            "  cassia agent auto obj.rds --out runs/adaptive --strategy adaptive \\\n"
            "    --model gpt-6-astra --reasoning-effort high --max-topology-edits 4\n"
            "  cassia agent auto obj.rds --out runs/preview --dry-run"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    auto.add_argument("rds", help="Input Seurat .rds.")
    auto.add_argument("--out", required=True, help="New or resumable run directory.")
    auto.add_argument(
        "--strategy",
        choices=AUTOMATION_STRATEGIES,
        default="conservative",
        help="Topology strategy (default conservative).",
    )
    auto.add_argument(
        "--backend",
        default="codex-cli",
        choices=["codex-cli", "claude-cli", "cursor-agent", "opencode", "shell"],
        help="Coding-agent backend (default codex-cli).",
    )
    auto.add_argument(
        "--model",
        help="Agent model. Codex defaults to gpt-6-astra for this command.",
    )
    auto.add_argument(
        "--reasoning-effort",
        choices=["minimal", "low", "medium", "high", "xhigh", "max"],
        help="Supported agent reasoning effort. Codex defaults to high.",
    )
    auto.add_argument(
        "--command-template",
        help="Shell command template; required only with --backend shell.",
    )
    auto.add_argument(
        "--initial",
        choices=["existing", "compute"],
        default="existing",
        help="Use existing identities or compute a Seurat partition (default existing).",
    )
    auto.add_argument("--species", default="human", help="Species context (default human).")
    auto.add_argument("--tissue", default="tissue blind", help="Tissue context (default tissue blind).")
    auto.add_argument("--additional-info", help="Optional dataset context for the agent.")
    auto.add_argument("--n-var", type=int, default=2000, help="Variable features for --initial compute.")
    auto.add_argument("--n-pcs", type=int, default=30, help="Principal components for --initial compute.")
    auto.add_argument("--resolution", type=float, default=0.5, help="Initial Seurat resolution.")
    auto.add_argument("--seed", type=int, default=17, help="Initial clustering seed.")
    auto.add_argument(
        "--algorithm",
        type=int,
        choices=[1, 2, 3, 4],
        default=1,
        help="Seurat FindClusters algorithm code.",
    )
    auto.add_argument("--skip-scale", action="store_true", help="Skip ScaleData when computing the initial partition.")
    auto.add_argument("--autozyme", action="store_true", help="Activate autozyme patches in R transactions.")
    auto.add_argument(
        "--max-topology-edits",
        type=int,
        default=4,
        help="Maximum merge/subcluster operations; fixed forces 0, conservative caps at 2.",
    )
    auto.add_argument(
        "--max-agent-commands",
        type=int,
        default=80,
        help="Controller command budget stated in the agent contract (default 80).",
    )
    auto.add_argument(
        "--min-child-cells",
        type=int,
        default=5,
        help="Reject a split if any child has fewer cells (default 5).",
    )
    auto.add_argument(
        "--max-children-per-split",
        type=int,
        default=4,
        help="Reject overly fragmented splits; conservative caps this at 3.",
    )
    auto.add_argument("--timeout", type=int, default=3600, help="Agent timeout in seconds.")
    auto.add_argument("--resume", action="store_true", help="Restore recovery/baseline state and retry the same run.")
    auto.add_argument("--dry-run", action="store_true", help="Write the manifest and prompt without starting R or an agent.")
    auto.set_defaults(func=cmd_auto)

    compare = subparsers.add_parser(
        "compare",
        help="Compare completed automation runs from their audited artifacts.",
        description=(
            "Summarize saved runs without calling an LLM or rerunning Seurat. "
            "Metrics include committed topology edits, cluster fragmentation, "
            "QA, and labeled-cell coverage."
        ),
    )
    compare.add_argument("runs", nargs="+", help="Automation run directories.")
    compare.add_argument("--out", help="Optional Markdown output path.")
    compare.add_argument("--json", action="store_true", help="Print structured metrics as JSON.")
    compare.set_defaults(func=cmd_compare)


__all__ = ["add_automation_subparsers", "cmd_auto", "cmd_compare"]
