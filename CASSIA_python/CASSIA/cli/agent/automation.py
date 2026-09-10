"""Automated clustering-and-annotation orchestration for ``cassia agent``.

Versioned on-disk checkpoints own scientific state and enforce mutation
budgets.  A coding agent may investigate marker evidence and propose bounded
merge/subcluster operations, but it can only change the object through audited
one-shot CLI transactions.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from ..backends import AgentCLIBackend, is_agent_backend
from . import direct
from .experiment import (
    read_audit_entries,
    render_automation_comparison,
    summarize_automation_audit,
    summarize_automation_run,
)
from .policy import AUTOMATION_STRATEGIES, strategy_policy
from .prompting import build_automation_prompt
from .state import ensure_workdir


AUTOMATION_SCHEMA_VERSION = "cassia.agent.automation.v1"


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _checkpoint(workdir: Path, path: Path) -> Dict[str, Any]:
    return direct.call(
        workdir, "checkpoint", path=str(path.resolve()), timeout=1800.0
    )["result"]


def _export_partition(workdir: Path, out_dir: Path, partition_id: str) -> Dict[str, Any]:
    return direct.call(
        workdir,
        "export",
        out_dir=str(out_dir.resolve()),
        partition_id=partition_id,
        timeout=600.0,
    )["result"]


def _build_manifest(args: Any, rds_path: Path, out_dir: Path) -> Dict[str, Any]:
    policy = strategy_policy(
        args.strategy,
        args.max_topology_edits,
        args.min_child_cells,
        args.max_children_per_split,
    )
    return {
        "schema_version": AUTOMATION_SCHEMA_VERSION,
        "created_at": _utc_now(),
        "updated_at": _utc_now(),
        "status": "running",
        "input": {
            "path": str(rds_path),
            "sha256": _sha256_file(rds_path),
            "size_bytes": rds_path.stat().st_size,
        },
        "out_dir": str(out_dir),
        "strategy": args.strategy,
        "policy": policy,
        "backend": args.backend,
        "model": args.model,
        "reasoning_effort": args.reasoning_effort,
        "parameters": {
            "species": args.species,
            "tissue": args.tissue,
            "additional_info": args.additional_info,
            "initial": args.initial,
            "n_var": args.n_var,
            "n_pcs": args.n_pcs,
            "resolution": args.resolution,
            "seed": args.seed,
            "algorithm": args.algorithm,
            "skip_scale": bool(args.skip_scale),
            "autozyme": bool(args.autozyme),
            "max_topology_edits": args.max_topology_edits,
            "min_child_cells": args.min_child_cells,
            "max_children_per_split": args.max_children_per_split,
            "max_agent_commands": args.max_agent_commands,
            "timeout": args.timeout,
        },
        "attempts": [],
        "artifacts": {},
    }


def _validate_run_args(args: Any) -> tuple[Path, Path]:
    rds_path = Path(args.rds).expanduser().resolve()
    if not rds_path.exists():
        raise FileNotFoundError(f"Seurat RDS not found: {rds_path}")
    if rds_path.suffix.lower() != ".rds":
        raise ValueError("cassia agent auto currently requires a Seurat .rds input")
    out_dir = Path(args.out).expanduser().resolve()
    if not is_agent_backend(args.backend):
        raise ValueError("cassia agent auto requires an agent CLI backend")
    if args.backend == "shell" and not args.command_template and not args.dry_run:
        raise ValueError("--command-template is required with --backend shell")
    if args.backend == "codex-cli":
        # The integrated controller is benchmarked with Astra/high, while
        # keeping both values explicitly overridable for reproduction.
        args.model = args.model or "gpt-6-astra"
        args.reasoning_effort = args.reasoning_effort or "high"
    strategy_policy(
        args.strategy,
        args.max_topology_edits,
        args.min_child_cells,
        args.max_children_per_split,
    )
    if args.max_agent_commands < 1:
        raise ValueError("max_agent_commands must be at least 1")
    if args.timeout < 1:
        raise ValueError("timeout must be at least 1 second")
    if args.initial not in {"existing", "compute"}:
        raise ValueError("initial must be 'existing' or 'compute'")
    return rds_path, out_dir


def _validate_resume_manifest(manifest: Mapping[str, Any], args: Any) -> None:
    """Refuse a resume that would silently change the scientific run."""
    saved_parameters = manifest.get("parameters") or {}
    current = {
        "strategy": args.strategy,
        "backend": args.backend,
        "model": args.model,
        "reasoning_effort": args.reasoning_effort,
        "species": args.species,
        "tissue": args.tissue,
        "additional_info": args.additional_info,
        "initial": args.initial,
        "n_var": args.n_var,
        "n_pcs": args.n_pcs,
        "resolution": args.resolution,
        "seed": args.seed,
        "algorithm": args.algorithm,
        "skip_scale": bool(args.skip_scale),
        "autozyme": bool(args.autozyme),
        "max_topology_edits": args.max_topology_edits,
        "min_child_cells": args.min_child_cells,
        "max_children_per_split": args.max_children_per_split,
        "max_agent_commands": args.max_agent_commands,
    }
    saved = {
        "strategy": manifest.get("strategy"),
        "backend": manifest.get("backend"),
        "model": manifest.get("model"),
        "reasoning_effort": manifest.get("reasoning_effort"),
        **{key: saved_parameters.get(key) for key in current if key in saved_parameters},
    }
    mismatches = [
        key for key, value in current.items()
        if key in saved and saved.get(key) != value
    ]
    if mismatches:
        joined = ", ".join(sorted(mismatches))
        raise ValueError(
            f"Resume configuration differs from the saved run ({joined}); "
            "use a new --out directory for a different experiment"
        )


def run_auto_agent(args: Any) -> int:
    """Run one bounded coding-agent controller over transactional CLI state."""
    rds_path, out_dir = _validate_run_args(args)
    manifest_path = out_dir / "automation_manifest.json"

    if manifest_path.exists() and not args.resume:
        raise FileExistsError(
            f"Automation run already exists at {out_dir}; use --resume or a new --out"
        )
    if args.resume and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        expected_hash = manifest.get("input", {}).get("sha256")
        actual_hash = _sha256_file(rds_path)
        if expected_hash and expected_hash != actual_hash:
            raise ValueError("Input RDS changed since the saved automation run")
        _validate_resume_manifest(manifest, args)
        manifest["updated_at"] = _utc_now()
        manifest["status"] = "running"
    else:
        manifest = _build_manifest(args, rds_path, out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    prompt_dir = out_dir / "prompts"
    raw_dir = out_dir / "raw"
    checkpoint_dir = out_dir / "checkpoints"
    provenance_dir = out_dir / "provenance"
    output_dir = out_dir / "outputs"
    for directory in (prompt_dir, raw_dir, checkpoint_dir, provenance_dir, output_dir):
        directory.mkdir(parents=True, exist_ok=True)

    workdir = ensure_workdir(out_dir)
    prompt = build_automation_prompt(
        workdir=workdir,
        strategy=args.strategy,
        species=args.species,
        tissue=args.tissue,
        additional_info=args.additional_info,
        max_topology_edits=args.max_topology_edits,
        max_agent_commands=args.max_agent_commands,
        min_child_cells=args.min_child_cells,
        max_children_per_split=args.max_children_per_split,
    )
    prompt_path = prompt_dir / "automation.md"
    prompt_path.write_text(prompt, encoding="utf-8")
    manifest["artifacts"]["prompt"] = str(prompt_path)

    if args.dry_run:
        manifest["status"] = "dry-run"
        manifest["updated_at"] = _utc_now()
        _atomic_write_json(manifest_path, manifest)
        print(f"Wrote {manifest_path}")
        print(f"Wrote {prompt_path}")
        return 0

    _atomic_write_json(manifest_path, manifest)
    baseline_path = checkpoint_dir / "baseline.rds"
    recovery_path = checkpoint_dir / "recovery.rds"
    attempt_started = _utc_now()
    try:
        restored_from: Optional[Path] = None
        if args.resume and direct.has_state(workdir):
            current = direct.state_status(workdir).get("current_state")
            restored_from = Path(str(current)) if current else None
        elif args.resume:
            for candidate in (recovery_path, baseline_path):
                if candidate.exists():
                    direct.restore_checkpoint(workdir, candidate, timeout=1800.0)
                    restored_from = candidate
                    break

        if restored_from is None:
            init_result = direct.call(
                workdir,
                "init",
                rds_path=str(rds_path),
                autozyme=bool(args.autozyme),
                timeout=1800.0,
            )["result"]
            if args.initial == "compute":
                direct.call(
                    workdir,
                    "preprocess",
                    n_var=args.n_var,
                    n_pcs=args.n_pcs,
                    resolution=args.resolution,
                    seed=args.seed,
                    algorithm=args.algorithm,
                    skip_scale=bool(args.skip_scale),
                    timeout=3600.0,
                )
            elif init_result.get("phase") != "ORIENT":
                raise ValueError(
                    "The input has no usable existing partition; rerun with --initial compute"
                )
            baseline = _checkpoint(workdir, baseline_path)
            baseline_export = _export_partition(
                workdir, provenance_dir, "baseline"
            )
            manifest["artifacts"]["baseline_checkpoint"] = baseline
            manifest["artifacts"]["baseline_partition"] = baseline_export

        policy = strategy_policy(
            args.strategy,
            args.max_topology_edits,
            args.min_child_cells,
            args.max_children_per_split,
        )
        direct.call(
            workdir,
            "policy",
            enabled=True,
            allowed_topology_actions=policy["allowed_topology_actions"],
            max_topology_edits=policy["max_topology_edits"],
            min_child_cells=policy["min_child_cells"],
            max_children_per_split=policy["max_children_per_split"],
            reset_counter=restored_from is None,
            lock_preprocess=True,
        )

        audit_path = workdir / "audit.jsonl"
        audit_start = len(read_audit_entries(audit_path))
        raw_path = raw_dir / "agent_final.txt"
        backend = AgentCLIBackend(
            args.backend,
            command_template=args.command_template,
            timeout_seconds=args.timeout,
            model=args.model,
            reasoning_effort=args.reasoning_effort,
            sandbox_mode="workspace-write" if args.backend == "codex-cli" else None,
        )
        response = backend.run(
            prompt=prompt,
            prompt_file=prompt_path,
            cwd=out_dir,
            context={
                "input": str(rds_path),
                "out": str(out_dir),
                "cluster": "all",
                "agent_output_file": str(raw_path),
            },
        )
        raw_path.write_text(response + ("\n" if response and not response.endswith("\n") else ""), encoding="utf-8")

        audit_summary = summarize_automation_audit(
            read_audit_entries(audit_path), start_index=audit_start
        )
        qa = direct.call(workdir, "qa")["result"]
        final_payload: Optional[Dict[str, Any]] = None
        if qa.get("pass"):
            final_payload = direct.call(
                workdir,
                "finalize",
                force=False,
                out_rds=str((output_dir / "annotated.rds").resolve()),
                out_tsv=str((output_dir / "annotation.tsv").resolve()),
                out_md=str((output_dir / "report.md").resolve()),
                timeout=1800.0,
            )["result"]

        recovery = _checkpoint(workdir, recovery_path)
        final_export = _export_partition(workdir, provenance_dir, "final")
        attempt = {
            "started_at": attempt_started,
            "completed_at": _utc_now(),
            "status": "completed" if qa.get("pass") else "needs_attention",
            "restored_from": str(restored_from) if restored_from else None,
            "audit": audit_summary,
            "qa": qa,
            "agent_metadata": backend.last_run_metadata,
            "raw_response": str(raw_path),
        }
        manifest["attempts"].append(attempt)
        manifest["artifacts"]["recovery_checkpoint"] = recovery
        manifest["artifacts"]["final_partition"] = final_export
        if final_payload is not None:
            manifest["artifacts"]["final_outputs"] = final_payload
        manifest["status"] = attempt["status"]
        manifest["updated_at"] = _utc_now()
        _atomic_write_json(manifest_path, manifest)

        print(f"Wrote {manifest_path}")
        print(f"Agent commands: {audit_summary['command_count']}")
        print(f"Topology edits: {audit_summary['topology_edit_count']}")
        print(f"QA: {'PASS' if qa.get('pass') else 'NEEDS ATTENTION'}")
        if final_payload:
            print(f"Wrote {final_payload.get('out_rds')}")
            print(f"Wrote {final_payload.get('out_tsv')}")
            print(f"Wrote {final_payload.get('out_md')}")
        return 0 if qa.get("pass") else 1
    except Exception as exc:
        if direct.has_state(workdir):
            try:
                manifest["artifacts"]["recovery_checkpoint"] = _checkpoint(
                    workdir, recovery_path
                )
            except Exception as checkpoint_exc:
                manifest["recovery_error"] = str(checkpoint_exc)
        manifest["attempts"].append(
            {
                "started_at": attempt_started,
                "completed_at": _utc_now(),
                "status": "failed",
                "error": str(exc),
            }
        )
        manifest["status"] = "failed"
        manifest["updated_at"] = _utc_now()
        _atomic_write_json(manifest_path, manifest)
        raise


__all__ = [
    "AUTOMATION_SCHEMA_VERSION",
    "AUTOMATION_STRATEGIES",
    "build_automation_prompt",
    "render_automation_comparison",
    "run_auto_agent",
    "summarize_automation_audit",
    "summarize_automation_run",
]
