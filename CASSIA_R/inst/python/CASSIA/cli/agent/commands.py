"""Subcommand handlers for `cassia agent ...`."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from . import client, direct, launcher
from .automation_cli import add_automation_subparsers
from .state import ensure_workdir, read_session, workdir_for


# ── pretty printers ─────────────────────────────────────────────────────

def _print_json(obj: Any) -> None:
    print(json.dumps(obj, indent=2, ensure_ascii=False, default=str))


def _print_clusters_table(clusters: List[Dict[str, Any]]) -> None:
    if not clusters:
        print("  (no clusters)")
        return
    print(f"  {'uuid':<10} {'seurat':<8} {'cells':>10}  alive  origin")
    for c in clusters:
        print(
            f"  {c.get('uuid','?'):<10} "
            f"{str(c.get('seurat_id') or '-'):<8} "
            f"{c.get('n_cells', 0):>10,}  "
            f"{('yes' if c.get('alive') else 'no'):<5}  "
            f"{c.get('origin','?')}"
        )


def _print_markers(payload: Dict[str, Any]) -> None:
    scope = payload.get("scope")
    markers = payload.get("markers") or []
    if scope == "all":
        if not markers:
            print("  (no markers)")
            return
        # Group by cluster
        by_cluster: Dict[str, List[Dict[str, Any]]] = {}
        for row in markers:
            key = f"{row.get('cluster_uuid','?')} (id {row.get('cluster','?')})"
            by_cluster.setdefault(key, []).append(row)
        for key, rows in by_cluster.items():
            gene_list = "  ".join(r["gene"] for r in rows)
            print(f"  {key}: {gene_list}")
        return
    # single-cluster table
    ident1 = payload.get("ident_1_uuid") or payload.get("ident_1")
    ident2 = payload.get("ident_2", "rest")
    print(f"  cluster={ident1}  vs={ident2}  n={payload.get('n')}")
    print(f"  {'gene':<14} {'avg_log2FC':>10}  {'pct.1':>6}  {'pct.2':>6}  p_val_adj")
    for r in markers:
        print(
            f"  {r.get('gene',''):<14} "
            f"{float(r.get('avg_log2FC',0)):>10.3f}  "
            f"{float(r.get('pct.1',0)):>6.3f}  "
            f"{float(r.get('pct.2',0)):>6.3f}  "
            f"{float(r.get('p_val_adj',0)):.3e}"
        )


# ── handlers ────────────────────────────────────────────────────────────

def _live_daemon(workdir: Path) -> bool:
    session = read_session(workdir)
    return bool(session and launcher._is_daemon_alive(session))


def _call(workdir: Path, op: str, *, timeout: float = 600.0, **kwargs: Any) -> Dict[str, Any]:
    """Use an explicitly started daemon, otherwise one-shot transactions."""
    if _live_daemon(workdir):
        return client.call(workdir, op, timeout=timeout, **kwargs)
    return direct.call(workdir, op, timeout=timeout, **kwargs)

def cmd_init(args: argparse.Namespace) -> int:
    rds_path = Path(args.rds).resolve()
    if not rds_path.exists():
        raise SystemExit(f"rds not found: {rds_path}")
    wd = ensure_workdir(args.workdir or rds_path.parent)
    print(f"[cassia agent] workdir: {wd}")
    if args.daemon:
        print(f"[cassia agent] spawning optional daemon accelerator...")
        session = launcher.spawn(wd)
        print(f"  daemon pid={session['pid']} port={session['port']}")
    else:
        print("[cassia agent] execution: one-shot transactions (no daemon)")
    print(f"[cassia agent] loading {rds_path}...")
    resp = _call(wd, "init", rds_path=str(rds_path),
                 autozyme=bool(args.autozyme), timeout=1800.0)
    payload = resp["result"]
    if args.json:
        _print_json(payload)
    else:
        print(f"  loaded: {payload.get('n_cells'):,} cells × "
              f"{payload.get('n_features'):,} features, "
              f"{payload.get('n_clusters_alive')} alive cluster(s), "
              f"phase={payload.get('phase')}")
        if payload.get("phase") == "INIT":
            print("  [warn] no clusters found — run `cassia agent preprocess` "
                  "(not yet in iteration 1) or load a preprocessed .rds")
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    resp = _call(wd, "status")
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    print(f"phase:           {payload.get('phase')}")
    print(f"step:            {payload.get('step')}")
    print(f"rds:             {payload.get('rds_path')}")
    print(f"cells:           {payload.get('n_cells'):,}")
    print(f"clusters alive:  {payload.get('n_clusters_alive')}")
    print(f"seurat memory:   {payload.get('seurat_mb')} MB")
    print(f"autozyme active: {payload.get('autozyme_active')}")
    return 0


def cmd_checkpoint(args: argparse.Namespace) -> int:
    """Export the complete persisted Seurat/annotation state."""
    wd = workdir_for(args.workdir or Path.cwd())
    requested = str(Path(args.out).resolve()) if args.out else None
    resp = _call(wd, "checkpoint", path=requested, timeout=1800.0)
    payload = resp["result"]
    if args.json:
        _print_json(payload)
    else:
        print(f"checkpoint: {payload.get('path')}")
        print(
            f"  step={payload.get('step')} phase={payload.get('phase')} "
            f"clusters={payload.get('n_clusters_alive')} md5={payload.get('md5')}"
        )
    return 0


def cmd_restore(args: argparse.Namespace) -> int:
    """Restore a complete checkpoint as a new persisted state version."""
    wd = ensure_workdir(args.workdir or Path.cwd())
    if _live_daemon(wd):
        resp = client.call(
            wd,
            "restore",
            path=str(Path(args.checkpoint).resolve()),
            timeout=1800.0,
        )
    else:
        resp = direct.restore_checkpoint(
            wd, Path(args.checkpoint).resolve(), timeout=1800.0
        )
    payload = resp["result"]
    if args.json:
        _print_json(payload)
    else:
        print(f"restored: {payload.get('path')}")
        print(
            f"  step={payload.get('step')} phase={payload.get('phase')} "
            f"clusters={payload.get('n_clusters_alive')}"
        )
    return 0


def cmd_export(args: argparse.Namespace) -> int:
    """Export versioned per-cell memberships and the alive cluster registry."""
    wd = workdir_for(args.workdir or Path.cwd())
    out_dir = Path(args.out).resolve() if args.out else wd / "artifacts"
    resp = _call(
        wd,
        "export",
        out_dir=str(out_dir),
        partition_id=args.partition_id,
        timeout=600.0,
    )
    payload = resp["result"]
    if args.json:
        _print_json(payload)
    else:
        print(f"partition: {payload.get('partition_id')}")
        print(f"  wrote: {payload.get('memberships_path')}")
        print(f"  wrote: {payload.get('clusters_path')}")
        print(
            f"  cells={payload.get('n_cells')} clusters={payload.get('n_clusters')} "
            f"md5={payload.get('memberships_md5')}"
        )
    return 0


def cmd_clusters(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    resp = _call(wd, "clusters", alive_only=not args.all)
    clusters = resp["result"]
    if args.json:
        _print_json(clusters)
        return 0
    print(f"clusters ({len(clusters)}):")
    _print_clusters_table(clusters)
    return 0


def cmd_markers(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    resp = _call(
        wd, "markers",
        cluster=args.cluster,
        vs=args.vs,
        top_n=args.top,
        only_pos=args.only_pos,
    )
    elapsed = resp.get("elapsed_s", 0.0)
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    print(f"(elapsed {elapsed:.2f}s)")
    _print_markers(payload)
    return 0


def cmd_preprocess(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    print(f"[cassia agent] preprocessing (Normalize→VarFeat→Scale→PCA→Neighbors→Clusters)...")
    resp = _call(
        wd, "preprocess",
        n_var=args.n_var,
        n_pcs=args.n_pcs,
        resolution=args.resolution,
        seed=args.seed,
        algorithm=args.algorithm,
        skip_scale=args.skip_scale,
        timeout=1800.0,
    )
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    print(f"  done in {resp.get('elapsed_s', 0):.1f}s")
    print(f"  resolution={payload.get('resolution')} → {payload.get('n_clusters')} clusters")
    print(f"  phase={payload.get('phase')}")
    return 0


def cmd_merge(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    resp = _call(wd, "merge", uuids=args.uuids, reason=args.reason)
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    print(f"  [merged] {len(payload.get('merged_uuids', []))} clusters → "
          f"new uuid {payload.get('new_uuid')} (n={payload.get('n_cells'):,})")
    print(f"  → {payload.get('note')}")
    return 0


def cmd_subcluster(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    resolutions = args.resolution_ladder if args.auto_resolution else [args.resolution]
    attempts: List[Dict[str, Any]] = []
    resp: Optional[Dict[str, Any]] = None
    last_error: Optional[Exception] = None
    for resolution in resolutions:
        try:
            resp = _call(
                wd, "subcluster",
                uuid=args.uuid,
                resolution=resolution,
                graph=args.graph,
                algorithm=args.algorithm,
                reason=args.reason,
                timeout=600.0,
            )
            attempts.append({"resolution": resolution, "accepted": True})
            break
        except client.DaemonError as exc:
            attempts.append(
                {"resolution": resolution, "accepted": False, "reason": str(exc)}
            )
            last_error = exc
            if not args.auto_resolution or "transaction rolled back" not in str(exc):
                raise
    if resp is None:
        tried = ", ".join(str(item["resolution"]) for item in attempts)
        raise client.DaemonError(
            f"auto-resolution exhausted without an acceptable split ({tried}): {last_error}"
        )
    payload = resp["result"]
    payload["resolution_attempts"] = attempts
    if args.json:
        _print_json(payload)
        return 0
    children = payload.get("children") or []
    print(f"  [subclustered] {args.uuid} → {len(children)} children "
          f"(res={payload.get('resolution')})")
    if args.auto_resolution:
        rejected = [item for item in attempts if not item["accepted"]]
        print(
            f"  auto-resolution: {len(rejected)} rejected, "
            f"accepted {payload.get('resolution')}"
        )
    for ch in children:
        print(f"    {ch.get('uuid'):<10} id={ch.get('seurat_id'):<10} "
              f"n={ch.get('n_cells'):,}")
    print(f"  → {payload.get('note')}")
    return 0


def cmd_gene(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    # Build a single comma-separated string from positional genes; the R worker
    # also accepts comma/whitespace-separated tokens within any entry.
    raw_genes = " ".join(args.genes)
    resp = _call(wd, "gene", genes=raw_genes)
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    if payload.get("cache_warning"):
        print(f"  [warn] {payload['cache_warning']}")
    cached = payload.get("n_clusters_with_cache", 0)
    results = payload.get("results") or []
    if not results:
        print("  (no genes parsed)")
        return 0
    for i, r in enumerate(results):
        gene = r.get("gene", "?")
        hits = r.get("hits") or []
        if not hits:
            print(f"  {gene}: no hits in {cached} cached cluster(s)")
            if i < len(results) - 1:
                print()
            continue
        print(f"  {gene}: enriched in {len(hits)}/{cached} cached cluster(s)")
        print(f"  {'uuid':<10} {'id':<5} {'cells':>8}  {'rank':>4}  {'log2FC':>7}  {'pct.1':>6}  {'pct.2':>6}  label")
        for h in hits:
            lbl = h.get("cluster_label") or "-"
            if lbl in (None, "NA"):
                lbl = "-"
            print(
                f"  {h.get('cluster_uuid','?'):<10} "
                f"{str(h.get('seurat_id') or '-'):<5} "
                f"{h.get('n_cells', 0):>8,}  "
                f"{h.get('rank', 0):>4}  "
                f"{float(h.get('avg_log2FC', 0) or 0):>7.3f}  "
                f"{float(h.get('pct.1', 0) or 0):>6.3f}  "
                f"{float(h.get('pct.2', 0) or 0):>6.3f}  "
                f"{lbl}"
            )
        if i < len(results) - 1:
            print()
    return 0


def cmd_prompt(args: argparse.Namespace) -> int:
    from pathlib import Path as _P
    here = _P(__file__).resolve().parent
    prompt_path = here / "system_prompt.md"
    if not prompt_path.exists():
        raise SystemExit(f"system prompt template missing: {prompt_path}")
    sys.stdout.write(prompt_path.read_text(encoding="utf-8"))
    return 0


def cmd_label(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    markers = _parse_marker_list(args.markers)
    resp = _call(
        wd, "label",
        uuid=args.uuid,
        name=args.name,
        confidence=args.confidence,
        markers=markers,
        reason=args.reason,
        skip=args.skip,
    )
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    if args.skip:
        print(f"  [skipped] {args.uuid}: {args.reason}")
    else:
        print(f"  [labeled] {payload.get('uuid')}: \"{payload.get('name')}\" "
              f"(conf={payload.get('confidence')}, {payload.get('n_markers')} markers)")
    return 0


def cmd_unlabel(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    resp = _call(wd, "unlabel", uuid=args.uuid)
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    if payload.get("removed"):
        prev = payload.get("previous") or {}
        print(f"  [unlabeled] {args.uuid} (was: {prev.get('name')})")
    else:
        print(f"  [no-op] {args.uuid}: {payload.get('msg')}")
    return 0


def cmd_qa(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    resp = _call(wd, "qa")
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0 if payload.get("pass") else 1
    n_alive = payload.get("n_alive", 0)
    n_lab = payload.get("n_labeled", 0)
    n_skip = payload.get("n_skipped", 0)
    n_err = payload.get("n_errors", 0)
    n_warn = payload.get("n_warnings", 0)
    print(f"qa: {'PASS' if payload.get('pass') else 'FAIL'}")
    print(f"  alive={n_alive}  labeled={n_lab}  skipped={n_skip}  unresolved={n_alive - n_lab - n_skip}")
    print(f"  errors={n_err}  warnings={n_warn}")
    for f in payload.get("findings") or []:
        print(f"  [{f.get('severity','?').upper():<5}] {f.get('code','?'):<24} {f.get('msg','')}")
    return 0 if payload.get("pass") else 1


def cmd_finalize(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    kw = {"force": args.force}
    if args.out:
        out = Path(args.out).resolve()
        kw["out_rds"] = str(out)
    resp = _call(wd, "finalize", **kw, timeout=600.0)
    payload = resp["result"]
    if args.json:
        _print_json(payload)
        return 0
    qa = payload.get("qa") or {}
    print(f"qa: {'PASS' if qa.get('pass') else 'FORCED'} "
          f"(errors={qa.get('n_errors', 0)}, warnings={qa.get('n_warnings', 0)})")
    print(f"  wrote: {payload.get('out_rds')}")
    print(f"  wrote: {payload.get('out_tsv')}")
    print(f"  wrote: {payload.get('out_md')}")
    print(f"  labeled={payload.get('n_labeled')}  skipped={payload.get('n_skipped')}")
    return 0


def _parse_marker_list(s: Optional[str]) -> List[str]:
    if not s:
        return []
    parts: List[str] = []
    for piece in s.replace(";", ",").split(","):
        piece = piece.strip()
        if piece:
            parts.append(piece)
    return parts


def _parse_resolution_ladder(value: str) -> List[float]:
    try:
        resolutions = [float(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("resolution ladder must contain numbers") from exc
    if not resolutions or any(item <= 0 for item in resolutions):
        raise argparse.ArgumentTypeError("resolution ladder values must be positive")
    if resolutions != sorted(set(resolutions)):
        raise argparse.ArgumentTypeError("resolution ladder must be unique and increasing")
    return resolutions


def cmd_history(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    audit_path = wd / "audit.jsonl"
    if not audit_path.exists():
        raise SystemExit(f"no audit log at {audit_path}")
    entries: List[Dict[str, Any]] = []
    with audit_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    # Filter: only op events by default (skip daemon_start / daemon_stop).
    if not args.all:
        entries = [e for e in entries if e.get("event") == "op"]
    if args.last:
        entries = entries[-args.last:]
    if args.json:
        _print_json(entries)
        return 0
    for e in entries:
        ev = e.get("event")
        if ev != "op":
            print(f"  [{e.get('ts','?'):<24}] {ev}: {json.dumps({k:v for k,v in e.items() if k not in ('event','ts')}, default=str)}")
            continue
        ts = e.get("ts", "?")
        op = e.get("op", "?")
        ok = "OK " if e.get("ok") else "ERR"
        elapsed = e.get("elapsed_s", 0.0)
        op_args = e.get("args") or {}
        summary = e.get("summary") or {}
        arg_str = _fmt_args(op, op_args)
        sum_str = _fmt_summary(op, summary, e.get("error"))
        print(f"  [{ts:<24}] {ok} {elapsed:>6.2f}s  {op:<10}  {arg_str}")
        if sum_str:
            print(f"      → {sum_str}")
    return 0


def _fmt_args(op: str, args: Dict[str, Any]) -> str:
    if not args:
        return ""
    if op == "markers":
        parts = []
        if args.get("cluster"):
            parts.append(f"cluster={args['cluster']}")
        if args.get("vs"):
            parts.append(f"vs={args['vs']}")
        if args.get("top_n"):
            parts.append(f"top={args['top_n']}")
        if args.get("only_pos"):
            parts.append("only_pos")
        return " ".join(parts)
    if op == "init":
        return f"rds={args.get('rds_path','?')}"
    if op == "clusters":
        return "alive_only" if args.get("alive_only") else "all"
    return json.dumps(args, default=str)


def _fmt_summary(op: str, summary: Dict[str, Any], error: Optional[str]) -> str:
    if error:
        return f"error: {error}"
    if not summary:
        return ""
    if op == "markers":
        scope = summary.get("scope")
        if scope == "all":
            n = summary.get("n_rows", 0)
            return f"all clusters, {n} rows"
        n = summary.get("n_rows", 0)
        top = summary.get("top_genes") or []
        top_str = ", ".join(top[:5])
        ident1 = summary.get("ident_1_uuid") or summary.get("ident_1")
        vs = summary.get("ident_2", "rest")
        return f"{ident1} vs {vs}: {n} rows; top: {top_str}"
    if op == "init":
        return (f"{summary.get('n_cells'):,} cells, "
                f"{summary.get('n_clusters_alive')} clusters, "
                f"phase={summary.get('phase')}")
    if op == "status":
        return (f"phase={summary.get('phase')}, "
                f"step={summary.get('step')}, "
                f"clusters_alive={summary.get('n_clusters_alive')}")
    if op == "clusters":
        return f"{summary.get('n_returned')} clusters"
    return json.dumps(summary, default=str)


def cmd_close(args: argparse.Namespace) -> int:
    wd = workdir_for(args.workdir or Path.cwd())
    session = read_session(wd)
    if session is None:
        if direct.has_state(wd):
            print(f"no daemon is running; persisted state retained at {wd}")
        else:
            print(f"no active state at {wd}")
        return 0
    ok = launcher.shutdown(wd)
    print(f"daemon shut down: {ok}")
    return 0


# ── argparse wiring ─────────────────────────────────────────────────────

def add_agent_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    agent = subparsers.add_parser(
        "agent",
        help="Transactional agent CLI for iterative clustering + annotation.",
        description=(
            "Drive versioned one-shot R transactions for iterative clustering + "
            "marker investigation. State persists in <workdir>/.cassia/.\n\n"
            "Typical workflow:\n"
            "  cassia agent init obj.rds            # load versioned state\n"
            "  cassia agent markers all --top 15    # populate marker cache\n"
            "  cassia agent gene CD3D CD8A NKG7     # batch reverse-lookup\n"
            "  cassia agent label <uuid> \"<name>\" --markers M1,M2,M3 --confidence H --reason \"...\"\n"
            "  cassia agent qa                      # rule-checks\n"
            "  cassia agent finalize                # write annotated.rds + report\n"
            "  cassia agent close                   # stop optional accelerator\n\n"
            "Read `cassia agent prompt` for the full system prompt the agent should follow."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = agent.add_subparsers(dest="agent_command", required=True)
    add_automation_subparsers(sub)

    p_init = sub.add_parser(
        "init",
        help="Load a Seurat .rds into versioned state.",
        description=(
            "Load a Seurat object into transactional on-disk state. No background "
            "process is started unless --daemon is requested.\n\n"
            "Examples:\n"
            "  cassia agent init obj.rds\n"
            "  cassia agent init obj.rds --workdir runs/my_session\n"
            "  cassia agent init obj.rds --autozyme\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_init.add_argument("rds", help="Path to Seurat .rds.")
    p_init.add_argument("--workdir", help="Workdir for .cassia/ state. Defaults to the .rds parent.")
    p_init.add_argument("--autozyme", action="store_true", help="Activate autozyme patches.")
    p_init.add_argument(
        "--daemon",
        action="store_true",
        help="Keep R resident as an optional accelerator; default is one-shot transactions.",
    )
    p_init.add_argument("--json", action="store_true", help="Print init result as JSON.")
    p_init.set_defaults(func=cmd_init)

    p_status = sub.add_parser("status", help="Report current persisted state.")
    p_status.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_status.add_argument("--json", action="store_true", help="Print status as JSON.")
    p_status.set_defaults(func=cmd_status)

    p_checkpoint = sub.add_parser(
        "checkpoint",
        help="Persist the complete Seurat, cluster, marker-cache, and label state.",
    )
    p_checkpoint.add_argument("--out", help="Checkpoint .rds path. Defaults to .cassia/checkpoints/step_N.rds.")
    p_checkpoint.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_checkpoint.add_argument("--json", action="store_true", help="Print checkpoint metadata as JSON.")
    p_checkpoint.set_defaults(func=cmd_checkpoint)

    p_restore = sub.add_parser(
        "restore",
        help="Restore a complete checkpoint as a new state version.",
    )
    p_restore.add_argument("checkpoint", help="Checkpoint .rds path.")
    p_restore.add_argument("--workdir", help="Workdir for .cassia/ state. Defaults to the current directory.")
    p_restore.add_argument("--json", action="store_true", help="Print restore metadata as JSON.")
    p_restore.set_defaults(func=cmd_restore)

    p_export = sub.add_parser(
        "export",
        help="Export versioned cell memberships and the alive cluster registry.",
    )
    p_export.add_argument("--out", help="Artifact directory. Defaults to .cassia/artifacts.")
    p_export.add_argument("--partition-id", help="Stable partition identifier. Defaults to the current step.")
    p_export.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_export.add_argument("--json", action="store_true", help="Print export metadata as JSON.")
    p_export.set_defaults(func=cmd_export)

    p_cl = sub.add_parser("clusters", help="List clusters and UUIDs.")
    p_cl.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_cl.add_argument("--all", action="store_true", help="Include dead clusters.")
    p_cl.add_argument("--json", action="store_true", help="Print clusters as JSON.")
    p_cl.set_defaults(func=cmd_clusters)

    p_mk = sub.add_parser(
        "markers",
        help="Query markers for one cluster or all clusters.",
        description=(
            "Find differentially-expressed markers. Caches top-50 per cluster\n"
            "so subsequent `gene` and `qa` calls don't re-run FindMarkers.\n\n"
            "Examples:\n"
            "  cassia agent markers all --top 15            # orient phase\n"
            "  cassia agent markers <uuid> --top 20         # one cluster vs rest\n"
            "  cassia agent markers <uuid_A> --vs <uuid_B>  # pairwise comparison\n"
            "  cassia agent markers <uuid> --only-pos       # positive markers only\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_mk.add_argument("cluster", help="Cluster UUID, raw seurat id, or 'all'.")
    p_mk.add_argument("--vs", help="Compare against cluster UUID(s) (comma-separated). Defaults to 'rest'.")
    p_mk.add_argument("--top", type=int, default=20, help="Top N markers (per cluster if 'all').")
    p_mk.add_argument("--only-pos", action="store_true", help="Only positive markers.")
    p_mk.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_mk.add_argument("--json", action="store_true", help="Print result as JSON.")
    p_mk.set_defaults(func=cmd_markers)

    p_pre = sub.add_parser(
        "preprocess",
        help="Run NormalizeData → FindVariableFeatures → ScaleData → PCA → Neighbors → Clusters.",
        description=(
            "Run a default Seurat preprocessing pipeline on the loaded object.\n"
            "Use this when the .rds is raw (no PCA / no clusters). Refuses if\n"
            "any cluster already has a label (those would be invalidated).\n\n"
            "Examples:\n"
            "  cassia agent preprocess\n"
            "  cassia agent preprocess --resolution 0.8 --n-pcs 20\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_pre.add_argument("--n-var", type=int, default=2000, help="N variable features (default 2000).")
    p_pre.add_argument("--n-pcs", type=int, default=30, help="N PCs (default 30).")
    p_pre.add_argument("--resolution", type=float, default=0.5, help="Cluster resolution (default 0.5).")
    p_pre.add_argument("--seed", type=int, default=17, help="Random seed passed to Seurat clustering (default 17).")
    p_pre.add_argument(
        "--algorithm",
        type=int,
        choices=[1, 2, 3, 4],
        default=1,
        help="Seurat FindClusters algorithm code (default 1).",
    )
    p_pre.add_argument("--skip-scale", action="store_true", help="Skip ScaleData (faster, lower-quality PCA).")
    p_pre.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_pre.add_argument("--json", action="store_true", help="Print result as JSON.")
    p_pre.set_defaults(func=cmd_preprocess)

    p_merge = sub.add_parser(
        "merge",
        help="Merge two or more clusters into one new cluster.",
        description=(
            "Combine cells from multiple clusters into a single new cluster.\n"
            "The merged clusters become dead; a new UUID is born with all combined\n"
            "cells. Refuses if any merged cluster has a label (unlabel first).\n\n"
            "Examples:\n"
            "  cassia agent merge <uuid_A> <uuid_B> --reason \"near-duplicates, only cycling differs\"\n"
            "  cassia agent merge <uuid_A> <uuid_B> <uuid_C> --reason \"three subsets of same lineage\"\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_merge.add_argument("uuids", nargs="+", help="≥ 2 cluster UUIDs to merge.")
    p_merge.add_argument("--reason", required=True, help="Why merge (recorded in audit).")
    p_merge.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_merge.add_argument("--json", action="store_true", help="Print result as JSON.")
    p_merge.set_defaults(func=cmd_merge)

    p_sub = sub.add_parser(
        "subcluster",
        help="Split a cluster into sub-clusters using FindSubCluster.",
        description=(
            "Re-run clustering on cells of a single cluster only, then split it\n"
            "into children. Parent dies, children are born with UUIDs and\n"
            "seurat_ids of the form '<parent>_0', '<parent>_1', etc.\n\n"
            "Examples:\n"
            "  cassia agent subcluster <uuid> --resolution 0.3 --reason \"T+NK mix in top markers\"\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_sub.add_argument("uuid", help="Cluster UUID to subcluster.")
    p_sub.add_argument("--resolution", type=float, default=0.3,
                       help="Subcluster resolution (default 0.3, conservative).")
    p_sub.add_argument(
        "--auto-resolution",
        action="store_true",
        help="Try an increasing resolution ladder and commit the first policy-valid split.",
    )
    p_sub.add_argument(
        "--resolution-ladder",
        type=_parse_resolution_ladder,
        default=[0.3, 0.6, 0.8, 1.0],
        metavar="R1,R2,...",
        help="Ladder used by --auto-resolution (default 0.3,0.6,0.8,1.0).",
    )
    p_sub.add_argument("--graph", help="Seurat graph name. Defaults to the first *_snn graph.")
    p_sub.add_argument(
        "--algorithm",
        type=int,
        choices=[1, 2, 3, 4],
        default=1,
        help="Seurat FindSubCluster algorithm code (default 1).",
    )
    p_sub.add_argument("--reason", required=True, help="Why subcluster (recorded in audit).")
    p_sub.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_sub.add_argument("--json", action="store_true", help="Print result as JSON.")
    p_sub.set_defaults(func=cmd_subcluster)

    p_gene = sub.add_parser(
        "gene",
        help="Reverse-query: where are these genes enriched across clusters?",
        description=(
            "Look up one or many genes in the per-cluster cached FindMarkers tables.\n"
            "Pass multiple genes as separate args, comma-separated, or whitespace-separated:\n"
            "  cassia agent gene CD3D CD8A NKG7\n"
            "  cassia agent gene CD3D,CD8A,NKG7\n"
            "  cassia agent gene \"CD3D CD8A NKG7\"\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_gene.add_argument("genes", nargs="+", help="Gene symbol(s); comma- or whitespace-separated tokens are split.")
    p_gene.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_gene.add_argument("--json", action="store_true", help="Print results as JSON.")
    p_gene.set_defaults(func=cmd_gene)

    p_prompt = sub.add_parser("prompt", help="Print the recommended agent system prompt to stdout.")
    p_prompt.set_defaults(func=cmd_prompt)

    p_label = sub.add_parser(
        "label",
        help="Assign a cell-type label to a cluster (or --skip with reason).",
        description=(
            "Assign a label with cited evidence markers, or skip with reason.\n"
            "Confidence levels:  H ≥ 3 markers, M ≥ 2, L ≥ 1.\n"
            "QA later checks every cited marker against the cluster's top-50 cache.\n\n"
            "Examples:\n"
            "  cassia agent label <uuid> \"CD8 effector T cell\" \\\n"
            "      --markers GZMH,NKG7,GNLY,CD8A --confidence H \\\n"
            "      --reason \"GZMH+ NKG7+ GNLY+ CD8A+ canonical CD8 Teff\"\n"
            "  cassia agent label <uuid> --skip --reason \"epithelial contamination (KRT5+)\"\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_label.add_argument("uuid", help="Cluster UUID.")
    p_label.add_argument("name", nargs="?", help="Cell-type label (omit when using --skip).")
    p_label.add_argument("--markers", help="Comma-separated evidence markers. H needs ≥3, M ≥2, L ≥1.")
    p_label.add_argument("--confidence", default="M", choices=["H", "M", "L", "h", "m", "l"],
                         help="Label confidence (default M).")
    p_label.add_argument("--reason", default="", help="Free-text rationale. Required for L or --skip.")
    p_label.add_argument("--skip", action="store_true",
                         help="Mark this cluster as deliberately unresolved instead of labeling.")
    p_label.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_label.add_argument("--json", action="store_true", help="Print result as JSON.")
    p_label.set_defaults(func=cmd_label)

    p_unlabel = sub.add_parser("unlabel", help="Remove a previously assigned label.")
    p_unlabel.add_argument("uuid", help="Cluster UUID.")
    p_unlabel.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_unlabel.add_argument("--json", action="store_true", help="Print result as JSON.")
    p_unlabel.set_defaults(func=cmd_unlabel)

    p_qa = sub.add_parser("qa", help="Run rule-based QA checks on current labels.")
    p_qa.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_qa.add_argument("--json", action="store_true", help="Print findings as JSON.")
    p_qa.set_defaults(func=cmd_qa)

    p_final = sub.add_parser(
        "finalize",
        help="Write annotated .rds + report. Refuses unless qa passes.",
        description=(
            "Run qa; if it passes, write outputs/annotated.rds (with cassia_label\n"
            "+ cassia_cluster_uuid metadata columns), annotation.tsv, and report.md.\n\n"
            "Examples:\n"
            "  cassia agent finalize\n"
            "  cassia agent finalize --out outputs/pbmc_v1.rds\n"
            "  cassia agent finalize --force            # bypass qa errors (loudly!)\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_final.add_argument("--out", help="Output annotated .rds path (default outputs/annotated.rds).")
    p_final.add_argument("--force", action="store_true", help="Skip qa pass requirement.")
    p_final.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_final.add_argument("--json", action="store_true", help="Print result as JSON.")
    p_final.set_defaults(func=cmd_finalize)

    p_hist = sub.add_parser("history", help="Show the persisted audit log of every operation.")
    p_hist.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_hist.add_argument("--last", type=int, help="Show only the last N op entries.")
    p_hist.add_argument("--all", action="store_true", help="Include non-operation events too.")
    p_hist.add_argument("--json", action="store_true", help="Print full entries as JSON.")
    p_hist.set_defaults(func=cmd_history)

    p_close = sub.add_parser("close", help="Shut down the optional daemon; persisted state remains.")
    p_close.add_argument("--workdir", help="Workdir for .cassia/ state.")
    p_close.set_defaults(func=cmd_close)

    return agent
