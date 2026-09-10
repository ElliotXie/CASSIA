"""Subclustering utilities for the CASSIA CLI."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import pandas as pd

from .backends import AgentCLIBackend, is_agent_backend
from .runner import MarkerCluster, extract_json_object, load_marker_clusters, utc_now

try:
    from CASSIA.reports.generate_reports import generate_subclustering_report
except Exception:
    generate_subclustering_report = None

try:
    from CASSIA.agents.subclustering.subclustering import build_subcluster_reference_context
except Exception:
    build_subcluster_reference_context = None


def default_subcluster_dir() -> Path:
    """Return the default output directory for a subclustering CLI run."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("cassia_runs") / f"subcluster_{stamp}"


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _as_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    if isinstance(value, tuple):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.replace(";", ",").split(",") if item.strip()]


def _stringify_markers(value: Any) -> str:
    values = _as_list(value)
    return ", ".join(values) if values else str(value or "").strip()


def build_subcluster_prompt(
    clusters: Sequence[MarkerCluster],
    major_cluster_info: str,
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    additional_context: Optional[str] = None,
) -> str:
    """Build the prompt for agent-native subcluster annotation."""
    context_lines = [f"Parent cluster context: {major_cluster_info}"]
    if tissue:
        context_lines.append(f"Tissue: {tissue}")
    if species:
        context_lines.append(f"Species: {species}")
    if additional_context:
        context_lines.append(f"Additional context: {additional_context}")

    cluster_lines = []
    for cluster in clusters:
        marker_csv = ", ".join(cluster.markers)
        cluster_lines.append(f"Subcluster {cluster.cluster_id}: {marker_csv}")

    return f"""You are a careful computational biologist specializing in single-cell subtype and cell-state annotation.

Annotate subclusters from one parent cluster. Treat the parent cluster context as a constraint, but call out contamination or a distinct lineage when the marker evidence is strong.

{chr(10).join(context_lines)}

Marker genes are ranked from strongest to weakest within each subcluster.

{chr(10).join(cluster_lines)}

For every provided subcluster:
1. Identify decisive positive markers and any markers that argue against the parent identity.
2. Assign a broad cell type or state in "main_cell_type".
3. Assign a more specific subtype/state in "sub_cell_type". If the same label is most appropriate, repeat it.
4. Put only marker genes from the provided lists in "key_markers".
5. Keep "reason" concise and evidence-based.

Return only one valid JSON object with this exact schema:
{{
  "parent_cluster": "{major_cluster_info}",
  "subclusters": [
    {{
      "cluster_id": "exact subcluster id from the input",
      "main_cell_type": "broad cell type or state",
      "sub_cell_type": "specific subtype or state",
      "key_markers": ["GENE1", "GENE2"],
      "reason": "concise marker-based rationale"
    }}
  ]
}}

Rules:
- Include every subcluster exactly once.
- Preserve the exact input subcluster IDs.
- Do not include markdown fences, commentary, or text outside the JSON object.

Anti-hallucination rules (strict — applies to "key_markers", "sub_cell_type", and "reason"):
- Every gene symbol you mention MUST be present in that subcluster's provided marker list. Treat the marker list as the only allowed gene vocabulary for that subcluster.
- Reference documents describe canonical marker programs (e.g., "TOX/CXCL13 exhausted CD8", "SPP1/APOE TAM"). They are templates, not facts about this dataset. If the canonical marker for a state is not in the input list, do NOT name that marker — neither in the label nor in the rationale. Use only the markers actually provided to justify the state.
- Do not extrapolate from a partial match. If you see only some markers of a canonical program, state which markers ARE present and qualify the label accordingly (e.g., "exhaustion-leaning CD8 (LAYN+, CTLA4+)" rather than introducing absent companion genes).
- If you must compare against a canonical program from memory or reference docs, restrict the comparison to genes that appear in the input.
"""


def _extract_subcluster_items(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    items = payload.get("subclusters")
    if items is None:
        items = payload.get("clusters")
    if not isinstance(items, list):
        raise ValueError("Subcluster JSON must contain a 'subclusters' list")
    return [item for item in items if isinstance(item, dict)]


def normalize_subcluster_result(payload: Dict[str, Any], clusters: Sequence[MarkerCluster]) -> pd.DataFrame:
    """Normalize agent JSON into the existing CASSIA subclustering CSV format."""
    items = _extract_subcluster_items(payload)
    expected_ids = [str(cluster.cluster_id) for cluster in clusters]
    marker_map = {str(cluster.cluster_id): ", ".join(cluster.markers) for cluster in clusters}

    rows: List[Dict[str, str]] = []
    for idx, item in enumerate(items):
        top_types = _as_list(item.get("most_likely_top2_cell_types") or item.get("top2_cell_types"))
        cluster_id = str(
            item.get("cluster_id")
            or item.get("id")
            or item.get("Result ID")
            or (expected_ids[idx] if idx < len(expected_ids) else "")
        ).strip()
        main_cell_type = str(
            item.get("main_cell_type")
            or item.get("celltype1")
            or item.get("cell_type_1")
            or item.get("primary_cell_type")
            or (top_types[0] if top_types else "")
        ).strip()
        sub_cell_type = str(
            item.get("sub_cell_type")
            or item.get("subtype")
            or item.get("celltype2")
            or item.get("cell_type_2")
            or item.get("secondary_cell_type")
            or (top_types[1] if len(top_types) > 1 else main_cell_type)
        ).strip()
        reason = str(item.get("reason") or item.get("evidence") or item.get("explanation") or "").strip()
        key_markers = _stringify_markers(item.get("key_markers") or item.get("markers"))
        if not key_markers:
            key_markers = marker_map.get(cluster_id, "")

        rows.append({
            "Result ID": cluster_id,
            "main_cell_type": main_cell_type,
            "sub_cell_type": sub_cell_type,
            "key_markers": key_markers,
            "reason": reason,
        })

    returned_ids = [row["Result ID"] for row in rows]
    expected_set = set(expected_ids)
    returned_set = set(returned_ids)
    if expected_set != returned_set:
        if len(rows) == len(expected_ids):
            for row, expected_id in zip(rows, expected_ids):
                row["Result ID"] = expected_id
                if not row["key_markers"]:
                    row["key_markers"] = marker_map.get(expected_id, "")
        else:
            missing = sorted(expected_set - returned_set)
            unexpected = sorted(returned_set - expected_set)
            raise ValueError(
                "Subcluster JSON did not match input cluster IDs. "
                f"Missing: {missing or 'none'}; unexpected: {unexpected or 'none'}"
            )

    row_map = {row["Result ID"]: row for row in rows}
    ordered_rows = []
    for cluster_id in expected_ids:
        row = row_map[cluster_id]
        if not row["main_cell_type"]:
            raise ValueError(f"Subcluster '{cluster_id}' is missing main_cell_type")
        if not row["sub_cell_type"]:
            row["sub_cell_type"] = row["main_cell_type"]
        if not row["key_markers"]:
            row["key_markers"] = marker_map.get(cluster_id, "")
        ordered_rows.append(row)

    return pd.DataFrame(ordered_rows, columns=["Result ID", "main_cell_type", "sub_cell_type", "key_markers", "reason"])


def _build_inline_reference_pointer(args: Any) -> str:
    """Build an inline reference-library pointer for agent-CLI backends.

    Agent CLIs (cursor-agent / codex-cli / claude-cli) have file-reading tools,
    so instead of pre-generating a reference brief via separate LLM calls, we
    just tell the agent where the reference library lives and let it read what
    it needs while annotating. Collapses the 2-step plan + synthesize flow into
    a single annotation call.
    """
    try:
        from CASSIA.agents.reference_agent.utils import get_references_dir
    except ImportError:
        return ""

    references_dir = get_references_dir()
    if not references_dir.exists():
        return ""

    top_level = sorted(
        entry.name for entry in references_dir.iterdir()
        if entry.is_dir() and not entry.name.startswith("_")
    )
    hint = args.reference_cell_type_hint or args.major_cluster_info or ""
    return (
        "<reference_library>\n"
        "A curated subtype reference library is available on the local filesystem "
        f"at: {references_dir}\n\n"
        f"Available top-level lineages: {', '.join(top_level) if top_level else '(none)'}\n\n"
        f"Parent-lineage hint: {hint}\n\n"
        "Use your file-reading tools to consult this library before annotating. "
        "Suggested workflow:\n"
        "  1. List the relevant top-level lineage folder to see what is covered.\n"
        "  2. Read its _overview.md first to get the routing for that lineage, "
        "then one or two specific subtype docs that match the marker pattern.\n"
        "  3. Use the literature evidence to discriminate between similar "
        "subtype programs (for example IFN-gamma response vs type-I IFN, "
        "lipid-phagolysosomal TAM vs resident-like TAM).\n\n"
        "Do not dump the reference contents back into your output. Use them "
        "only as evidence to inform main_cell_type, sub_cell_type, key_markers, "
        "and reason for each subcluster.\n"
        "</reference_library>"
    )


def _build_boost_query_pointer(args: Any, clusters: Sequence[MarkerCluster]) -> str:
    """Emit a prompt block telling an agent how to reverse-query the full DE table.

    Top-30 prompts hide markers that sit at lower ranks but are still
    diagnostic (e.g. SPP1 at rank 87 in SPP1AREGMac). When the user passes
    ``--full-markers``, the agent gets a path to a cassia-boost-compatible
    DE CSV plus the exact ``cassia boost query`` command, so it can confirm
    or rule out lower-ranked markers without us pre-dumping the whole table
    into the prompt.

    Off-loads the work onto a clean CLI rather than asking the agent to
    write its own pandas/grep ad-hoc.
    """
    full_path = getattr(args, "full_markers", None)
    if not full_path:
        return ""
    abs_path = Path(full_path).resolve()
    if not abs_path.exists():
        return ""

    cluster_ids = ", ".join(str(cluster.cluster_id) for cluster in clusters[:8])
    if len(clusters) > 8:
        cluster_ids += ", ..."

    return (
        "<full_marker_query>\n"
        "A full positive differential-expression table is available at:\n"
        f"  {abs_path}\n\n"
        f"Cluster column uses the same subcluster IDs as the prompt above (e.g. {cluster_ids}).\n\n"
        "If a reference document or a candidate subtype name implies a marker "
        "that you do not see in this case's top-N marker list, you can confirm "
        "or rule out that marker by querying the full table via the cassia CLI:\n\n"
        "  cassia boost query --markers '" f"{abs_path}" "' \\\n"
        "      --cluster <case_id> --genes GENE1,GENE2 --format json\n\n"
        "Each result includes avg_log2FC, pct.1, pct.2 and p_val_adj. Treat a "
        "marker as supportive only when avg_log2FC > 0 and p_val_adj < 0.05.\n\n"
        "Use this sparingly — only for decisive markers a reference doc names "
        "that would change the subtype call (for example: confirming SPP1 in an "
        "SPP1AREGMac candidate, or CXCR4 vs CXCR6 to distinguish migratory vs "
        "tissue-resident states). Do not query routine markers that are already "
        "in the top-N list.\n\n"
        "After any boost query, record what you looked up and what you found in "
        "the per-cluster ``reason`` field of the output JSON.\n"
        "</full_marker_query>"
    )


def _build_reference_context(args: Any, clusters: Sequence[MarkerCluster]) -> str:
    backend = getattr(args, "backend", None)
    blocks: List[str] = []

    use_reference = getattr(args, "use_reference", False)
    if use_reference:
        # Agent-CLI path: skip the multi-step plan + synthesize brief LLM calls
        # entirely. The coding agent has filesystem access, so the annotation
        # call itself reads whatever references it needs. Whole reference +
        # annotation collapses to a single agent CLI invocation.
        if is_agent_backend(backend or ""):
            ref_block = _build_inline_reference_pointer(args)
            if ref_block:
                blocks.append(ref_block)
        else:
            ref_block = _build_legacy_reference_context(args, clusters)
            if ref_block:
                blocks.append(ref_block)

    # Annotation-boost mode: agent can reverse-query a full DE table via
    # ``cassia boost query``. Independent of --use-reference: usable alone
    # to give the agent the option to dig past top-N markers, or alongside
    # --use-reference so the agent can verify any reference-named marker.
    if getattr(args, "full_markers", None) and is_agent_backend(backend or ""):
        boost_block = _build_boost_query_pointer(args, clusters)
        if boost_block:
            blocks.append(boost_block)

    return "\n\n".join(blocks)


def _build_legacy_reference_context(args: Any, clusters: Sequence[MarkerCluster]) -> str:
    """Original API-provider reference brief (2 LLM calls). Used only when
    --backend is an HTTP API provider. Agent CLIs take the inline path."""

    if build_subcluster_reference_context is None:
        raise RuntimeError("Subcluster reference retrieval is not available in this installation")

    marker_df = pd.DataFrame({
        "cluster": [cluster.cluster_id for cluster in clusters],
        "markers": [", ".join(cluster.markers) for cluster in clusters],
    })
    context, info = build_subcluster_reference_context(
        marker=marker_df,
        major_cluster_info=args.major_cluster_info,
        provider=args.reference_provider or "openrouter",
        n_genes=args.n_genes,
        tissue=args.tissue,
        species=args.species,
        reference_provider=args.reference_provider,
        reference_model=args.reference_model,
        reference_cell_type_hint=args.reference_cell_type_hint,
        reference_depth=args.reference_depth,
        reference_max_content_length=args.reference_max_content_length,
        reference_max_context_length=args.reference_max_context_length,
    )
    if not info.get("reference_used"):
        return ""
    return context


def run_subcluster(args: Any) -> int:
    """Run one agent-native subcluster annotation pass."""
    out_dir = Path(args.out) if args.out else default_subcluster_dir()
    prompts_dir = out_dir / "prompts"
    raw_dir = out_dir / "raw"
    for directory in (prompts_dir, raw_dir):
        directory.mkdir(parents=True, exist_ok=True)

    clusters = load_marker_clusters(
        input_path=Path(args.markers),
        n_genes=args.n_genes,
        celltype_column=args.cluster_column,
        gene_column=args.gene_column,
        ranking_method=args.ranking_method,
        ascending=args.ascending,
        limit=args.limit,
    )
    if not clusters:
        raise ValueError("No subclusters were loaded from the marker table")

    additional_context = args.additional_context or ""
    reference_context = _build_reference_context(args, clusters)
    if reference_context:
        additional_context = f"{additional_context}\n\n{reference_context}".strip()

    prompt = build_subcluster_prompt(
        clusters=clusters,
        major_cluster_info=args.major_cluster_info,
        tissue=args.tissue,
        species=args.species,
        additional_context=additional_context or None,
    )
    prompt_path = prompts_dir / "subcluster_prompt.md"
    prompt_path.write_text(prompt, encoding="utf-8")

    manifest: Dict[str, Any] = {
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "status": "running",
        "backend": args.backend,
        "marker_table": str(Path(args.markers).resolve()),
        "out_dir": str(out_dir.resolve()),
        "major_cluster_info": args.major_cluster_info,
        "clusters": [{"cluster_id": cluster.cluster_id, "markers": cluster.markers} for cluster in clusters],
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items() if key != "func"},
    }
    _write_json(out_dir / "subcluster_manifest.json", manifest)

    if args.dry_run:
        manifest["status"] = "dry-run"
        manifest["updated_at"] = utc_now()
        manifest["prompt"] = str(prompt_path)
        _write_json(out_dir / "subcluster_manifest.json", manifest)
        print(f"Wrote {out_dir / 'subcluster_manifest.json'}")
        print(f"Wrote {prompt_path}")
        return 0

    raw_path = raw_dir / "subcluster_response.txt"
    backend = AgentCLIBackend(
        args.backend,
        command_template=args.command_template,
        timeout_seconds=args.timeout,
        model=getattr(args, "model", None),
        reasoning_effort=getattr(args, "reasoning_effort", None),
    )
    output_text = backend.run(
        prompt=prompt,
        prompt_file=prompt_path,
        cwd=Path.cwd(),
        context={
            "input": str(Path(args.markers)),
            "out": str(out_dir),
            "cluster": args.major_cluster_info,
            "agent_output_file": str(raw_path.resolve()),
        },
    )
    raw_path.write_text(output_text + ("\n" if not output_text.endswith("\n") else ""), encoding="utf-8")

    payload = extract_json_object(output_text)
    result_df = normalize_subcluster_result(payload, clusters)
    csv_path = out_dir / "subcluster_results.csv"
    result_df.to_csv(csv_path, index=False)

    html_path = out_dir / "subcluster_report.html"
    if generate_subclustering_report is not None:
        generate_subclustering_report(
            str(csv_path),
            html_report_path=str(html_path),
            model_name=args.backend,
        )

    manifest.update({
        "status": "completed",
        "updated_at": utc_now(),
        "prompt": str(prompt_path),
        "raw_response": str(raw_path),
        "result_csv": str(csv_path),
        "html_report": str(html_path) if html_path.exists() else "",
    })
    _write_json(out_dir / "subcluster_manifest.json", manifest)
    print(f"Wrote {out_dir / 'subcluster_manifest.json'}")
    print(f"Wrote {csv_path}")
    if html_path.exists():
        print(f"Wrote {html_path}")
    return 0
