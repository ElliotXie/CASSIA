"""Annotation boost utilities for the CASSIA CLI."""

from __future__ import annotations

import html
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from .backends import AgentCLIBackend
from .runner import extract_json_object, slugify, utc_now

try:
    from CASSIA.agents.annotation_boost.annotation_boost import (
        format_summary_to_html as _format_annotation_boost_summary_to_html,
    )
except Exception:
    _format_annotation_boost_summary_to_html = None


SCANPY_TO_SEURAT_COLUMNS = {
    "group": "cluster",
    "names": "gene",
    "logfoldchanges": "avg_log2FC",
    "pvals": "p_val",
    "pvals_adj": "p_val_adj",
}

DEFAULT_GENE_COLUMNS = ("gene", "Gene", "GENE", "names")
DEFAULT_CLUSTER_COLUMNS = ("cluster", "Cluster ID", "group", "ident", "seurat_clusters")
DISPLAY_COLUMNS = ("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "p_val", "status")


def parse_gene_args(gene_values: Sequence[str]) -> List[str]:
    """Parse comma/space/newline separated gene arguments while preserving order."""
    genes: List[str] = []
    seen = set()
    for value in gene_values:
        for gene in re.split(r"[,\s]+", value):
            cleaned = gene.strip()
            if cleaned and cleaned.upper() not in seen:
                genes.append(cleaned)
                seen.add(cleaned.upper())
    return genes


def load_marker_table(path: Path) -> pd.DataFrame:
    """Load and lightly normalize a raw marker table."""
    df = pd.read_csv(path)
    unnamed_cols = [col for col in df.columns if str(col).startswith("Unnamed:")]
    if unnamed_cols:
        df = df.drop(columns=unnamed_cols)
    return df.rename(columns=SCANPY_TO_SEURAT_COLUMNS)


def _resolve_column(df: pd.DataFrame, requested: Optional[str], candidates: Iterable[str], label: str) -> str:
    if requested:
        if requested not in df.columns:
            raise ValueError(f"{label} column '{requested}' was not found. Available columns: {list(df.columns)}")
        return requested
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    raise ValueError(f"Could not infer {label} column. Available columns: {list(df.columns)}")


def _ordered_columns(df: pd.DataFrame) -> List[str]:
    preferred = [col for col in DISPLAY_COLUMNS if col in df.columns]
    extra = [col for col in df.columns if col not in preferred]
    return preferred + extra


def query_marker_genes(
    marker_path: Path,
    genes: Sequence[str],
    cluster: Optional[str] = None,
    gene_column: Optional[str] = None,
    cluster_column: Optional[str] = None,
) -> pd.DataFrame:
    """Return marker statistics for requested genes, optionally within one cluster."""
    if not genes:
        raise ValueError("At least one gene must be provided")

    df = load_marker_table(marker_path)
    gene_col = _resolve_column(df, gene_column, DEFAULT_GENE_COLUMNS, "gene")
    resolved_cluster_col = None
    work = df.copy()

    if cluster is not None:
        resolved_cluster_col = _resolve_column(work, cluster_column, DEFAULT_CLUSTER_COLUMNS, "cluster")
        work = work[work[resolved_cluster_col].astype(str) == str(cluster)].copy()
        if work.empty:
            raise ValueError(f"No marker rows found for cluster '{cluster}'")
    elif cluster_column:
        resolved_cluster_col = _resolve_column(work, cluster_column, DEFAULT_CLUSTER_COLUMNS, "cluster")
    else:
        for candidate in DEFAULT_CLUSTER_COLUMNS:
            if candidate in work.columns:
                resolved_cluster_col = candidate
                break

    rows = []
    gene_lookup = work[gene_col].astype(str).str.upper()
    for gene in genes:
        matches = work[gene_lookup == gene.upper()].copy()
        if matches.empty:
            rows.append({
                "gene": gene,
                "cluster": cluster if cluster is not None else "",
                "status": "not_found",
            })
            continue

        matches = matches.copy()
        matches["status"] = "found"
        if gene_col != "gene":
            matches["gene"] = matches[gene_col].astype(str)
        if resolved_cluster_col:
            matches["cluster"] = matches[resolved_cluster_col].astype(str)
        elif not resolved_cluster_col and "cluster" not in matches.columns:
            matches["cluster"] = ""
        rows.extend(matches.to_dict(orient="records"))

    result = pd.DataFrame(rows)
    return result[_ordered_columns(result)]


def write_query_output(df: pd.DataFrame, output_format: str = "table", out: Optional[Path] = None) -> str:
    """Serialize query output and optionally write it to disk."""
    if output_format == "json":
        text = json.dumps(df.to_dict(orient="records"), indent=2, ensure_ascii=False)
    elif output_format == "csv":
        text = df.to_csv(index=False)
    else:
        text = df.to_string(index=False)

    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text + ("\n" if not text.endswith("\n") else ""), encoding="utf-8")
    return text


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _read_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def _sort_direction(ranking_method: str, ascending: Optional[bool]) -> bool:
    if ascending is not None:
        return ascending
    return ranking_method in {"p_val_adj", "p_val"}


def get_cluster_top_markers(
    marker_path: Path,
    cluster: str,
    n_genes: int = 50,
    gene_column: Optional[str] = None,
    cluster_column: Optional[str] = None,
    ranking_method: str = "avg_log2FC",
    ascending: Optional[bool] = None,
) -> Tuple[List[str], pd.DataFrame]:
    """Return top marker genes and rows for one cluster from a raw marker table."""
    df = load_marker_table(marker_path)
    gene_col = _resolve_column(df, gene_column, DEFAULT_GENE_COLUMNS, "gene")
    cluster_col = _resolve_column(df, cluster_column, DEFAULT_CLUSTER_COLUMNS, "cluster")
    cluster_df = df[df[cluster_col].astype(str) == str(cluster)].copy()
    if cluster_df.empty:
        raise ValueError(f"No marker rows found for cluster '{cluster}'")

    for column in ("avg_log2FC", "p_val_adj", "p_val", "pct.1", "pct.2", "Score"):
        if column in cluster_df.columns:
            cluster_df[column] = pd.to_numeric(cluster_df[column], errors="coerce")

    filtered_df = cluster_df.copy()
    if "p_val_adj" in filtered_df.columns:
        filtered_df = filtered_df[filtered_df["p_val_adj"] < 0.05]
    if "avg_log2FC" in filtered_df.columns:
        filtered_df = filtered_df[filtered_df["avg_log2FC"] > 0.25]
    if "pct.1" in filtered_df.columns and "pct.2" in filtered_df.columns:
        filtered_df = filtered_df[(filtered_df["pct.1"] >= 0.1) | (filtered_df["pct.2"] >= 0.1)]
    elif "pct.1" in filtered_df.columns:
        filtered_df = filtered_df[filtered_df["pct.1"] >= 0.1]

    if not filtered_df.empty:
        cluster_df = filtered_df

    if ranking_method in cluster_df.columns:
        cluster_df = cluster_df.sort_values(
            ranking_method,
            ascending=_sort_direction(ranking_method, ascending),
            na_position="last",
        )

    top_df = cluster_df.head(n_genes).copy()
    if gene_col != "gene":
        top_df["gene"] = top_df[gene_col].astype(str)
    if cluster_col != "cluster":
        top_df["cluster"] = top_df[cluster_col].astype(str)
    genes = top_df["gene"].astype(str).tolist()
    return genes, top_df[_ordered_columns(top_df)]


def load_annotation_context(run_dir: Path, cluster: str) -> Dict[str, Any]:
    """Load the original CASSIA annotation for a cluster from a run folder."""
    results_path = run_dir / "results.json"
    if results_path.exists():
        results = _read_json(results_path, {})
        lookup = cluster if cluster in results else str(cluster)
        if lookup in results:
            return {
                "source": str(results_path),
                "cluster_id": cluster,
                "annotation": results[lookup].get("analysis_result", results[lookup]),
            }

    summary_candidates = [run_dir / "summary.csv"] + sorted(run_dir.glob("*_summary.csv"))
    for summary_path in summary_candidates:
        if not summary_path.exists():
            continue
        df = pd.read_csv(summary_path)
        cluster_col = None
        for candidate in ("Cluster ID", "True Cell Type", "cluster"):
            if candidate in df.columns:
                cluster_col = candidate
                break
        if not cluster_col:
            continue
        match = df[df[cluster_col].astype(str) == str(cluster)]
        if match.empty:
            continue
        row = match.iloc[0].to_dict()
        return {
            "source": str(summary_path),
            "cluster_id": cluster,
            "annotation": row,
        }

    raise ValueError(f"Could not find annotation for cluster '{cluster}' in {run_dir}")


def format_annotation_context(context: Dict[str, Any]) -> str:
    """Format original annotation context for the boost prompt."""
    annotation = context.get("annotation", {})
    if isinstance(annotation, dict):
        return json.dumps(annotation, indent=2, ensure_ascii=False)
    return str(annotation)


def extract_check_genes(text: str, max_genes: int = 20) -> List[str]:
    """Extract requested genes from <check_genes> tags."""
    blocks = re.findall(r"<check_genes>\s*(.*?)\s*</check_genes>", text, flags=re.DOTALL | re.IGNORECASE)
    genes: List[str] = []
    for block in blocks:
        genes.extend(parse_gene_args([block]))
    return genes[:max_genes]


def build_boost_prompt(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    annotation_context: str,
    strategy: str = "breadth",
    additional_task: Optional[str] = None,
    max_genes_per_round: int = 20,
) -> str:
    """Build the initial annotation boost prompt for an agent CLI backend."""
    strategy_text = (
        "Use a depth-first strategy: investigate one leading hypothesis at a time, "
        "then go deeper into subtype/state if it is supported."
        if strategy == "depth"
        else "Use a breadth-first strategy: consider up to three plausible hypotheses, "
        "then choose decisive positive and negative markers to separate them."
    )
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    return f"""You are CASSIA annotation boost, a careful senior computational biologist called in to stress-test a single-cell annotation.

Cluster: {cluster}
Dataset context: {major_cluster_info}
Top ranked markers from the raw differential expression table:
{", ".join(top_markers)}

Original CASSIA annotation:
{annotation_context}
{task_text}
Workflow:
1. Briefly assess whether the original annotation is robust, ambiguous, mixed, or likely wrong.
2. {strategy_text}
3. If more evidence is needed, request marker checks using this exact tag format:
<check_genes>GENE1,GENE2,GENE3</check_genes>
4. Request no more than {max_genes_per_round} genes per round. Use official gene symbols only.
5. After CASSIA returns marker statistics, refine or pivot. Do not repeat already checked genes unless necessary.
6. When ready, return only one valid JSON object with this exact schema:
{{
  "final_cell_type": "general cell type",
  "final_sub_cell_type": "specific subtype or state if applicable",
  "confidence": "low|medium|high",
  "changed_from_original": false,
  "checked_genes": ["GENE1", "GENE2"],
  "supporting_markers": ["GENE1", "GENE2"],
  "refuting_markers": ["GENE3"],
  "alternatives": ["alternative if confidence is not high"],
  "evidence": "concise evidence based on original markers and queried marker statistics",
  "recommended_next_steps": "optional next validation step"
}}

Start with your assessment and the first <check_genes> request unless the original annotation is already fully supported by the top markers. Do not output final JSON until you have enough evidence.
"""


def build_boost_followup_prompt(transcript: str, query_text: Optional[str], is_final_round: bool) -> str:
    """Build a stateless follow-up prompt containing the transcript so far."""
    final_instruction = (
        "You are at the final round. Return only the final JSON object now."
        if is_final_round
        else "Continue the boost analysis. Either request another <check_genes> list or return the final JSON object."
    )
    query_section = f"\nLatest CASSIA marker query results:\n{query_text}\n" if query_text else ""
    return f"""Continue this CASSIA annotation boost session.

Transcript so far:
{transcript}
{query_section}
{final_instruction}

Remember:
- Use <check_genes>GENE1,GENE2</check_genes> for more local marker checks.
- Final output must be one valid JSON object only, with no markdown fences or commentary.
"""


def _append_transcript(messages: List[Dict[str, str]]) -> str:
    parts = []
    for msg in messages:
        parts.append(f"## {msg['role'].upper()}\n\n{msg['content']}")
    return ("\n\n" + ("=" * 80) + "\n\n").join(parts)


def normalize_boost_result(result: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize final boost JSON and accept common fallback field names."""
    normalized = dict(result)
    if "final_cell_type" not in normalized and "main_cell_type" in normalized:
        normalized["final_cell_type"] = normalized["main_cell_type"]
    if not normalized.get("final_cell_type"):
        raise ValueError("Final boost JSON is missing required field 'final_cell_type'")
    normalized.setdefault("final_sub_cell_type", "")
    normalized.setdefault("confidence", "")
    normalized.setdefault("changed_from_original", None)
    for key in ("checked_genes", "supporting_markers", "refuting_markers", "alternatives"):
        value = normalized.get(key)
        if value is None:
            normalized[key] = []
        elif isinstance(value, str):
            normalized[key] = [item.strip() for item in value.split(",") if item.strip()]
        elif not isinstance(value, list):
            normalized[key] = [str(value)]
    normalized.setdefault("evidence", "")
    normalized.setdefault("recommended_next_steps", "")
    return normalized


def write_boost_report(boost_dir: Path, manifest: Dict[str, Any], result: Optional[Dict[str, Any]], errors: List[str]) -> Path:
    """Write a compact Markdown report for a boost run."""
    lines = [
        "# CASSIA Boost Report",
        "",
        f"- Cluster: `{manifest.get('cluster')}`",
        f"- Backend: `{manifest.get('backend')}`",
        f"- Status: `{manifest.get('status')}`",
        f"- Checked genes: {', '.join(manifest.get('checked_genes', [])) or 'None'}",
        "",
    ]
    if result:
        lines.extend([
            "## Final Annotation",
            "",
            f"- Cell type: {result.get('final_cell_type', '')}",
            f"- Subtype/state: {result.get('final_sub_cell_type', '')}",
            f"- Confidence: {result.get('confidence', '')}",
            f"- Changed from original: {result.get('changed_from_original')}",
            "",
            "## Evidence",
            "",
            result.get("evidence", ""),
            "",
        ])
    if errors:
        lines.extend(["## Errors", ""])
        lines.extend(f"- {error}" for error in errors)
    report_path = boost_dir / "report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def _html_content(value: Any) -> str:
    if value is None or value == "":
        return "No information available"
    if isinstance(value, (dict, list)):
        text = json.dumps(value, indent=2, ensure_ascii=False)
    else:
        text = str(value)
    escaped = html.escape(text)
    return escaped.replace("\n", "<br>")


def _display_list(values: Sequence[Any]) -> str:
    cleaned = [str(value).strip() for value in values if str(value).strip()]
    return ", ".join(cleaned) if cleaned else "None"


def _row_value(row: pd.Series, column: str, default: Optional[str] = "N/A") -> Optional[str]:
    if column not in row or pd.isna(row[column]):
        return default
    return str(row[column])


def _query_frames_gene_stats(query_frames: Sequence[pd.DataFrame]) -> Dict[str, Dict[str, Optional[str]]]:
    gene_stats: Dict[str, Dict[str, Optional[str]]] = {}
    for frame in query_frames:
        if "gene" not in frame.columns:
            continue
        for _, row in frame.iterrows():
            gene = str(row.get("gene", "")).strip()
            if not gene:
                continue
            if str(row.get("status", "")).lower() == "not_found":
                continue
            gene_stats[gene.upper()] = {
                "avg_log2FC": _row_value(row, "avg_log2FC", "N/A"),
                "pct.1": _row_value(row, "pct.1", None),
                "pct.2": _row_value(row, "pct.2", None),
                "p_val_adj": _row_value(row, "p_val_adj", "N/A"),
            }
    return gene_stats


def _query_frame_genes(frame: pd.DataFrame) -> List[str]:
    if "gene" not in frame.columns:
        return []
    genes: List[str] = []
    seen = set()
    for gene in frame["gene"].dropna().astype(str):
        key = gene.upper()
        if key not in seen:
            genes.append(gene)
            seen.add(key)
    return genes


def _query_findings_html(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No marker query results returned."
    lines: List[str] = []
    for _, row in frame.iterrows():
        gene = str(row.get("gene", "")).strip() or "Unknown gene"
        if str(row.get("status", "")).lower() == "not_found":
            lines.append(f"{gene}: not found in the selected cluster marker table.")
            continue
        stats = []
        for column in ("avg_log2FC", "pct.1", "pct.2", "p_val_adj", "p_val"):
            if column in row and not pd.isna(row[column]):
                stats.append(f"{column}={row[column]}")
        lines.append(f"{gene}: {', '.join(stats) if stats else 'found in marker table'}")
    return "<br>".join(html.escape(line) for line in lines)


def _assistant_rounds(messages: Sequence[Dict[str, str]]) -> List[str]:
    return [msg.get("content", "") for msg in messages if msg.get("role") == "assistant"]


def _strip_check_gene_tags(text: str) -> str:
    cleaned = re.sub(
        r"<check_genes>\s*[\s\S]*?\s*</check_genes>",
        "",
        text,
        flags=re.IGNORECASE,
    )
    return cleaned.strip()


def _top_marker_gene_list(top_marker_rows: pd.DataFrame, limit: int = 25) -> List[str]:
    if "gene" not in top_marker_rows.columns:
        return []
    return top_marker_rows["gene"].dropna().astype(str).head(limit).tolist()


def build_boost_summary_tags(
    manifest: Dict[str, Any],
    messages: Sequence[Dict[str, str]],
    result: Optional[Dict[str, Any]],
    errors: Sequence[str],
    annotation_context: Dict[str, Any],
    query_frames: Sequence[pd.DataFrame],
    top_marker_rows: pd.DataFrame,
) -> str:
    """Build the tagged summary consumed by the original annotation boost HTML template."""
    checked_genes = manifest.get("checked_genes", [])
    cluster = manifest.get("cluster", "")
    overview = (
        f"Cluster {html.escape(str(cluster))} was reviewed with CASSIA annotation boost. "
        f"The CLI agent completed {len(query_frames)} marker-query round(s) and checked "
        f"{len(checked_genes)} unique gene(s)."
    )
    if result:
        overview += f" Final confidence: {html.escape(str(result.get('confidence', '')))}."
    elif errors:
        overview += " The boost run did not produce a valid final JSON result."

    initial_assessment = (
        f"Original annotation source: {html.escape(str(annotation_context.get('source', 'unknown')))}<br>"
        f"Original annotation:<br>{_html_content(annotation_context.get('annotation'))}<br>"
        f"Top ranked markers: {_html_content(_display_list(_top_marker_gene_list(top_marker_rows)))}"
    )

    assistant_messages = _assistant_rounds(messages)
    iteration_blocks: List[str] = []
    for idx, frame in enumerate(query_frames, start=1):
        raw_hypotheses = assistant_messages[idx - 1] if idx - 1 < len(assistant_messages) else ""
        hypotheses = _strip_check_gene_tags(raw_hypotheses) or "The agent requested local marker evidence for this round."
        genes_checked = _display_list(_query_frame_genes(frame))
        findings = _query_findings_html(frame)
        iteration_blocks.append(
            f"""<ITERATION_{idx}>
<HYPOTHESES>
{_html_content(hypotheses)}
</HYPOTHESES>
<GENES_CHECKED>
{genes_checked}
</GENES_CHECKED>
<KEY_FINDINGS>
{findings}
</KEY_FINDINGS>
</ITERATION_{idx}>"""
        )

    if result:
        final_lines = [
            f"Cell type: {result.get('final_cell_type', '')}",
            f"Subtype/state: {result.get('final_sub_cell_type', '')}",
            f"Confidence: {result.get('confidence', '')}",
            f"Changed from original: {result.get('changed_from_original')}",
            f"Evidence: {result.get('evidence', '')}",
        ]
        final_annotation = "<br>".join(html.escape(line) for line in final_lines if line)
        marker_summary = "<br>".join([
            f"Supporting markers: {_html_content(_display_list(result.get('supporting_markers', [])))}",
            f"Refuting markers: {_html_content(_display_list(result.get('refuting_markers', [])))}",
            f"All checked markers: {_html_content(_display_list(result.get('checked_genes', checked_genes)))}",
            f"Top ranked raw markers: {_html_content(_display_list(_top_marker_gene_list(top_marker_rows, limit=15)))}",
        ])
        recommendations = _html_content(result.get("recommended_next_steps") or "No additional validation step was requested by the boost agent.")
    else:
        final_annotation = "No final JSON annotation was produced."
        marker_summary = f"Checked markers before failure: {_html_content(_display_list(checked_genes))}"
        recommendations = _html_content("; ".join(errors) if errors else "Retry with one additional iteration or inspect the raw transcript.")

    return "\n".join([
        "<OVERVIEW>",
        overview,
        "</OVERVIEW>",
        "<INITIAL_ASSESSMENT>",
        initial_assessment,
        "</INITIAL_ASSESSMENT>",
        *iteration_blocks,
        "<FINAL_ANNOTATION>",
        final_annotation,
        "</FINAL_ANNOTATION>",
        "<MARKER_SUMMARY>",
        marker_summary,
        "</MARKER_SUMMARY>",
        "<RECOMMENDATIONS>",
        recommendations,
        "</RECOMMENDATIONS>",
    ])


def write_boost_html_report(
    boost_dir: Path,
    manifest: Dict[str, Any],
    messages: Sequence[Dict[str, str]],
    result: Optional[Dict[str, Any]],
    errors: Sequence[str],
    annotation_context: Dict[str, Any],
    query_frames: Sequence[pd.DataFrame],
    top_marker_rows: pd.DataFrame,
    strategy: str,
) -> Path:
    """Write the original annotation boost HTML report from CLI boost artifacts."""
    summary_text = build_boost_summary_tags(
        manifest=manifest,
        messages=messages,
        result=result,
        errors=errors,
        annotation_context=annotation_context,
        query_frames=query_frames,
        top_marker_rows=top_marker_rows,
    )
    tags_path = boost_dir / "summary_tags.txt"
    tags_path.write_text(summary_text + "\n", encoding="utf-8")

    html_path = boost_dir / "summary.html"
    gene_stats = _query_frames_gene_stats(query_frames)
    if _format_annotation_boost_summary_to_html is not None:
        returned_path = _format_annotation_boost_summary_to_html(
            summary_text,
            str(html_path),
            search_strategy=strategy,
            report_style="per_iteration",
            gene_stats=gene_stats,
        )
        return Path(returned_path)

    html_path.write_text(
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>CASSIA Boost Summary</title></head>"
        f"<body><pre>{html.escape(summary_text)}</pre></body></html>\n",
        encoding="utf-8",
    )
    return html_path


def run_boost(args: Any) -> int:
    """Run the CLI annotation boost agent loop for one cluster."""
    run_dir = Path(args.run)
    cluster_slug = slugify(args.cluster)
    boost_dir = Path(args.out) if args.out else run_dir / "boost" / cluster_slug
    prompts_dir = boost_dir / "prompts"
    raw_dir = boost_dir / "raw"
    queries_dir = boost_dir / "queries"
    for directory in (prompts_dir, raw_dir, queries_dir):
        directory.mkdir(parents=True, exist_ok=True)

    annotation_context = load_annotation_context(run_dir, args.cluster)
    top_markers, top_marker_rows = get_cluster_top_markers(
        marker_path=Path(args.markers),
        cluster=args.cluster,
        n_genes=args.n_genes,
        gene_column=args.gene_column,
        cluster_column=args.cluster_column,
        ranking_method=args.ranking_method,
        ascending=args.ascending,
    )
    top_marker_rows.to_csv(boost_dir / "top_markers.csv", index=False)

    manifest: Dict[str, Any] = {
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "status": "running",
        "cluster": args.cluster,
        "backend": args.backend,
        "run_dir": str(run_dir.resolve()),
        "boost_dir": str(boost_dir.resolve()),
        "marker_table": str(Path(args.markers).resolve()),
        "annotation_source": annotation_context.get("source"),
        "checked_genes": [],
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items() if key != "func"},
    }
    _write_json(boost_dir / "boost_manifest.json", manifest)

    prompt = build_boost_prompt(
        cluster=args.cluster,
        major_cluster_info=args.major_cluster_info,
        top_markers=top_markers,
        annotation_context=format_annotation_context(annotation_context),
        strategy=args.strategy,
        additional_task=args.additional_task,
        max_genes_per_round=args.max_genes_per_round,
    )

    if args.dry_run:
        (prompts_dir / "round_001.md").write_text(prompt, encoding="utf-8")
        manifest["status"] = "dry-run"
        manifest["updated_at"] = utc_now()
        _write_json(boost_dir / "boost_manifest.json", manifest)
        return 0

    backend = AgentCLIBackend(args.backend, command_template=args.command_template, timeout_seconds=args.timeout)
    messages: List[Dict[str, str]] = []
    errors: List[str] = []
    final_result: Optional[Dict[str, Any]] = None
    latest_query_text: Optional[str] = None
    query_frames: List[pd.DataFrame] = []
    checked_seen = set()

    for round_idx in range(1, args.iterations + 1):
        if round_idx > 1:
            prompt = build_boost_followup_prompt(
                transcript=_append_transcript(messages),
                query_text=latest_query_text,
                is_final_round=round_idx == args.iterations,
            )

        prompt_path = prompts_dir / f"round_{round_idx:03d}.md"
        raw_path = raw_dir / f"round_{round_idx:03d}.txt"
        prompt_path.write_text(prompt, encoding="utf-8")
        messages.append({"role": "user", "content": prompt})

        try:
            raw = backend.run(
                prompt,
                prompt_path,
                boost_dir,
                {
                    "input": str(Path(args.markers).resolve()),
                    "out": str(boost_dir.resolve()),
                    "cluster": args.cluster,
                    "agent_output_file": str(raw_path),
                },
            )
            raw_path.write_text(raw, encoding="utf-8")
            messages.append({"role": "assistant", "content": raw})

            try:
                parsed = extract_json_object(raw)
                final_result = normalize_boost_result(parsed)
                break
            except Exception:
                final_result = None

            genes = extract_check_genes(raw, max_genes=args.max_genes_per_round)
            genes = [gene for gene in genes if gene.upper() not in checked_seen]
            if not genes:
                latest_query_text = "No new <check_genes> request was found. Please return final JSON or request new genes."
                continue

            for gene in genes:
                checked_seen.add(gene.upper())
            query_df = query_marker_genes(
                marker_path=Path(args.markers),
                genes=genes,
                cluster=args.cluster,
                gene_column=args.gene_column,
                cluster_column=args.cluster_column,
            )
            query_frames.append(query_df.copy())
            query_df.to_csv(queries_dir / f"round_{round_idx:03d}.csv", index=False)
            (queries_dir / f"round_{round_idx:03d}.json").write_text(
                json.dumps(query_df.to_dict(orient="records"), indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            latest_query_text = query_df.to_string(index=False)
            manifest["checked_genes"] = list(checked_seen)
            manifest["updated_at"] = utc_now()
            _write_json(boost_dir / "boost_manifest.json", manifest)
        except Exception as exc:
            errors.append(str(exc))
            break

    if final_result is None and not errors:
        final_prompt = build_boost_followup_prompt(
            transcript=_append_transcript(messages),
            query_text=latest_query_text,
            is_final_round=True,
        )
        prompt_path = prompts_dir / "final.md"
        raw_path = raw_dir / "final.txt"
        prompt_path.write_text(final_prompt, encoding="utf-8")
        messages.append({"role": "user", "content": final_prompt})
        try:
            raw = backend.run(
                final_prompt,
                prompt_path,
                boost_dir,
                {
                    "input": str(Path(args.markers).resolve()),
                    "out": str(boost_dir.resolve()),
                    "cluster": args.cluster,
                    "agent_output_file": str(raw_path),
                },
            )
            raw_path.write_text(raw, encoding="utf-8")
            messages.append({"role": "assistant", "content": raw})
            final_result = normalize_boost_result(extract_json_object(raw))
        except Exception as exc:
            errors.append(str(exc))

    transcript = _append_transcript(messages)
    (boost_dir / "transcript.md").write_text(transcript + "\n", encoding="utf-8")
    if final_result:
        _write_json(boost_dir / "final.json", final_result)

    manifest["status"] = "completed" if final_result and not errors else "failed"
    manifest["updated_at"] = utc_now()
    manifest["checked_genes"] = list(checked_seen)
    manifest["final_json"] = str(boost_dir / "final.json") if final_result else None
    manifest["errors"] = errors
    report_path = write_boost_report(boost_dir, manifest, final_result, errors)
    html_path = write_boost_html_report(
        boost_dir=boost_dir,
        manifest=manifest,
        messages=messages,
        result=final_result,
        errors=errors,
        annotation_context=annotation_context,
        query_frames=query_frames,
        top_marker_rows=top_marker_rows,
        strategy=args.strategy,
    )
    manifest["markdown_report"] = str(report_path)
    manifest["html_report"] = str(html_path)
    manifest["summary_tags"] = str(boost_dir / "summary_tags.txt")
    _write_json(boost_dir / "boost_manifest.json", manifest)
    print(f"Wrote {boost_dir / 'boost_manifest.json'}")
    print(f"Wrote {boost_dir / 'transcript.md'}")
    print(f"Wrote {report_path}")
    print(f"Wrote {html_path}")
    return 0 if final_result and not errors else 1
