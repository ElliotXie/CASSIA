"""Run-folder orchestration for the CASSIA CLI."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd

from CASSIA import __version__
from CASSIA.core.marker_utils import get_top_markers, split_markers

from .backends import AgentCLIBackend, is_api_backend


@dataclass
class MarkerCluster:
    """A marker list prepared for one cluster."""

    cluster_id: str
    markers: List[str]


def utc_now() -> str:
    """Return an ISO timestamp with second precision."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def default_run_dir() -> Path:
    """Return the default run directory for a new CLI run."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("cassia_runs") / stamp


def slugify(value: str) -> str:
    """Make a stable filename fragment for a cluster ID."""
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return slug.strip("._") or "cluster"


def ensure_run_dirs(run_dir: Path) -> None:
    """Create the standard CASSIA CLI run-folder layout."""
    for child in ("prompts", "raw", "errors"):
        (run_dir / child).mkdir(parents=True, exist_ok=True)


def _looks_preformatted(df: pd.DataFrame, gene_column: Optional[str]) -> bool:
    if len(df.columns) == 2:
        return True
    if gene_column and gene_column in df.columns:
        sample = df[gene_column].dropna().astype(str).head(5)
        return sample.str.contains(r"[,; ]").any()
    return False


def load_marker_clusters(
    input_path: Path,
    n_genes: int = 50,
    celltype_column: Optional[str] = None,
    gene_column: Optional[str] = None,
    ranking_method: str = "avg_log2FC",
    ascending: Optional[bool] = None,
    limit: Optional[int] = None,
) -> List[MarkerCluster]:
    """Load marker data and return one marker list per cluster."""
    df = pd.read_csv(input_path)
    unnamed_cols = [col for col in df.columns if str(col).startswith("Unnamed:")]
    if unnamed_cols:
        df = df.drop(columns=unnamed_cols)

    if len(df.columns) < 2:
        raise ValueError("Marker input must have at least two columns")

    if _looks_preformatted(df, gene_column):
        prepared = df.copy()
    else:
        prepared = get_top_markers(
            df,
            n_genes=n_genes,
            ranking_method=ranking_method,
            ascending=ascending,
        )

    cluster_col = celltype_column or prepared.columns[0]
    marker_col = gene_column or prepared.columns[1]

    if cluster_col not in prepared.columns:
        raise ValueError(f"Cluster column '{cluster_col}' was not found")
    if marker_col not in prepared.columns:
        raise ValueError(f"Gene/marker column '{marker_col}' was not found")

    clusters: List[MarkerCluster] = []
    for _, row in prepared.iterrows():
        cluster_id = str(row[cluster_col])
        markers = split_markers(str(row[marker_col]))
        clusters.append(MarkerCluster(cluster_id=cluster_id, markers=markers[:n_genes]))

    if limit is not None:
        clusters = clusters[:limit]

    return clusters


def build_annotation_prompt(
    cluster: MarkerCluster,
    tissue: str,
    species: str,
    additional_info: Optional[str] = None,
) -> str:
    """Create the prompt sent to local agent CLI backends."""
    marker_list = ", ".join(cluster.markers)
    marker_lines = "\n".join(f"{idx + 1}. {marker}" for idx, marker in enumerate(cluster.markers))
    extra = additional_info.strip() if additional_info else "None"
    is_tissue_blind = tissue.lower() in ["none", "tissue blind"] if tissue else True
    tissue_instruction = (
        "The tissue of origin is not specified. Consider multiple plausible tissues "
        "and include a ranked possible_tissues list."
        if is_tissue_blind
        else f"The tissue of origin is {tissue}. Prioritize annotations consistent with this tissue."
    )
    possible_tissues_schema = (
        '  "possible_tissues": ["ranked possible tissue 1", "ranked possible tissue 2"],\n'
        if is_tissue_blind
        else ""
    )

    return f"""You are CASSIA, a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).

Annotate one cluster from a {species} dataset.
{tissue_instruction}
Cluster ID: {cluster.cluster_id}
Additional context: {extra}

The marker genes are ranked by expression intensity from high to low.
Ranked marker genes as a comma-separated list:
{marker_list}

Ranked marker genes with rank numbers:
{marker_lines}

Analyze the marker list systematically before writing the final JSON:
1. Identify key functional or pathway markers and what they imply.
2. Identify key cell type markers and the cell types they support.
3. Cross-check the marker pattern against established scRNA-seq knowledge, marker databases, and relevant literature you know.
4. Determine the most probable general cell type.
5. Rank the top three most probable sub-cell types or states from most likely to least likely.
6. Consider mixed-cell evidence only when multiple distinct cell types are strongly supported by several high-ranking markers.
7. Prefer specific, biologically standard labels over vague labels, but do not over-specify when marker evidence is weak.

Return only one valid JSON object with this exact schema:
{{
  "main_cell_type": "general cell type",
  "sub_cell_types": ["most likely subtype/state", "second most likely", "third most likely"],
  "possible_mixed_cell_types": [],
{possible_tissues_schema}  "key_functional_markers": [
    {{"genes": ["GENE1", "GENE2"], "interpretation": "brief functional/pathway interpretation"}}
  ],
  "key_cell_type_markers": [
    {{"genes": ["GENE3", "GENE4"], "supports": "cell type or state", "interpretation": "brief marker interpretation"}}
  ],
  "confidence": "low|medium|high",
  "evidence": "concise marker-based rationale without step-by-step hidden reasoning"
}}

Rules:
- Do not include markdown fences, commentary, tool logs, or any text outside the JSON object.
- Use only marker genes that are present in the provided ranked marker list when citing evidence.
- If the main cell type and most likely subtype are the same, still include it in sub_cell_types.
- Keep evidence concise and defensible for an expert user.
"""


def extract_json_object(text: str) -> Dict[str, Any]:
    """Extract the first valid JSON object from agent output."""
    stripped = text.strip()
    if not stripped:
        raise ValueError("Agent returned an empty response")

    try:
        value = json.loads(stripped)
        if isinstance(value, dict):
            return value
    except json.JSONDecodeError:
        pass

    fence_matches = re.findall(r"```(?:json)?\s*(\{.*?\})\s*```", stripped, flags=re.DOTALL)
    for candidate in fence_matches:
        try:
            value = json.loads(candidate)
            if isinstance(value, dict):
                return value
        except json.JSONDecodeError:
            continue

    starts = [match.start() for match in re.finditer(r"\{", stripped)]
    for start in starts:
        depth = 0
        in_string = False
        escape = False
        for idx in range(start, len(stripped)):
            char = stripped[idx]
            if in_string:
                if escape:
                    escape = False
                elif char == "\\":
                    escape = True
                elif char == '"':
                    in_string = False
                continue
            if char == '"':
                in_string = True
            elif char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    candidate = stripped[start:idx + 1]
                    try:
                        value = json.loads(candidate)
                        if isinstance(value, dict):
                            return value
                    except json.JSONDecodeError:
                        break

    raise ValueError("Could not find a valid JSON object in agent output")


def normalize_annotation_result(result: Dict[str, Any], cluster: MarkerCluster) -> Dict[str, Any]:
    """Normalize local agent JSON into the CASSIA summary shape."""
    if not result.get("main_cell_type"):
        raise ValueError("Agent JSON is missing required field 'main_cell_type'")

    normalized = dict(result)
    for key in ("sub_cell_types", "possible_mixed_cell_types"):
        value = normalized.get(key)
        if value is None:
            normalized[key] = []
        elif isinstance(value, str):
            normalized[key] = [item.strip() for item in value.split(",") if item.strip()]
        elif not isinstance(value, list):
            normalized[key] = [str(value)]

    normalized["num_markers"] = len(cluster.markers)
    normalized["marker_list"] = cluster.markers
    return normalized


def write_json(path: Path, payload: Any) -> None:
    """Write formatted JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def read_json(path: Path, default: Any) -> Any:
    """Read JSON, returning a default when the file does not exist."""
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def write_summary_csv(run_dir: Path, results: Dict[str, Dict[str, Any]]) -> Path:
    """Write a CASSIA-like summary CSV for local agent CLI output."""
    path = run_dir / "summary.csv"
    headers = [
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
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for cluster_id in sorted(results):
            row = results[cluster_id]
            annotation = row["analysis_result"]
            writer.writerow([
                cluster_id,
                annotation.get("main_cell_type", ""),
                ", ".join(annotation.get("sub_cell_types") or []),
                ", ".join(annotation.get("possible_mixed_cell_types") or []),
                annotation.get("num_markers", ""),
                ", ".join(annotation.get("marker_list") or []),
                row.get("backend", ""),
                annotation.get("confidence", ""),
                annotation.get("evidence", ""),
            ])
    return path


def generate_markdown_report(run_dir: Path) -> Path:
    """Generate a compact Markdown report for a CASSIA CLI run folder."""
    results = read_json(run_dir / "results.json", {})
    errors = read_json(run_dir / "errors.json", [])
    manifest = read_json(run_dir / "run_manifest.json", {})
    summary_rows: List[Dict[str, str]] = []
    if not results:
        summary_candidates = [run_dir / "summary.csv"] + sorted(run_dir.glob("*_summary.csv"))
        for candidate in summary_candidates:
            if candidate.exists():
                with candidate.open(newline="", encoding="utf-8") as handle:
                    summary_rows = list(csv.DictReader(handle))
                break

    lines = [
        "# CASSIA CLI Report",
        "",
        f"- Run directory: `{run_dir}`",
        f"- Backend: `{manifest.get('backend', 'unknown')}`",
        f"- Status: `{manifest.get('status', 'unknown')}`",
        f"- Completed clusters: {len(results) or len(summary_rows)}",
        f"- Failed clusters: {len(errors)}",
        "",
        "| Cluster ID | General Cell Type | Detailed Cell Type | Confidence |",
        "| --- | --- | --- | --- |",
    ]
    for cluster_id in sorted(results):
        annotation = results[cluster_id]["analysis_result"]
        detail = ", ".join(annotation.get("sub_cell_types") or [])
        lines.append(
            f"| {cluster_id} | {annotation.get('main_cell_type', '')} | "
            f"{detail} | {annotation.get('confidence', '')} |"
        )
    for row in summary_rows:
        lines.append(
            f"| {row.get('Cluster ID', '')} | "
            f"{row.get('Predicted General Cell Type', '')} | "
            f"{row.get('Predicted Detailed Cell Type', '')} | "
            f"{row.get('Confidence', '')} |"
        )

    if errors:
        lines.extend(["", "## Errors", ""])
        for error in errors:
            lines.append(f"- `{error.get('cluster_id', 'unknown')}`: {error.get('error', '')}")

    report_path = run_dir / "report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def _manifest_from_args(args: Any, run_dir: Path, mode: str) -> Dict[str, Any]:
    parameters = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in vars(args).items()
        if key not in {"func"}
    }
    return {
        "cassia_version": __version__,
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "status": "running",
        "mode": mode,
        "backend": args.backend,
        "input": str(Path(args.input).resolve()),
        "run_dir": str(run_dir.resolve()),
        "parameters": parameters,
    }


def run_api_annotation(args: Any) -> int:
    """Run annotation through the existing CASSIA API-backed batch function."""
    from CASSIA.engine.tools_function import runCASSIA_batch

    run_dir = Path(args.out) if args.out else default_run_dir()
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = _manifest_from_args(args, run_dir, mode="api")
    write_json(run_dir / "run_manifest.json", manifest)

    output_base = run_dir / "cassia_results.json"
    try:
        runCASSIA_batch(
            marker=args.input,
            output_name=str(output_base),
            n_genes=args.n_genes,
            model=args.model,
            temperature=args.temperature,
            tissue=args.tissue,
            species=args.species,
            additional_info=args.additional_info,
            celltype_column=args.celltype_column,
            gene_column_name=args.gene_column,
            max_workers=args.max_workers,
            provider=args.backend,
            max_retries=args.max_retries,
            ranking_method=args.ranking_method,
            ascending=args.ascending,
            validator_involvement=args.validator_involvement,
            reasoning=args.reasoning,
            use_reference=args.use_reference,
            reference_model=args.reference_model,
            reference_cell_type_hint=args.reference_cell_type_hint,
            validate_api_key_before_start=not args.skip_api_key_validation,
            verbose=not args.quiet,
        )
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["updated_at"] = utc_now()
        manifest["error"] = str(exc)
        write_json(run_dir / "run_manifest.json", manifest)
        raise

    manifest["status"] = "completed"
    manifest["updated_at"] = utc_now()
    manifest["outputs"] = {
        "summary_csv": str(run_dir / "cassia_results_summary.csv"),
        "conversations_json": str(run_dir / "cassia_results_conversations.json"),
        "html_report": str(run_dir / "cassia_results_report.html"),
    }
    write_json(run_dir / "run_manifest.json", manifest)
    return 0


def run_agent_annotation(args: Any) -> int:
    """Run annotation through a local agent CLI backend."""
    run_dir = Path(args.out) if args.out else default_run_dir()
    ensure_run_dirs(run_dir)

    manifest_path = run_dir / "run_manifest.json"
    if args.resume and manifest_path.exists():
        manifest = read_json(manifest_path, {})
        manifest["updated_at"] = utc_now()
        manifest["status"] = "running"
    else:
        manifest = _manifest_from_args(args, run_dir, mode="agent-cli")
    write_json(manifest_path, manifest)

    clusters = load_marker_clusters(
        input_path=Path(args.input),
        n_genes=args.n_genes,
        celltype_column=args.celltype_column,
        gene_column=args.gene_column,
        ranking_method=args.ranking_method,
        ascending=args.ascending,
        limit=args.limit,
    )

    results: Dict[str, Dict[str, Any]] = read_json(run_dir / "results.json", {})
    errors: List[Dict[str, Any]] = read_json(run_dir / "errors.json", [])
    backend = AgentCLIBackend(
        args.backend,
        command_template=args.command_template,
        timeout_seconds=args.timeout,
    )

    completed = set(results)
    for cluster in clusters:
        if args.resume and cluster.cluster_id in completed:
            continue

        prompt = build_annotation_prompt(
            cluster=cluster,
            tissue=args.tissue,
            species=args.species,
            additional_info=args.additional_info,
        )
        prompt_path = run_dir / "prompts" / f"{slugify(cluster.cluster_id)}.md"
        prompt_path.write_text(prompt, encoding="utf-8")

        if args.dry_run:
            continue

        context = {
            "input": str(Path(args.input).resolve()),
            "out": str(run_dir.resolve()),
            "cluster": cluster.cluster_id,
            "agent_output_file": str(run_dir / "raw" / f"{slugify(cluster.cluster_id)}.txt"),
        }
        try:
            raw = backend.run(prompt, prompt_path, run_dir, context)
            raw_path = run_dir / "raw" / f"{slugify(cluster.cluster_id)}.txt"
            raw_path.write_text(raw, encoding="utf-8")
            parsed = extract_json_object(raw)
            annotation = normalize_annotation_result(parsed, cluster)
            results[cluster.cluster_id] = {
                "analysis_result": annotation,
                "backend": args.backend,
                "prompt_file": str(prompt_path),
                "raw_output_file": str(raw_path),
                "completed_at": utc_now(),
            }
            write_json(run_dir / "results.json", results)
            write_summary_csv(run_dir, results)
        except Exception as exc:
            error_record = {
                "cluster_id": cluster.cluster_id,
                "error": str(exc),
                "prompt_file": str(prompt_path),
                "failed_at": utc_now(),
            }
            errors.append(error_record)
            write_json(run_dir / "errors.json", errors)
            (run_dir / "errors" / f"{slugify(cluster.cluster_id)}.txt").write_text(
                str(exc),
                encoding="utf-8",
            )
            if not args.keep_going:
                break

    if args.dry_run:
        manifest["status"] = "dry-run"
    elif errors:
        manifest["status"] = "completed_with_errors" if results else "failed"
    else:
        manifest["status"] = "completed"
    manifest["updated_at"] = utc_now()
    manifest["cluster_count"] = len(clusters)
    manifest["completed_count"] = len(results)
    manifest["error_count"] = len(errors)
    write_json(manifest_path, manifest)

    if results:
        write_summary_csv(run_dir, results)
        generate_markdown_report(run_dir)
    return 1 if errors and not args.dry_run else 0


def resume_run(run_dir: Path) -> int:
    """Resume an agent CLI run using parameters from its manifest."""
    from argparse import Namespace

    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"No run_manifest.json found in {run_dir}")

    manifest = read_json(manifest_path, {})
    if manifest.get("mode") != "agent-cli":
        raise ValueError("Only agent-cli runs can be resumed by the CLI MVP")

    parameters = dict(manifest.get("parameters", {}))
    parameters["resume"] = True
    parameters["out"] = str(run_dir)
    parameters.setdefault("dry_run", False)
    parameters.setdefault("keep_going", True)
    return run_agent_annotation(Namespace(**parameters))


def dispatch_annotation(args: Any) -> int:
    """Dispatch annotate to either API or agent CLI execution."""
    if is_api_backend(args.backend):
        return run_api_annotation(args)
    return run_agent_annotation(args)
