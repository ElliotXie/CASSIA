"""Run-folder orchestration for the CASSIA CLI."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd

from CASSIA import __version__
from CASSIA.core.marker_utils import get_top_markers, split_markers
from CASSIA.engine.main_function_code import (
    construct_prompt,
    final_annotation_system_v1,
    final_annotation_system_v2,
)

from .agent_validation import run_validated_annotation
from .backends import AgentCLIBackend, is_api_backend
from .result_schema import ANNOTATION_SCHEMA_VERSION, normalize_annotation_payload


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
    for child in ("prompts", "raw", "validation", "errors"):
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

    is_preformatted = _looks_preformatted(df, gene_column)
    if is_preformatted:
        prepared = df.copy()
    else:
        prepared = get_top_markers(
            df,
            n_genes=n_genes,
            ranking_method=ranking_method,
            ascending=ascending,
        )

    if is_preformatted:
        cluster_col = celltype_column or prepared.columns[0]
        marker_col = gene_column or prepared.columns[1]
    else:
        cluster_col = "cluster" if "cluster" in prepared.columns else prepared.columns[0]
        marker_col = "markers" if "markers" in prepared.columns else prepared.columns[1]

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


def build_annotation_prompt_v1(
    cluster: MarkerCluster,
    tissue: str,
    species: str,
    additional_info: Optional[str] = None,
) -> str:
    """Build the benchmarked one-shot prompt from the original CASSIA prompt.

    The annotation system prompt and user task come from the core engine.  The
    only addition is a final JSON instruction so an agent CLI can return a
    machine-readable result without a separate LLM formatting call.
    """
    blind = not tissue or tissue.strip().lower() in {"none", "tissue blind"}
    system_source = final_annotation_system_v2 if blind else final_annotation_system_v1
    # Historical benchmark prompts removed source-code-only trailing spaces.
    # Normalize them so v1 is byte-for-byte stable across package builds.
    system = "\n".join(line.rstrip() for line in system_source.splitlines()).strip()
    user_data: Dict[str, Any] = {
        "species": species,
        "tissue_type": tissue,
        "marker_list": list(cluster.markers),
    }
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info
    user_prompt = construct_prompt(user_data)
    possible_tissues = (
        ', "possible_tissues": ["<ranked tissue 1>", "<ranked tissue 2>"]'
        if blind
        else ""
    )
    confidence_instruction = (
        "confidence is an INTEGER 0-10 = your calibrated certainty in the rank-1 subtype "
        "(0 = pure guess, 10 = certain). Be honest and well-spread: reserve 9-10 for "
        "unambiguous canonical calls, use 4-6 when markers fit several subtypes, and 0-3 "
        "when the markers do not pin a specific state."
    )
    parse_tail = (
        "After completing your step-by-step analysis above, output exactly ONE JSON object as the\n"
        "very last line of your response so it can be parsed automatically (do not add anything after it).\n"
        f"Here {confidence_instruction}\n"
        '{"main_cell_type": "<general cell type>", "sub_cell_types": ["<most likely subtype>", "<second>", "<third>"], '
        f'"possible_mixed_cell_types": []{possible_tissues}, "confidence": <integer 0-10>, '
        '"evidence": "<concise marker-based rationale>"}'
    )
    return f"{system}\n\n{user_prompt}\n\n{parse_tail}\n"


def select_annotation_prompt(
    prompt_version: str,
    cluster: MarkerCluster,
    tissue: str,
    species: str,
    additional_info: Optional[str] = None,
) -> str:
    """Build a versioned agent-CLI annotation prompt."""
    if prompt_version == "v1":
        return build_annotation_prompt_v1(cluster, tissue, species, additional_info)
    if prompt_version == "v2":
        return build_annotation_prompt(cluster, tissue, species, additional_info)
    raise ValueError(f"Unknown annotation prompt version: {prompt_version}")


def extract_json_object(text: str) -> Dict[str, Any]:
    """Extract a valid JSON object from agent output.

    Returns the *largest* valid JSON object found in the text. This handles
    common agent-CLI patterns where the model writes a short malformed first
    attempt, narrates a correction, and then writes a complete final JSON —
    in which case the first object is incomplete and the second is the real
    answer. Picking the largest object also handles the simpler case where
    only one object is present.

    Lookup order:
      1. The whole stripped text as a single JSON object.
      2. Code-fenced ``` ```json ``` blocks.
      3. Every balanced ``{ ... }`` span in the text.

    Whichever path produces the most keys (and longest serialization as a
    tie-breaker) wins.
    """
    stripped = text.strip()
    if not stripped:
        raise ValueError("Agent returned an empty response")

    candidates: List[Dict[str, Any]] = []

    try:
        value = json.loads(stripped, strict=False)
        if isinstance(value, dict):
            candidates.append(value)
    except json.JSONDecodeError:
        pass

    fence_matches = re.findall(r"```(?:json)?\s*(\{.*?\})\s*```", stripped, flags=re.DOTALL)
    for candidate in fence_matches:
        try:
            value = json.loads(candidate, strict=False)
            if isinstance(value, dict):
                candidates.append(value)
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
                        value = json.loads(candidate, strict=False)
                        if isinstance(value, dict):
                            candidates.append(value)
                    except json.JSONDecodeError:
                        pass
                    break

    if not candidates:
        raise ValueError("Could not find a valid JSON object in agent output")

    def _score(obj: Dict[str, Any]) -> Tuple[int, int]:
        # Prefer the LONGEST serialized form first. A complete outer wrapper
        # ({"judgments": [...lots of items...]}) is always longer than any
        # individual inner item, and a corrected JSON with N items is always
        # longer than a broken first attempt with fewer items. Use top-level
        # key count as the tie-breaker so two equal-length objects pick the
        # richer schema.
        try:
            serialized_len = len(json.dumps(obj, ensure_ascii=False))
        except (TypeError, ValueError):
            serialized_len = 0
        return (serialized_len, len(obj))

    candidates.sort(key=_score, reverse=True)
    return candidates[0]


def normalize_annotation_result(result: Dict[str, Any], cluster: MarkerCluster) -> Dict[str, Any]:
    """Normalize local agent JSON into the CASSIA summary shape."""
    return normalize_annotation_payload(
        result,
        cluster_id=cluster.cluster_id,
        markers=cluster.markers,
        annotation_mode="standard",
    )


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


def generate_html_report(run_dir: Path) -> Path:
    """Render agent-CLI results with CASSIA's existing batch HTML formatter.

    This is a deterministic presentation step: it reads saved JSON/raw history
    and never calls an LLM, so report generation cannot alter annotations.
    """
    from CASSIA.reports.generate_batch_report import generate_batch_html_report_from_data

    results = read_json(run_dir / "results.json", {})
    if not results:
        raise ValueError("No agent-CLI results are available for HTML reporting")
    manifest = read_json(run_dir / "run_manifest.json", {})
    parameters = manifest.get("parameters", {})
    rows: List[Dict[str, Any]] = []

    for cluster_id in sorted(results):
        result = results[cluster_id]
        annotation = result["analysis_result"]
        annotations: List[str] = []
        validations: List[str] = []
        history_file = result.get("validation_history_file")
        if history_file and Path(history_file).exists():
            history = read_json(Path(history_file), {}).get("history", [])
            annotations = [
                str(item.get("response", ""))
                for item in history
                if item.get("stage") == "annotation"
            ]
            validations = [
                str(item.get("response", ""))
                for item in history
                if item.get("stage") == "validation"
            ]
        else:
            raw_file = result.get("raw_output_file")
            if raw_file and Path(raw_file).exists():
                annotations = [Path(raw_file).read_text(encoding="utf-8")]

        rows.append({
            "Cluster ID": cluster_id,
            "Predicted General Cell Type": annotation.get("main_cell_type", ""),
            "Predicted Detailed Cell Type": ", ".join(annotation.get("sub_cell_types") or []),
            "Possible Mixed Cell Types": ", ".join(
                annotation.get("possible_mixed_cell_types") or []
            ),
            "Marker Number": annotation.get("num_markers", ""),
            "Marker List": ", ".join(annotation.get("marker_list") or []),
            "Iterations": result.get("validation_attempts") or 1,
            "Model": result.get("model") or parameters.get("model") or "default",
            "Provider": result.get("backend") or manifest.get("backend", ""),
            "Tissue": parameters.get("tissue", ""),
            "Species": parameters.get("species", ""),
            "Additional Info": parameters.get("additional_info") or "None",
            "Conversation History": {
                "annotations": annotations,
                "validations": validations,
                # The structured summary is already available. Serializing it
                # here replaces the old LLM formatting-agent call with code.
                "formatting": json.dumps(annotation, ensure_ascii=False),
                "scoring": "",
            },
        })

    report_path = run_dir / "report.html"
    generate_batch_html_report_from_data(
        rows,
        str(report_path),
        report_title="CASSIA Agent CLI Annotation Report",
    )
    return report_path


def _manifest_from_args(args: Any, run_dir: Path, mode: str) -> Dict[str, Any]:
    parameters = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in vars(args).items()
        if key not in {"func"}
    }
    manifest = {
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
    if mode == "agent-cli":
        manifest["result_schema_version"] = ANNOTATION_SCHEMA_VERSION
    return manifest


def run_api_annotation(args: Any) -> int:
    """Run annotation through the existing CASSIA API-backed batch function."""
    from CASSIA.engine.tools_function import runCASSIA_batch

    run_dir = Path(args.out) if args.out else default_run_dir()
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = _manifest_from_args(args, run_dir, mode="api")
    write_json(run_dir / "run_manifest.json", manifest)

    output_base = run_dir / "cassia_results.json"
    marker_input: Any = args.input
    celltype_column = args.celltype_column
    gene_column = args.gene_column
    if args.limit is not None:
        if args.limit < 1:
            raise ValueError("--limit must be at least 1")
        limited_clusters = load_marker_clusters(
            Path(args.input),
            n_genes=args.n_genes,
            celltype_column=args.celltype_column,
            gene_column=args.gene_column,
            ranking_method=args.ranking_method,
            ascending=args.ascending,
            limit=args.limit,
        )
        marker_input = pd.DataFrame(
            {
                "cluster": [cluster.cluster_id for cluster in limited_clusters],
                "markers": [", ".join(cluster.markers) for cluster in limited_clusters],
            }
        )
        celltype_column = "cluster"
        gene_column = "markers"
    try:
        runCASSIA_batch(
            marker=marker_input,
            output_name=str(output_base),
            n_genes=args.n_genes,
            model=args.model,
            temperature=args.temperature,
            tissue=args.tissue,
            species=args.species,
            additional_info=args.additional_info,
            celltype_column=celltype_column,
            gene_column_name=gene_column,
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
        model=getattr(args, "model", None),
        reasoning_effort=getattr(args, "reasoning_effort", None),
    )
    workflow = getattr(args, "workflow", "one-shot")
    prompt_version = getattr(args, "prompt_version", "v1")
    validation_max_attempts = getattr(args, "validation_max_attempts", 3)

    completed = set(results)
    for cluster in clusters:
        if args.resume and cluster.cluster_id in completed:
            continue

        prompt = select_annotation_prompt(
            prompt_version=prompt_version,
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
            raw_path = run_dir / "raw" / f"{slugify(cluster.cluster_id)}.txt"
            validation_history_path: Optional[Path] = None
            validation_passed_value: Optional[bool] = None
            validation_attempts = 0

            if workflow == "validated":
                cluster_slug = slugify(cluster.cluster_id)

                def call_stage(stage: str, attempt: int, stage_prompt: str) -> str:
                    stage_prompt_path = (
                        prompt_path
                        if stage == "annotation" and attempt == 1
                        else run_dir / "prompts" / f"{cluster_slug}_{stage}_{attempt:03d}.md"
                    )
                    stage_raw_path = run_dir / "raw" / f"{cluster_slug}_{stage}_{attempt:03d}.txt"
                    stage_prompt_path.write_text(stage_prompt, encoding="utf-8")
                    stage_context = dict(context)
                    stage_context["agent_output_file"] = str(stage_raw_path)
                    response = backend.run(stage_prompt, stage_prompt_path, run_dir, stage_context)
                    stage_raw_path.write_text(response, encoding="utf-8")
                    return response

                validated = run_validated_annotation(
                    initial_prompt=prompt,
                    marker_list=cluster.markers,
                    tissue=args.tissue,
                    call_agent=call_stage,
                    additional_info=args.additional_info,
                    involvement=getattr(args, "validator_involvement", "v1"),
                    max_attempts=validation_max_attempts,
                    species=args.species,
                )
                raw = validated.final_response
                validation_passed_value = validated.validation_passed
                validation_attempts = validated.validation_attempts
                validation_history_path = run_dir / "validation" / f"{cluster_slug}.json"
                write_json(validation_history_path, {
                    "cluster_id": cluster.cluster_id,
                    "passed": validated.validation_passed,
                    "attempts": validated.validation_attempts,
                    "history": validated.history,
                })
            else:
                raw = backend.run(prompt, prompt_path, run_dir, context)

            raw_path.write_text(raw, encoding="utf-8")
            parsed = extract_json_object(raw)
            annotation = normalize_annotation_result(parsed, cluster)
            results[cluster.cluster_id] = {
                "analysis_result": annotation,
                "backend": args.backend,
                "model": getattr(args, "model", None),
                "workflow": workflow,
                "prompt_version": prompt_version,
                "prompt_file": str(prompt_path),
                "raw_output_file": str(raw_path),
                "validation_passed": validation_passed_value,
                "validation_attempts": validation_attempts,
                "validation_history_file": (
                    str(validation_history_path) if validation_history_path else None
                ),
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
    manifest["formatter"] = "deterministic-python (no LLM call)"
    if results:
        manifest["outputs"] = {
            "results_json": str(run_dir / "results.json"),
            "summary_csv": str(run_dir / "summary.csv"),
            "markdown_report": str(run_dir / "report.md"),
            "html_report": str(run_dir / "report.html"),
        }
    if workflow == "validated":
        manifest["validation_passed_count"] = sum(
            row.get("validation_passed") is True for row in results.values()
        )
        manifest["validation_failed_count"] = sum(
            row.get("validation_passed") is False for row in results.values()
        )
    write_json(manifest_path, manifest)

    if results:
        write_summary_csv(run_dir, results)
        generate_markdown_report(run_dir)
        generate_html_report(run_dir)
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
