"""Annotation boost utilities for the CASSIA CLI."""

from __future__ import annotations

import html
import json
import os
import re
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from CASSIA.engine.main_function_code import final_annotation_system_v1

from .backends import AgentCLIBackend
from .result_schema import ANNOTATION_SCHEMA_VERSION, normalize_fused_boost_payload
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
DEFAULT_AUTO_CONFIDENCE = ("low", "medium", "unknown")
AUTO_KEYWORDS = (
    "ambiguous",
    "unclear",
    "uncertain",
    "mixed",
    "doublet",
    "contaminat",
    "conflict",
    "inconsistent",
    "low confidence",
    "weak evidence",
    "cannot determine",
    "likely wrong",
)


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


def load_marker_table(path: Any) -> pd.DataFrame:
    """Load and lightly normalize a raw marker table or DataFrame."""
    df = path.copy() if isinstance(path, pd.DataFrame) else pd.read_csv(path)
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


def _record_agent_run_metadata(manifest: Dict[str, Any], backend: AgentCLIBackend) -> None:
    """Accumulate structured CLI usage without changing model-visible content."""
    metadata = backend.last_run_metadata or {}
    manifest["llm_calls"] = int(manifest.get("llm_calls", 0) or 0) + 1
    usage = metadata.get("usage") or {}
    if not usage:
        return
    totals = manifest.setdefault("usage", {})
    for key in (
        "input_tokens",
        "cached_input_tokens",
        "output_tokens",
        "reasoning_output_tokens",
    ):
        totals[key] = int(totals.get(key, 0) or 0) + int(usage.get(key, 0) or 0)
    totals["total_tokens"] = int(totals.get("input_tokens", 0)) + int(
        totals.get("output_tokens", 0)
    )


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


def _safe_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.api.types.is_scalar(value) and pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _as_list(value: Any) -> List[str]:
    if value is None:
        return []
    try:
        if pd.api.types.is_scalar(value) and pd.isna(value):
            return []
    except (TypeError, ValueError):
        pass
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    if isinstance(value, tuple):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("[") and text.endswith("]"):
        try:
            parsed = json.loads(text)
            if isinstance(parsed, list):
                return [str(item).strip() for item in parsed if str(item).strip()]
        except Exception:
            pass
    return [item.strip() for item in re.split(r"[,;|]+", text) if item.strip()]


def _normalize_confidence(value: Any) -> str:
    text = _safe_text(value).lower()
    if "high" in text:
        return "high"
    if "medium" in text or "moderate" in text:
        return "medium"
    if "low" in text:
        return "low"
    return "unknown"


def _annotation_search_text(annotation: Dict[str, Any], cluster: str) -> str:
    fields = [
        cluster,
        annotation.get("main_cell_type"),
        annotation.get("final_cell_type"),
        annotation.get("Predicted General Cell Type"),
        annotation.get("final_sub_cell_type"),
        annotation.get("Predicted Detailed Cell Type"),
        annotation.get("confidence"),
        annotation.get("Confidence"),
        annotation.get("evidence"),
        annotation.get("Evidence"),
        annotation.get("recommended_next_steps"),
    ]
    for key in ("sub_cell_types", "possible_mixed_cell_types", "alternatives", "marker_list", "Marker List"):
        fields.extend(_as_list(annotation.get(key)))
    return " ".join(_safe_text(field) for field in fields if _safe_text(field)).lower()


def _candidate_from_annotation(cluster: str, annotation: Dict[str, Any], source: str) -> Dict[str, Any]:
    mixed = _as_list(annotation.get("possible_mixed_cell_types"))
    if not mixed:
        mixed = _as_list(annotation.get("Possible Mixed Cell Types"))
    return {
        "cluster": cluster,
        "source": source,
        "confidence": _normalize_confidence(annotation.get("confidence", annotation.get("Confidence"))),
        "main_cell_type": _safe_text(annotation.get("main_cell_type", annotation.get("Predicted General Cell Type"))),
        "sub_cell_types": _as_list(annotation.get("sub_cell_types", annotation.get("Predicted Detailed Cell Type"))),
        "mixed_cell_types": mixed,
        "evidence": _safe_text(annotation.get("evidence", annotation.get("Evidence"))),
        "annotation": annotation,
    }


def load_boost_candidates(run_dir: Path) -> List[Dict[str, Any]]:
    """Load boost candidate metadata from a CASSIA CLI run directory."""
    results_path = run_dir / "results.json"
    if results_path.exists():
        results = _read_json(results_path, {})
        candidates = []
        for cluster, details in results.items():
            annotation = details.get("analysis_result", details) if isinstance(details, dict) else {}
            if isinstance(annotation, dict):
                candidates.append(_candidate_from_annotation(str(cluster), annotation, str(results_path)))
        if candidates:
            return candidates

    summary_candidates = [run_dir / "summary.csv"] + sorted(run_dir.glob("*_summary.csv"))
    for summary_path in summary_candidates:
        if not summary_path.exists():
            continue
        df = pd.read_csv(summary_path)
        cluster_col = None
        for candidate in ("Cluster ID", "cluster", "True Cell Type"):
            if candidate in df.columns:
                cluster_col = candidate
                break
        if not cluster_col:
            continue
        candidates = []
        for _, row in df.iterrows():
            annotation = {key: _safe_text(value) for key, value in row.to_dict().items()}
            candidates.append(_candidate_from_annotation(str(row[cluster_col]), annotation, str(summary_path)))
        if candidates:
            return candidates

    raise ValueError(f"Could not find results.json or summary.csv in {run_dir}")


def parse_target_terms(values: Optional[Sequence[str]]) -> List[str]:
    """Parse comma-separated target lineage terms without splitting multi-word labels."""
    terms: List[str] = []
    seen = set()
    for value in values or []:
        for term in str(value).split(","):
            cleaned = term.strip()
            key = cleaned.lower()
            if cleaned and key not in seen:
                terms.append(cleaned)
                seen.add(key)
    return terms


def score_boost_candidate(
    candidate: Dict[str, Any],
    confidence_levels: Sequence[str] = DEFAULT_AUTO_CONFIDENCE,
    target_terms: Optional[Sequence[str]] = None,
    select_all: bool = False,
    only_low_confidence: bool = False,
) -> Optional[Dict[str, Any]]:
    """Score one annotation for boost auto-selection."""
    confidence = candidate.get("confidence") or "unknown"
    text = _annotation_search_text(candidate.get("annotation", {}), candidate.get("cluster", ""))
    target_terms = [term.lower() for term in (target_terms or []) if term.strip()]
    reasons: List[str] = []
    score = 0

    if target_terms and not any(term in text for term in target_terms):
        return None

    if only_low_confidence and confidence != "low":
        return None

    if confidence in confidence_levels:
        reasons.append(f"confidence={confidence}")
        score += {"low": 100, "medium": 60, "unknown": 35, "high": 10}.get(confidence, 20)

    if candidate.get("mixed_cell_types"):
        reasons.append("possible mixed cell types")
        score += 50

    matched_keywords = [keyword for keyword in AUTO_KEYWORDS if keyword in text]
    if matched_keywords:
        reasons.append(f"keyword match: {', '.join(matched_keywords[:3])}")
        score += min(45, 15 * len(matched_keywords))

    if target_terms:
        reasons.append(f"target lineage: {', '.join(target_terms)}")
        score += 30

    if select_all and not reasons:
        reasons.append("selected by --all")
        score += 1

    if not reasons:
        return None

    selected = dict(candidate)
    selected["score"] = score
    selected["reasons"] = reasons
    return selected


def select_boost_candidates(
    run_dir: Path,
    confidence_levels: Sequence[str] = DEFAULT_AUTO_CONFIDENCE,
    target_terms: Optional[Sequence[str]] = None,
    max_clusters: Optional[int] = None,
    select_all: bool = False,
    only_low_confidence: bool = False,
) -> List[Dict[str, Any]]:
    """Return scored boost candidates sorted by priority."""
    candidates = []
    for candidate in load_boost_candidates(run_dir):
        scored = score_boost_candidate(
            candidate,
            confidence_levels=confidence_levels,
            target_terms=target_terms,
            select_all=select_all,
            only_low_confidence=only_low_confidence,
        )
        if scored:
            candidates.append(scored)

    candidates.sort(key=lambda item: (-int(item.get("score", 0)), str(item.get("cluster", ""))))
    if max_clusters is not None:
        candidates = candidates[:max_clusters]
    return candidates


def format_annotation_context(context: Dict[str, Any]) -> str:
    """Format original annotation context for the boost prompt."""
    annotation = context.get("annotation", {})
    if isinstance(annotation, dict):
        return json.dumps(annotation, indent=2, ensure_ascii=False)
    return str(annotation)


def extract_check_genes(text: str, max_genes: Optional[int] = None) -> List[str]:
    """Extract requested genes from <check_genes> tags."""
    # Disallow nested angle brackets inside a request. This prevents an explanatory
    # mention such as ``the previous `<check_genes>` panel`` from swallowing all
    # prose up to a later, real closing tag.
    blocks = re.findall(
        r"<check_genes>\s*([^<>]*?)\s*</check_genes>",
        text,
        flags=re.DOTALL | re.IGNORECASE,
    )
    genes: List[str] = []
    for block in blocks:
        genes.extend(
            gene
            for gene in parse_gene_args([block])
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9._-]*", gene)
        )
    if max_genes is not None and max_genes > 0:
        return genes[:max_genes]
    return genes


def extract_candidate_set(text: str) -> List[Dict[str, Any]]:
    """Extract the auditable hypothesis slate emitted by Candidate Boost."""
    match = re.search(
        r"<candidate_set>\s*([\s\S]*?)\s*</candidate_set>",
        text,
        flags=re.IGNORECASE,
    )
    if not match:
        return []
    try:
        payload = json.loads(match.group(1))
    except json.JSONDecodeError:
        return []
    candidates = payload.get("candidates") if isinstance(payload, dict) else payload
    if not isinstance(candidates, list):
        return []
    normalized: List[Dict[str, Any]] = []
    for index, candidate in enumerate(candidates, start=1):
        if isinstance(candidate, str):
            candidate = {"label": candidate}
        if not isinstance(candidate, dict):
            continue
        label = _safe_text(candidate.get("label") or candidate.get("cell_type"))
        if not label:
            continue
        normalized.append({
            "rank": index,
            "label": label,
            "broad_lineage": _safe_text(candidate.get("broad_lineage")),
            "why_plausible": _safe_text(candidate.get("why_plausible")),
            "positive_markers": _as_list(candidate.get("positive_markers")),
            "exclusion_markers": _as_list(candidate.get("exclusion_markers")),
        })
    return normalized


def build_candidate_boost_prompt(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    candidate_count: int = 5,
    additional_task: Optional[str] = None,
    max_genes_per_round: Optional[int] = None,
) -> str:
    """Build the minimal-candidate-first active-evidence experiment prompt."""
    if candidate_count not in {3, 5}:
        raise ValueError("Candidate Boost requires candidate_count to be 3 or 5")
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    gene_rule = (
        f"Request no more than {max_genes_per_round} genes per round."
        if max_genes_per_round is not None and max_genes_per_round > 0
        else "There is no numerical gene cap. Query every marker needed for a decisive comparison, "
        "but keep the panel tied to the candidate slate."
    )
    return f"""You are CASSIA Candidate Boost, an active-evidence single-cell annotator.

Cluster: {cluster}
Dataset context: {major_cluster_info}
Top ranked positive markers from the target-vs-rest differential-expression table:
{", ".join(top_markers)}
{task_text}
This is a candidate-first experiment. Do not start from a long narrative annotation and do
not assume any candidate is correct.

Phase 1 — minimal candidate slate:
1. From only the context and ranked markers above, propose exactly {candidate_count} genuinely
   distinct, conventional cell identities ranked by prior plausibility. Cover the leading
   near-neighbor alternatives; do not fill the slate with cosmetic state variants.
2. For every candidate, name a coherent positive identity program and reciprocal markers that
   would weaken or exclude it. Shared activation, interferon, stress, cell-cycle, ribosomal, and
   mitochondrial programs are not sufficient identity evidence by themselves.
3. Emit the slate in this machine-readable block:
<candidate_set>{{"candidates":[{{"label":"...","broad_lineage":"...","why_plausible":"...","positive_markers":["..."],"exclusion_markers":["..."]}}]}}</candidate_set>
4. Immediately request one discriminating panel using exactly:
<check_genes>GENE1,GENE2,GENE3</check_genes>
5. {gene_rule} Include evidence for the strongest alternatives, not only the current favorite.
6. Evidence boundary: never inspect files/workspaces or run tools. Only use markers in this
   prompt and statistics explicitly returned by CASSIA. Stop after the gene request and wait.

Phase 2 — evidence tournament:
- Compare every candidate against returned enrichment, target prevalence, reference prevalence,
  and coherent multi-gene programs. Absence under dropout is weak unless a program is jointly
  absent/depleted. Expression without enrichment may be shared or ambient.
- You may query another discriminating panel if the top candidates remain unresolved.
- Before finalizing, mark every original candidate supported, weakened, refuted, or unresolved.
  A new candidate may replace the slate only if returned evidence exposes a missed coherent
  identity program; record why.
- The final primary label is not required to equal the initial rank 1.

Final output: return only one JSON object, with no markdown:
{{
  "final_cell_type": "conventional broad identity",
  "final_sub_cell_type": "most likely specific subtype or state",
  "ranked_sub_cell_types": ["most likely", "second", "third"],
  "possible_mixed_cell_types": [],
  "confidence": "low|medium|high",
  "changed_from_original": null,
  "checked_genes": ["GENE1", "GENE2"],
  "supporting_markers": ["GENE1", "GENE2"],
  "refuting_markers": ["GENE3"],
  "alternatives": ["strongest viable alternative"],
  "candidate_audit": [{{"label":"...","status":"supported|weakened|refuted|unresolved","decisive_evidence":"..."}}],
  "evidence": "concise evidence explaining why the winner beat the candidate slate",
  "recommended_next_steps": "optional validation step"
}}

Begin Phase 1 now. Do not output final JSON before CASSIA returns marker-query results.
"""


def build_branch_search_fused_prompt(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    additional_task: Optional[str] = None,
    max_genes_per_round: Optional[int] = None,
) -> str:
    """Build the iterative breadth-then-depth hypothesis-branch experiment."""
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    gene_rule = (
        f"Request no more than {max_genes_per_round} genes per round."
        if max_genes_per_round is not None and max_genes_per_round > 0
        else "There is no numerical gene cap. Use the smallest panel that gives every active "
        "branch a fair positive and reciprocal test; do not query an unfocused marker catalog."
    )
    faithful_annotation_prompt = "\n".join(
        line.rstrip() for line in final_annotation_system_v1.splitlines()
    ).strip()
    return f"""{faithful_annotation_prompt}

CASSIA ACTIVE-EVIDENCE EXTENSION — ITERATIVE BRANCH SEARCH

You are the primary annotator in a fused Annotation Boost session. No prior annotation is
available or trusted. Preserve the original CASSIA functional-marker, cell-type-marker,
general-type, and ranked-subtype reasoning, but organize active evidence as a mutable
breadth-then-depth hypothesis search. The hypotheses are search branches, not answers.

Cluster: {cluster}
Dataset context: {major_cluster_info}
Top ranked positive markers from the target-vs-rest differential-expression table:
{", ".join(top_markers)}
{task_text}
Search protocol:
1. Start with exactly three genuinely distinct conventional identity hypotheses. Prefer
   meaningful lineage or sibling-subtype competitors; do not fill the slate with cosmetic
   activation/state variants.
2. For each branch, specify:
   - the coherent positive identity program that would support it,
   - reciprocal or sibling markers that would weaken it,
   - which visible ranked markers motivated the branch.
3. BREADTH ROUND: request one combined panel that fairly tests all three branches. Every branch
   must contribute positive markers and at least one useful reciprocal discriminator. Emit:
<branch_ledger>{{"branches":[{{"label":"...","visible_basis":["..."],"positive_markers":["..."],"reciprocal_markers":["..."]}}]}}</branch_ledger>
<check_genes>GENE1,GENE2,GENE3</check_genes>
4. {gene_rule} Use official gene symbols. Stop after the request and wait for CASSIA results.
5. After every returned panel, update every branch as supported, weakened, refuted, or unresolved.
   Do not equate one absent marker with refutation under dropout; judge coherent programs using
   enrichment, target prevalence, reference prevalence, and reciprocal evidence.
6. OPEN BRANCH RULE: inspect the ranked and returned evidence for a coherent program unexplained
   by all current branches. A new branch may enter only from such positive unexplained evidence;
   record which old branch it replaces and why. This prevents the initial three from becoming a
   closed world without encouraging generic marker fishing.
7. DEPTH ROUND: after broad comparison, select the leading branch and its strongest surviving
   rival. Request a second focused panel that tests stable subtype/lineage discriminators between
   them and attempts to falsify the leader. Separate identity from activation, interferon, stress,
   cell cycle, maturation, and anatomical state.
8. Complete both the breadth and depth evidence rounds before finalizing. Additional rounds are
   allowed only when a specific unresolved branch or newly exposed coherent program justifies them.
9. Final head-to-head gate: choose the conventional broad identity and rank-1 subtype that best
   explain all evidence. Unsupported precision loses to a broader supported label. Mixed/doublet
   calls require two coherent incompatible identity programs.
10. Evidence boundary: never inspect files/workspaces or run tools yourself. Only use ranked
   markers above and statistics explicitly returned by CASSIA.

Final output: return only one JSON object, with no markdown:
{{
  "final_cell_type": "conventional broad identity",
  "final_sub_cell_type": "most likely specific subtype or state",
  "ranked_sub_cell_types": ["most likely", "second", "third"],
  "possible_mixed_cell_types": [],
  "confidence": "low|medium|high",
  "changed_from_original": null,
  "checked_genes": ["GENE1", "GENE2"],
  "supporting_markers": ["GENE1", "GENE2"],
  "refuting_markers": ["GENE3"],
  "alternatives": ["strongest surviving alternative"],
  "branch_audit": [{{"label":"...","status":"supported|weakened|refuted|unresolved","decisive_evidence":"..."}}],
  "evidence": "concise breadth-then-depth evidence explaining the winner",
  "recommended_next_steps": "optional validation step"
}}

Begin with the three-branch ledger and breadth-round gene request. Do not output final JSON yet.
"""


def _experimental_fused_output_schema() -> str:
    """Return the shared parseable schema for answer-agnostic fused experiments."""
    return """{
  "final_cell_type": "conventional broad identity",
  "final_sub_cell_type": "most likely specific subtype or state",
  "ranked_sub_cell_types": ["most likely", "second", "third"],
  "possible_mixed_cell_types": [],
  "confidence": "low|medium|high",
  "changed_from_original": null,
  "checked_genes": ["GENE1", "GENE2"],
  "supporting_markers": ["GENE1", "GENE2"],
  "refuting_markers": ["GENE3"],
  "alternatives": ["strongest viable alternative"],
  "evidence": "concise evidence from ranked and queried marker statistics",
  "recommended_next_steps": "optional validation step"
}"""


def build_open_world_falsification_prompt(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    additional_task: Optional[str] = None,
    max_genes_per_round: Optional[int] = None,
) -> str:
    """Build an anti-anchoring, open-world counterfactual fused experiment."""
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    gene_rule = (
        f"Request no more than {max_genes_per_round} genes per round."
        if max_genes_per_round is not None and max_genes_per_round > 0
        else "There is no numerical gene cap; use a compact hypothesis-driven panel."
    )
    return f"""You are CASSIA Open-World Falsification, an active-evidence single-cell annotator.

Cluster: {cluster}
Dataset context: {major_cluster_info}
Top ranked positive markers from the target-vs-rest differential-expression table:
{", ".join(top_markers)}
{task_text}
Objective: identify the conventional broad lineage and rank-1 subtype while preventing the
initial shortlist from anchoring the final answer. No prior annotation exists.

First response:
1. Form up to three provisional, genuinely distinct identity hypotheses. They are search handles,
   not a closed candidate list and not a ranking commitment.
2. For each, identify a coherent positive program and reciprocal evidence that would refute it.
3. Add an explicit OPEN-WORLD probe: genes capable of revealing a coherent identity outside all
   provisional hypotheses. This probe must be biologically motivated by unexplained ranked markers,
   tissue context, or a plausible competing lineage; do not query a generic catalog.
4. Request one combined discriminating panel using exactly:
<check_genes>GENE1,GENE2,GENE3</check_genes>
5. {gene_rule} Include positive and reciprocal markers. Never inspect files or use tools yourself.
   Stop after the gene request and wait for statistics returned by CASSIA.

After each result:
- Treat enrichment, target prevalence, reference prevalence, and multi-gene coherence as evidence.
  A missing single marker under dropout is weak; expression without enrichment may be shared/ambient.
- Run a mandatory null-slate reconstruction: temporarily ignore the original hypothesis ranking and
  ask which identity best explains all ranked plus queried evidence. A final identity outside the
  provisional set has exactly the same burden of proof as one inside it.
- Explicitly seek the strongest counterexample to the current leader. Query another focused panel
  when the winner is supported only by shared state markers or a viable alternative remains.
- Separate stable identity from activation, interferon, stress, cell cycle, location, and maturation.
  Prefer a conventional broader subtype over unsupported precision.

When evidence is sufficient, return only this JSON shape with no markdown:
{_experimental_fused_output_schema()}

Begin with the provisional hypotheses, the open-world probe rationale, and the first gene request.
Do not output final JSON before CASSIA returns marker-query results.
"""


def build_program_first_fused_prompt(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    additional_task: Optional[str] = None,
    max_genes_per_round: Optional[int] = None,
) -> str:
    """Build a program-first experiment that delays cell-type naming until evidence returns."""
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    gene_rule = (
        f"Request no more than {max_genes_per_round} genes per round."
        if max_genes_per_round is not None and max_genes_per_round > 0
        else "There is no numerical gene cap; query every gene needed to complete or reject the programs."
    )
    return f"""You are CASSIA Program-First Inference, an active-evidence single-cell annotator.

Cluster: {cluster}
Dataset context: {major_cluster_info}
Top ranked positive markers from the target-vs-rest differential-expression table:
{", ".join(top_markers)}
{task_text}
The first round is phenotype-blind: do not emit or rank any cell-type, lineage, or subtype names
before CASSIA returns the first marker-query statistics. This prevents label-first confirmation.

First response:
1. Partition the observed markers into molecular programs without naming a phenotype: stable
   identity/structure, effector/function, signaling/state, proliferation/stress, housekeeping,
   possible ambient signal, and any mutually incompatible program.
2. Identify which stable program is incomplete or ambiguous. Select completion genes expected to
   co-enrich if it is real, plus reciprocal exclusion genes for the most plausible competing program.
   Famous single markers are insufficient; test coherent modules.
3. Emit a short machine-readable plan without phenotype labels:
<program_plan>{{"observed_programs":["..."],"unexplained_markers":["..."],"decision_tests":["..."]}}</program_plan>
4. Request the completion/exclusion panel using exactly:
<check_genes>GENE1,GENE2,GENE3</check_genes>
5. {gene_rule} Never inspect files or use tools yourself. Stop and wait for CASSIA.

After the first result:
- Only now map coherent stable programs to broad lineage hypotheses, then subtype siblings. Keep
  transient state separate from identity.
- Weight multi-gene enrichment and target prevalence; down-weight isolated expression, ambient
  genes, shared activation, and absence of a single dropout-prone marker.
- Attempt to falsify the leading mapping with reciprocal markers. If subtype siblings remain
  unresolved, request one second focused panel rather than manufacturing a precise label.
- Mixed/doublet calls require two coherent incompatible identity programs.

When evidence is sufficient, return only this JSON shape with no markdown:
{_experimental_fused_output_schema()}

Begin with the phenotype-blind program plan and first gene request. Do not name a cell type and do
not output final JSON before CASSIA returns marker-query results.
"""


def build_fused_boost_prompt(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    strategy: str = "breadth",
    additional_task: Optional[str] = None,
    max_genes_per_round: Optional[int] = None,
    prompt_variant: str = "v2-compact",
    candidate_count: int = 5,
) -> str:
    """Fuse the faithful annotation prompt with active full-marker querying."""
    if prompt_variant == "v2-compact":
        prompt = build_fused_boost_prompt(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            strategy=strategy,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
            prompt_variant="v2",
            candidate_count=candidate_count,
        )
        extension = "CASSIA ACTIVE-EVIDENCE EXTENSION\n"
        _, separator, suffix = prompt.partition(extension)
        if not separator:
            raise ValueError("Could not isolate the Fused v2 active-evidence extension")
        return (
            (extension + suffix)
            .replace(
                "Perform the original CASSIA annotation reasoning above, while actively querying",
                "Perform a preliminary annotation from the supplied markers, while actively querying",
                1,
            )
            .replace(
                "Perform the original CASSIA functional-marker, cell-type-marker, general-type, and top-three-subtype analysis.",
                "Analyze functional markers, cell-type markers, general type, and top-three subtypes.",
                1,
            )
        )
    if prompt_variant in {"gsea_tool", "ucell_tool"}:
        base = build_fused_boost_prompt(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            strategy=strategy,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
            prompt_variant="v2",
            candidate_count=candidate_count,
        )
        if prompt_variant == "gsea_tool":
            tool_name = "weighted preranked GSEA"
            request_tag = "gsea_request"
            evidence_description = (
                "CASSIA will score each requested signature against the complete signed "
                "target-vs-rest log-fold-change ranking and return ES, NES, a deterministic "
                "gene-set-permutation p value, and leading-edge genes."
            )
        else:
            tool_name = "UCell"
            request_tag = "ucell_request"
            evidence_description = (
                "CASSIA will calculate the published per-cell Mann-Whitney rank score "
                "(maxRank=1500) in anonymously sampled target and reference cells and return "
                "score distributions plus probability of superiority."
            )
        return base + f"""

CASSIA SIGNATURE-EVIDENCE TOOL — {tool_name.upper()}

This experiment keeps the Fused v2 annotation task and output schema unchanged, but gives you
one deterministic signature-level evidence tool. For the first evidence request, use this tool
instead of <check_genes>. Propose 2-4 genuinely competing, answer-agnostic signatures with 3-50
official gene symbols each. Include coherent programs for the leading identity and its strongest
sibling/lineage alternative; do not create a signature from only the visible top markers.

Request exactly one JSON block:
<{request_tag}>{{"signatures":[{{"name":"short hypothesis name","genes":["GENE1","GENE2","GENE3"]}}]}}</{request_tag}>

{evidence_description}

After receiving results, use signature evidence as an aid rather than an automatic label:
- Coherent multi-gene separation matters more than one famous marker.
- State programs must not replace stable identity programs.
- A non-enriched signature may be incomplete or dropout-sensitive; compare alternatives.
- You may then request a focused <check_genes> panel or one more signature-tool round if needed.
- Never infer the hidden target label; the tool does not disclose it.

Override the earlier startup sentence for this experiment: begin with the mandatory
<{request_tag}> request, stop, and wait. Do not output final JSON in the first response.
"""
    if prompt_variant == "candidate":
        return build_candidate_boost_prompt(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            candidate_count=candidate_count,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
        )
    if prompt_variant == "branch_search":
        return build_branch_search_fused_prompt(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
        )
    if prompt_variant == "open_world":
        return build_open_world_falsification_prompt(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
        )
    if prompt_variant == "program_first":
        return build_program_first_fused_prompt(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
        )
    if prompt_variant == "v3":
        return build_fused_boost_prompt_v3(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            strategy=strategy,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
        )
    if prompt_variant not in {"v2", "v14"}:
        raise ValueError(f"Unknown fused Boost prompt variant: {prompt_variant}")
    strategy_text = (
        "Use a depth-first strategy: investigate one leading hypothesis at a time, "
        "then go deeper into subtype/state if it is supported."
        if strategy == "depth"
        else "Use a breadth-first strategy: maintain up to three plausible cell-type or state "
        "hypotheses, then choose decisive positive and negative markers that separate them."
    )
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    gene_rule = (
        f"Request no more than {max_genes_per_round} genes per round."
        if max_genes_per_round is not None and max_genes_per_round > 0
        else "There is no numerical gene cap. Request every gene needed for a decisive comparison, "
        "while keeping each panel hypothesis-driven rather than exhaustive."
    )
    faithful_annotation_prompt = "\n".join(
        line.rstrip() for line in final_annotation_system_v1.splitlines()
    ).strip()
    return f"""{faithful_annotation_prompt}

CASSIA ACTIVE-EVIDENCE EXTENSION

You are the primary annotator in a fused Annotation Boost session. No prior annotation is available or trusted.
Perform the original CASSIA annotation reasoning above, while actively querying the full
target-vs-rest marker table whenever the supplied top-ranked markers do not decisively distinguish
the leading hypotheses.

Cluster: {cluster}
Dataset context: {major_cluster_info}
Top ranked markers from the raw differential expression table:
{", ".join(top_markers)}
{task_text}
Active-evidence workflow:
1. Perform the original CASSIA functional-marker, cell-type-marker, general-type, and top-three-subtype analysis.
2. {strategy_text}
3. Consider mixed populations, doublets, transitional states, and ambient RNA only when supported by coherent evidence.
4. Request local marker statistics using exactly:
<check_genes>GENE1,GENE2,GENE3</check_genes>
5. {gene_rule} Use official gene symbols and include both confirming and refuting markers where useful.
6. Evidence boundary: do not inspect the filesystem or workspace, run shell/browser tools, or invent query results. The only valid evidence is the ranked markers in this prompt and statistics explicitly returned by CASSIA after a <check_genes> request. After requesting genes, stop and wait; never simulate the results in the same response.
7. Always complete at least one marker-query round before finalizing. After each result, refine, retain, or pivot; do not repeat genes unnecessarily.
8. When the evidence is sufficient, return only one valid JSON object with this exact schema:
{{
  "final_cell_type": "general cell type",
  "final_sub_cell_type": "most likely specific subtype or state",
  "ranked_sub_cell_types": ["most likely", "second", "third"],
  "possible_mixed_cell_types": [],
  "confidence": "low|medium|high",
  "changed_from_original": null,
  "checked_genes": ["GENE1", "GENE2"],
  "supporting_markers": ["GENE1", "GENE2"],
  "refuting_markers": ["GENE3"],
  "alternatives": ["alternative if confidence is not high"],
  "evidence": "concise evidence based on ranked and queried marker statistics",
  "recommended_next_steps": "optional next validation step"
}}

Start with the original CASSIA preliminary analysis and the first targeted <check_genes> request.
Do not output final JSON before receiving marker-query results.
"""


def build_fused_boost_prompt_v3(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    strategy: str = "breadth",
    additional_task: Optional[str] = None,
    max_genes_per_round: Optional[int] = None,
) -> str:
    """Build the hierarchy-first, evidence-calibrated fused Boost prompt."""
    strategy_text = (
        "Use depth-first search after the broad lineage is secure: test the leading lineage, "
        "then resolve subtype or state within it."
        if strategy == "depth"
        else "Maintain up to three genuinely distinct hypotheses until a targeted panel "
        "separates them; do not create cosmetic variants of the same hypothesis."
    )
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    gene_rule = (
        f"Request no more than {max_genes_per_round} genes per round."
        if max_genes_per_round is not None and max_genes_per_round > 0
        else "There is no numerical gene cap. Keep panels hypothesis-driven and request every "
        "gene needed for the decision, but do not query broad unfocused gene catalogs."
    )
    faithful_annotation_prompt = "\n".join(
        line.rstrip() for line in final_annotation_system_v1.splitlines()
    ).strip()
    return f"""{faithful_annotation_prompt}

CASSIA ACTIVE-EVIDENCE EXTENSION — HIERARCHICAL CALIBRATION

You are the primary annotator in a fused Annotation Boost session. No prior annotation is
available or trusted. Solve this cluster independently from its molecular evidence; do not
infer a desired answer from benchmark conventions or optimize for any named evaluation set.

Cluster: {cluster}
Dataset context: {major_cluster_info}
Top ranked positive markers from the target-vs-rest differential-expression table:
{", ".join(top_markers)}
{task_text}
Decision protocol:
1. Separate three levels that must not be conflated:
   A. broad lineage/general cell type,
   B. stable subtype identity,
   C. transient state, activation, stress, cell cycle, or location.
   Shared effector, interferon, stress, ribosomal, mitochondrial, and cell-cycle programs are
   state evidence unless accompanied by a coherent identity program.
2. {strategy_text}
3. First establish or challenge the broad lineage. Query coherent positive programs for the
   leading hypotheses plus reciprocal exclusion markers. Only then spend evidence on fine
   subtype/state resolution. A single famous marker is not a program.
4. Request local marker statistics using exactly:
<check_genes>GENE1,GENE2,GENE3</check_genes>
5. {gene_rule} Use official gene symbols. Include decisive markers for the strongest alternative,
   not only confirmatory markers for the current favorite.
6. Interpret returned statistics carefully:
   - Strong positive identity evidence combines enrichment, meaningful target prevalence, and
     multiple biologically coherent markers.
   - A marker's absence is weak under dropout. Treat negative evidence as decisive only when
     several expected program members are absent/depleted or the marker is detectably expressed
     in the reference population.
   - Expression without target enrichment may reflect a shared program or ambient RNA.
   - Mixed/doublet calls require two coherent, incompatible identity programs; one stray marker
     is insufficient.
7. Evidence boundary: never inspect files/workspaces or run tools yourself. Only use the ranked
   markers above and statistics explicitly returned by CASSIA. After a gene request, stop and wait.
8. Complete at least one marker-query round. If the strongest alternative remains viable after
   the first result, query a second discriminating panel instead of forcing certainty.
9. Before finalizing, run this counterfactual decision gate internally:
   - What is the strongest alternative?
   - Which returned evidence distinguishes the primary call from it?
   - Is the proposed subtype supported by identity markers, or only by a shared state program?
   If exact subtype evidence is not discriminative, keep the correct broad identity and use a
   conventional broader subtype with medium/low confidence. Do not invent a composite subtype
   from unrelated adjectives merely to sound precise.
10. Return only one JSON object:
{{
  "final_cell_type": "conventional broad identity supported by a coherent program",
  "final_sub_cell_type": "specific identity only when discriminative evidence supports it",
  "ranked_sub_cell_types": ["most likely", "second", "third"],
  "possible_mixed_cell_types": [],
  "confidence": "low|medium|high",
  "changed_from_original": null,
  "checked_genes": ["GENE1", "GENE2"],
  "supporting_markers": ["GENE1", "GENE2"],
  "refuting_markers": ["GENE3"],
  "alternatives": ["strongest viable alternative"],
  "evidence": "concise lineage-then-subtype evidence using ranked and queried statistics",
  "recommended_next_steps": "optional next validation step"
}}

Start with the lineage-level hypotheses and the first discriminating <check_genes> request.
Do not output final JSON before receiving marker-query results.
"""


def build_boost_prompt(
    cluster: str,
    major_cluster_info: str,
    top_markers: Sequence[str],
    annotation_context: Optional[str] = None,
    strategy: str = "breadth",
    additional_task: Optional[str] = None,
    max_genes_per_round: Optional[int] = None,
    mode: str = "review",
    prompt_variant: str = "v2-compact",
    candidate_count: int = 5,
) -> str:
    """Build the initial annotation boost prompt for an agent CLI backend."""
    if mode == "fused":
        return build_fused_boost_prompt(
            cluster=cluster,
            major_cluster_info=major_cluster_info,
            top_markers=top_markers,
            strategy=strategy,
            additional_task=additional_task,
            max_genes_per_round=max_genes_per_round,
            prompt_variant=prompt_variant,
            candidate_count=candidate_count,
        )
    if mode != "review":
        raise ValueError(f"Unknown Annotation Boost mode: {mode}")
    if not annotation_context:
        raise ValueError("Review mode requires an original annotation context")
    strategy_text = (
        "Use a depth-first strategy: investigate one leading hypothesis at a time, "
        "then go deeper into subtype/state if it is supported."
        if strategy == "depth"
        else "Use a breadth-first strategy: consider up to three plausible hypotheses, "
        "then choose decisive positive and negative markers to separate them."
    )
    task_text = f"\nAdditional task: {additional_task}\n" if additional_task else ""
    gene_rule = (
        f"Request no more than {max_genes_per_round} genes per round."
        if max_genes_per_round is not None and max_genes_per_round > 0
        else "There is no numerical gene cap. Request every gene needed for a decisive comparison, "
        "while keeping each panel hypothesis-driven rather than exhaustive."
    )
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
4. {gene_rule} Use official gene symbols only.
5. Evidence boundary: do not inspect the filesystem or workspace, run shell/browser tools, or invent query results. Only use evidence in this prompt and marker statistics explicitly returned by CASSIA. After requesting genes, stop and wait for the results.
6. After CASSIA returns marker statistics, refine or pivot. Do not repeat already checked genes unless necessary.
7. When ready, return only one valid JSON object with this exact schema:
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


def build_boost_followup_prompt(
    transcript: str,
    query_text: Optional[str],
    is_final_round: bool,
    prompt_variant: str = "v2-compact",
) -> str:
    """Build a stateless follow-up prompt containing the transcript so far."""
    final_instruction = (
        "You are at the final round. Return only the final JSON object now."
        if is_final_round
        else "Continue the boost analysis. Either request another <check_genes> list or return the final JSON object."
    )
    query_section = f"\nLatest CASSIA marker query results:\n{query_text}\n" if query_text else ""
    calibration = (
        "- Re-apply the hierarchy: broad lineage first, stable subtype second, transient state third.\n"
        "- Before finalizing, identify the strongest alternative and the returned evidence that distinguishes it.\n"
        "- Prefer a conventional broader subtype over an unsupported precise or composite label.\n"
        if prompt_variant == "v3"
        else (
            "- Interpret GSEA at the program level: prioritize coherent positive NES and leading-edge genes, not the name you gave the signature.\n"
            "- You may request another <gsea_request> JSON block or a focused <check_genes> panel if alternatives remain unresolved.\n"
            "- Do not let a shared pathway or state signature replace stable identity evidence.\n"
            if prompt_variant == "gsea_tool"
            else (
                "- Interpret UCell distributions and target-vs-reference separation; do not treat a score as an automatic cell-type label.\n"
                "- You may request another <ucell_request> JSON block or a focused <check_genes> panel if alternatives remain unresolved.\n"
                "- Prefer coherent stable-identity signatures over shared activation/state programs.\n"
                if prompt_variant == "ucell_tool"
                else (
            "- Reconstruct the answer from all evidence without privileging the initial shortlist.\n"
            "- Give an identity outside the shortlist the same evidentiary burden as one inside it.\n"
            "- Seek the strongest counterexample; do not let shared state markers decide identity.\n"
            if prompt_variant == "open_world"
            else (
                "- The first returned statistics now permit phenotype naming: map coherent stable programs to lineage first, then subtype.\n"
                "- Falsify the leading mapping with reciprocal markers and keep state separate from identity.\n"
                "- If subtype siblings remain unresolved, query a focused second panel or use a conventional broader label.\n"
                if prompt_variant == "program_first"
                else (
            "- Preserve the original candidate slate as an explicit tournament. Compare every candidate "
            "with the returned statistics and do not simply defend the initial rank 1.\n"
            "- Before finalizing, include candidate_audit entries for the original slate and explain the "
            "decisive evidence that selected the winner.\n"
            "- A candidate outside the slate is allowed only when returned evidence reveals a coherent "
            "missed identity program.\n"
            if prompt_variant == "candidate"
            else (
            "- Maintain a mutable branch ledger: update every branch as supported, weakened, refuted, or unresolved.\n"
            "- After the breadth panel, run a distinct depth panel on the leader versus its strongest surviving rival.\n"
            "- Admit a new branch only when a coherent positive program is unexplained by the current branches.\n"
            "- Do not finalize until both breadth and depth marker-query rounds are complete.\n"
            if prompt_variant == "branch_search"
            else (
            "- Final sufficiency gate: preserve the stable identity supported by the direct ranked-marker program; "
            "pivot only when returned queries provide a coherent reciprocal program, not an isolated post-hoc marker.\n"
            "- A precise segment, maturation/state, or anatomical qualifier needs discriminative evidence against "
            "its strongest sibling. Do not introduce an anatomical region outside the supplied dataset context.\n"
            "- If identity evidence is mixed, use a conventional supported identity as the primary label and put "
            "uncertain state or location wording in evidence/alternatives.\n"
            if prompt_variant == "v14"
            else ""
            )
            )
                )
            )
                )
            )
        )
    )
    return f"""Continue this CASSIA annotation boost session.

Transcript so far:
{transcript}
{query_section}
{final_instruction}

Remember:
- Use <check_genes>GENE1,GENE2</check_genes> for more local marker checks.
- Do not inspect files or use tools yourself; only use marker statistics explicitly returned in this prompt.
- If requesting genes, stop after the <check_genes> request and wait for CASSIA results.
{calibration}- Final output must be one valid JSON object only, with no markdown fences or commentary.
"""


def _append_transcript(messages: List[Dict[str, str]]) -> str:
    parts = []
    for msg in messages:
        parts.append(f"## {msg['role'].upper()}\n\n{msg['content']}")
    return ("\n\n" + ("=" * 80) + "\n\n").join(parts)


def normalize_boost_result(result: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize final boost JSON into both boost and canonical CASSIA fields."""
    return normalize_fused_boost_payload(result)


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
    """Write the mode-appropriate deterministic HTML report from boost artifacts."""
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
    if manifest.get("mode") == "fused":
        from CASSIA.reports.fused_boost_report import write_fused_boost_html_report

        return write_fused_boost_html_report(
            output_path=html_path,
            manifest=manifest,
            messages=messages,
            result=result,
            errors=errors,
            query_frames=query_frames,
            top_marker_rows=top_marker_rows,
            quality_assessment=manifest.get("quality_assessment"),
        )

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


def regenerate_boost_reports(boost_dir: Path) -> Dict[str, str]:
    """Rebuild Markdown and HTML reports from saved boost artifacts only."""
    boost_dir = Path(boost_dir)
    manifest_path = boost_dir / "boost_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"No boost_manifest.json found in {boost_dir}")

    manifest = _read_json(manifest_path, {})
    final_result = _read_final_result(boost_dir)
    if final_result:
        final_result = normalize_boost_result(final_result)
    errors = [str(item) for item in manifest.get("errors", [])]

    messages: List[Dict[str, str]] = []
    prompt_files = sorted((boost_dir / "prompts").glob("round_*.md"))
    for prompt_path in prompt_files:
        messages.append({
            "role": "user",
            "content": prompt_path.read_text(encoding="utf-8"),
        })
        raw_path = boost_dir / "raw" / f"{prompt_path.stem}.txt"
        if raw_path.exists():
            messages.append({
                "role": "assistant",
                "content": raw_path.read_text(encoding="utf-8"),
            })
    final_prompt = boost_dir / "prompts" / "final.md"
    if final_prompt.exists():
        messages.append({"role": "user", "content": final_prompt.read_text(encoding="utf-8")})
        final_raw = boost_dir / "raw" / "final.txt"
        if final_raw.exists():
            messages.append({"role": "assistant", "content": final_raw.read_text(encoding="utf-8")})

    query_frames = [
        pd.read_csv(path)
        for path in sorted((boost_dir / "queries").glob("round_*.csv"))
    ]
    top_markers_path = boost_dir / "top_markers.csv"
    top_marker_rows = (
        pd.read_csv(top_markers_path)
        if top_markers_path.exists()
        else pd.DataFrame(columns=["gene"])
    )

    cluster = str(manifest.get("cluster", ""))
    if manifest.get("mode") == "fused":
        annotation_context = {
            "source": None,
            "cluster_id": cluster,
            "annotation": {"mode": "fused_primary_annotation", "prior_annotation": None},
        }
    else:
        try:
            annotation_context = load_annotation_context(
                Path(manifest.get("run_dir", boost_dir)),
                cluster,
            )
        except Exception:
            annotation_context = {
                "source": manifest.get("annotation_source"),
                "cluster_id": cluster,
                "annotation": "Original annotation artifact is unavailable.",
            }

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
        strategy=manifest.get("parameters", {}).get("strategy", "breadth"),
    )
    manifest["markdown_report"] = str(report_path)
    manifest["html_report"] = str(html_path)
    manifest["summary_tags"] = str(boost_dir / "summary_tags.txt")
    manifest["outputs"] = {
        "final_json": str(boost_dir / "final.json") if final_result else None,
        "markdown_report": str(report_path),
        "html_report": str(html_path),
        "summary_tags": str(boost_dir / "summary_tags.txt"),
        "transcript": str(boost_dir / "transcript.md"),
    }
    manifest["updated_at"] = utc_now()
    _write_json(manifest_path, manifest)
    return {"markdown": str(report_path), "html": str(html_path)}


def _auto_result_row(result: Dict[str, Any]) -> Dict[str, Any]:
    final = result.get("final_result") or {}
    candidate = result.get("candidate") or {}
    return {
        "cluster": candidate.get("cluster", ""),
        "status": result.get("status", ""),
        "score": candidate.get("score", ""),
        "original_confidence": candidate.get("confidence", ""),
        "selection_reasons": "; ".join(candidate.get("reasons", [])),
        "original_cell_type": candidate.get("main_cell_type", ""),
        "boost_cell_type": final.get("final_cell_type", ""),
        "boost_sub_cell_type": final.get("final_sub_cell_type", ""),
        "boost_confidence": final.get("confidence", ""),
        "changed_from_original": final.get("changed_from_original", ""),
        "boost_dir": result.get("boost_dir", ""),
        "html_report": result.get("html_report", ""),
        "error": result.get("error", ""),
    }


def write_auto_report(auto_dir: Path, manifest: Dict[str, Any], run_results: Sequence[Dict[str, Any]]) -> Dict[str, str]:
    """Write aggregate CSV, Markdown, and HTML reports for boost auto."""
    auto_dir.mkdir(parents=True, exist_ok=True)
    rows = [_auto_result_row(result) for result in run_results]
    csv_path = auto_dir / "auto_summary.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    md_lines = [
        "# CASSIA Boost Auto Report",
        "",
        f"- Status: `{manifest.get('status')}`",
        f"- Selected clusters: {len(manifest.get('selected_candidates', []))}",
        f"- Completed: {sum(1 for row in rows if row['status'] == 'completed')}",
        f"- Failed: {sum(1 for row in rows if row['status'] == 'failed')}",
        f"- Skipped: {sum(1 for row in rows if row['status'] == 'skipped')}",
        "",
        "| Cluster | Status | Original Confidence | Reasons | Boost Annotation | Boost Confidence |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        boost_label = row["boost_cell_type"]
        if row["boost_sub_cell_type"]:
            boost_label = f"{boost_label} / {row['boost_sub_cell_type']}" if boost_label else row["boost_sub_cell_type"]
        md_lines.append(
            "| {cluster} | {status} | {original_confidence} | {reasons} | {boost_label} | {boost_confidence} |".format(
                cluster=row["cluster"],
                status=row["status"],
                original_confidence=row["original_confidence"],
                reasons=row["selection_reasons"],
                boost_label=boost_label,
                boost_confidence=row["boost_confidence"],
            )
        )
    md_path = auto_dir / "auto_report.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    html_rows = []
    for row in rows:
        link = row["html_report"]
        if link:
            try:
                link = os.path.relpath(Path(link), auto_dir)
            except ValueError:
                pass
        report_link = f'<a href="{html.escape(link)}">summary.html</a>' if link else ""
        html_rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['cluster']))}</td>"
            f"<td>{html.escape(str(row['status']))}</td>"
            f"<td>{html.escape(str(row['score']))}</td>"
            f"<td>{html.escape(str(row['original_confidence']))}</td>"
            f"<td>{html.escape(str(row['selection_reasons']))}</td>"
            f"<td>{html.escape(str(row['boost_cell_type']))}</td>"
            f"<td>{html.escape(str(row['boost_sub_cell_type']))}</td>"
            f"<td>{html.escape(str(row['boost_confidence']))}</td>"
            f"<td>{report_link}</td>"
            f"<td>{html.escape(str(row['error']))}</td>"
            "</tr>"
        )
    html_path = auto_dir / "auto_report.html"
    html_doc = """<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CASSIA Boost Auto Report</title>
  <style>
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 32px; color: #1f2937; }
    table { border-collapse: collapse; width: 100%; font-size: 14px; }
    th, td { border: 1px solid #d1d5db; padding: 8px; vertical-align: top; }
    th { background: #f3f4f6; text-align: left; }
    code { background: #f3f4f6; padding: 2px 4px; border-radius: 4px; }
  </style>
</head>
<body>
  <h1>CASSIA Boost Auto Report</h1>
  <p>Status: <code>__STATUS__</code>. Selected clusters: __SELECTED__.</p>
  <table>
    <thead>
      <tr>
        <th>Cluster</th><th>Status</th><th>Score</th><th>Original Confidence</th>
        <th>Reasons</th><th>Boost Cell Type</th><th>Boost Subtype</th>
        <th>Boost Confidence</th><th>Report</th><th>Error</th>
      </tr>
    </thead>
    <tbody>
      __ROWS__
    </tbody>
  </table>
</body>
</html>
"""
    html_doc = (
        html_doc
        .replace("__STATUS__", html.escape(str(manifest.get("status", ""))))
        .replace("__SELECTED__", str(len(manifest.get("selected_candidates", []))))
        .replace("__ROWS__", "\n".join(html_rows))
    )
    html_path.write_text(
        html_doc,
        encoding="utf-8",
    )
    return {
        "csv": str(csv_path),
        "markdown": str(md_path),
        "html": str(html_path),
    }


def _read_final_result(boost_dir: Path) -> Optional[Dict[str, Any]]:
    final_path = boost_dir / "final.json"
    if final_path.exists():
        return _read_json(final_path, {})
    return None


def _auto_additional_task(candidate: Dict[str, Any], user_task: Optional[str]) -> str:
    reasons = "; ".join(candidate.get("reasons", []))
    auto_task = (
        f"This cluster was auto-selected for annotation boost because: {reasons}. "
        "Stress-test the original annotation, query decisive positive and negative markers when useful, "
        "and clearly state whether the original annotation should change."
    )
    if user_task:
        return f"{auto_task}\n\nUser additional task: {user_task}"
    return auto_task


def run_boost_auto(args: Any) -> int:
    """Auto-select clusters from a run folder and run boost on each selected cluster."""
    run_dir = Path(args.run)
    auto_dir = Path(args.out) if args.out else run_dir / "boost" / "_auto"
    auto_dir.mkdir(parents=True, exist_ok=True)

    target_terms = parse_target_terms(args.target_lineage)
    confidence_levels = tuple(args.confidence or DEFAULT_AUTO_CONFIDENCE)
    selected = select_boost_candidates(
        run_dir=run_dir,
        confidence_levels=confidence_levels,
        target_terms=target_terms,
        max_clusters=args.max_clusters,
        select_all=args.all,
        only_low_confidence=args.only_low_confidence,
    )

    manifest: Dict[str, Any] = {
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "status": "planned",
        "run_dir": str(run_dir.resolve()),
        "auto_dir": str(auto_dir.resolve()),
        "marker_table": str(Path(args.markers).resolve()),
        "selected_candidates": selected,
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items() if key != "func"},
        "results": [],
    }
    _write_json(auto_dir / "boost_auto_manifest.json", manifest)

    if not selected:
        manifest["status"] = "no-candidates"
        manifest["updated_at"] = utc_now()
        reports = write_auto_report(auto_dir, manifest, [])
        manifest["reports"] = reports
        _write_json(auto_dir / "boost_auto_manifest.json", manifest)
        print(f"No boost candidates found. Wrote {auto_dir / 'boost_auto_manifest.json'}")
        return 0

    if args.plan_only:
        manifest["status"] = "plan-only"
        manifest["updated_at"] = utc_now()
        plan_results = [{"candidate": candidate, "status": "planned"} for candidate in selected]
        reports = write_auto_report(auto_dir, manifest, plan_results)
        manifest["results"] = plan_results
        manifest["reports"] = reports
        _write_json(auto_dir / "boost_auto_manifest.json", manifest)
        print(f"Wrote {auto_dir / 'boost_auto_manifest.json'}")
        print(f"Wrote {reports['csv']}")
        print(f"Wrote {reports['html']}")
        return 0

    run_results: List[Dict[str, Any]] = []
    failed = False
    for candidate in selected:
        cluster = str(candidate["cluster"])
        cluster_slug = slugify(cluster)
        cluster_out = (auto_dir / "clusters" / cluster_slug) if args.out else run_dir / "boost" / cluster_slug
        result_record: Dict[str, Any] = {
            "candidate": candidate,
            "boost_dir": str(cluster_out.resolve()),
            "status": "pending",
        }

        if (cluster_out / "final.json").exists() and not args.force and not args.dry_run:
            final_result = _read_final_result(cluster_out)
            html_report = cluster_out / "summary.html"
            result_record.update({
                "status": "skipped",
                "final_result": final_result,
                "html_report": str(html_report.resolve()) if html_report.exists() else "",
            })
            run_results.append(result_record)
            continue

        boost_args = SimpleNamespace(**vars(args))
        boost_args.cluster = cluster
        boost_args.out = str(cluster_out)
        boost_args.additional_task = _auto_additional_task(candidate, args.additional_task)
        boost_args.func = None
        try:
            code = run_boost(boost_args)
            final_result = _read_final_result(cluster_out)
            status = "dry-run" if args.dry_run else ("completed" if code == 0 else "failed")
            html_report = cluster_out / "summary.html"
            result_record.update({
                "status": status,
                "final_result": final_result,
                "html_report": str(html_report.resolve()) if html_report.exists() else "",
            })
            if code != 0:
                failed = True
        except Exception as exc:
            failed = True
            result_record.update({
                "status": "failed",
                "error": str(exc),
            })
        run_results.append(result_record)
        manifest["results"] = run_results
        manifest["updated_at"] = utc_now()
        _write_json(auto_dir / "boost_auto_manifest.json", manifest)
        if failed and args.fail_fast:
            break

    manifest["status"] = "dry-run" if args.dry_run else ("failed" if failed else "completed")
    manifest["updated_at"] = utc_now()
    manifest["results"] = run_results
    reports = write_auto_report(auto_dir, manifest, run_results)
    manifest["reports"] = reports
    _write_json(auto_dir / "boost_auto_manifest.json", manifest)
    print(f"Wrote {auto_dir / 'boost_auto_manifest.json'}")
    print(f"Wrote {reports['csv']}")
    print(f"Wrote {reports['markdown']}")
    print(f"Wrote {reports['html']}")
    return 1 if failed else 0


def run_boost(args: Any) -> int:
    """Run the CLI annotation boost agent loop for one cluster."""
    started_monotonic = time.monotonic()
    if args.iterations < 1:
        raise ValueError("--iterations must be at least 1")
    if args.n_genes < 1:
        raise ValueError("--n-genes must be at least 1")
    if args.max_genes_per_round is not None and args.max_genes_per_round < 1:
        raise ValueError("--max-genes-per-round must be at least 1 when provided")
    run_dir = Path(args.run)
    mode = getattr(args, "mode", "review")
    cluster_slug = slugify(args.cluster)
    boost_dir = Path(args.out) if args.out else run_dir / "boost" / cluster_slug
    prompts_dir = boost_dir / "prompts"
    raw_dir = boost_dir / "raw"
    queries_dir = boost_dir / "queries"
    for directory in (prompts_dir, raw_dir, queries_dir):
        directory.mkdir(parents=True, exist_ok=True)

    if mode == "fused":
        annotation_context = {
            "source": None,
            "cluster_id": args.cluster,
            "annotation": {
                "mode": "fused_primary_annotation",
                "prior_annotation": None,
            },
        }
    else:
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
        "mode": mode,
        "result_schema_version": ANNOTATION_SCHEMA_VERSION,
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
        mode=mode,
        prompt_variant=getattr(args, "fused_prompt_version", "v2-compact"),
        candidate_count=getattr(args, "candidate_count", 5),
    )

    if args.dry_run:
        (prompts_dir / "round_001.md").write_text(prompt, encoding="utf-8")
        manifest["status"] = "dry-run"
        manifest["updated_at"] = utc_now()
        _write_json(boost_dir / "boost_manifest.json", manifest)
        return 0

    backend = AgentCLIBackend(
        args.backend,
        command_template=args.command_template,
        timeout_seconds=args.timeout,
        model=getattr(args, "model", None),
        agent_mode="ask" if args.backend == "cursor-agent" else None,
    )
    backend_cwd = boost_dir
    if args.backend != "shell":
        backend_cwd = Path(tempfile.gettempdir()) / f"cassia_boost_cli_{cluster_slug}"
        backend_cwd.mkdir(parents=True, exist_ok=True)
    messages: List[Dict[str, str]] = []
    errors: List[str] = []
    final_result: Optional[Dict[str, Any]] = None
    latest_query_text: Optional[str] = None
    query_frames: List[pd.DataFrame] = []
    checked_seen = set()
    initial_candidates: List[Dict[str, Any]] = []

    for round_idx in range(1, args.iterations + 1):
        if round_idx > 1:
            prompt = build_boost_followup_prompt(
                transcript=_append_transcript(messages),
                query_text=latest_query_text,
                is_final_round=round_idx == args.iterations,
                prompt_variant=getattr(args, "fused_prompt_version", "v2-compact"),
            )

        prompt_path = prompts_dir / f"round_{round_idx:03d}.md"
        raw_path = raw_dir / f"round_{round_idx:03d}.txt"
        prompt_path.write_text(prompt, encoding="utf-8")
        messages.append({"role": "user", "content": prompt})

        try:
            raw = backend.run(
                prompt,
                prompt_path,
                backend_cwd,
                {
                    "input": str(Path(args.markers).resolve()),
                    "out": str(boost_dir.resolve()),
                    "cluster": args.cluster,
                    "agent_output_file": str(raw_path),
                },
            )
            _record_agent_run_metadata(manifest, backend)
            raw_path.write_text(raw, encoding="utf-8")
            messages.append({"role": "assistant", "content": raw})

            if (
                getattr(args, "fused_prompt_version", "v2-compact") == "candidate"
                and not initial_candidates
            ):
                initial_candidates = extract_candidate_set(raw)
                if initial_candidates:
                    _write_json(
                        boost_dir / "candidate_set.json",
                        {"candidates": initial_candidates},
                    )
                    manifest["initial_candidates"] = initial_candidates

            genes = extract_check_genes(raw, max_genes=args.max_genes_per_round)
            genes = [gene for gene in genes if gene.upper() not in checked_seen]
            if genes:
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
                continue

            try:
                parsed = extract_json_object(raw)
                candidate_result = normalize_boost_result(parsed)
                if mode == "fused" and not query_frames:
                    final_result = None
                    latest_query_text = (
                        "Fused mode requires at least one marker-query round before finalization. "
                        "Request a targeted <check_genes> panel now."
                    )
                    continue
                final_result = candidate_result
                break
            except Exception:
                final_result = None
                latest_query_text = "No new <check_genes> request was found. Please return final JSON or request new genes."
        except Exception as exc:
            errors.append(str(exc))
            break

    if final_result is None and not errors:
        final_prompt = build_boost_followup_prompt(
            transcript=_append_transcript(messages),
            query_text=latest_query_text,
            is_final_round=True,
            prompt_variant=getattr(args, "fused_prompt_version", "v2-compact"),
        )
        prompt_path = prompts_dir / "final.md"
        raw_path = raw_dir / "final.txt"
        prompt_path.write_text(final_prompt, encoding="utf-8")
        messages.append({"role": "user", "content": final_prompt})
        try:
            raw = backend.run(
                final_prompt,
                prompt_path,
                backend_cwd,
                {
                    "input": str(Path(args.markers).resolve()),
                    "out": str(boost_dir.resolve()),
                    "cluster": args.cluster,
                    "agent_output_file": str(raw_path),
                },
            )
            _record_agent_run_metadata(manifest, backend)
            raw_path.write_text(raw, encoding="utf-8")
            messages.append({"role": "assistant", "content": raw})
            if extract_check_genes(raw, max_genes=args.max_genes_per_round):
                raise ValueError("Agent requested additional genes after the final query round")
            candidate_result = normalize_boost_result(extract_json_object(raw))
            if mode == "fused" and not query_frames:
                raise ValueError("Fused mode cannot finalize without at least one marker-query round")
            final_result = candidate_result
        except Exception as exc:
            errors.append(str(exc))

    transcript = _append_transcript(messages)
    (boost_dir / "transcript.md").write_text(transcript + "\n", encoding="utf-8")
    if final_result:
        final_result["cluster_id"] = str(args.cluster)
        final_result["checked_genes"] = list(dict.fromkeys([
            *manifest.get("checked_genes", []),
            *final_result.get("checked_genes", []),
        ]))
        _write_json(boost_dir / "final.json", final_result)

    manifest["status"] = "completed" if final_result and not errors else "failed"
    manifest["updated_at"] = utc_now()
    manifest["checked_genes"] = list(checked_seen)
    manifest["initial_candidates"] = initial_candidates
    manifest["final_json"] = str(boost_dir / "final.json") if final_result else None
    manifest["errors"] = errors
    manifest["execution_time"] = round(time.monotonic() - started_monotonic, 3)
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
    manifest["outputs"] = {
        "final_json": str(boost_dir / "final.json") if final_result else None,
        "markdown_report": str(report_path),
        "html_report": str(html_path),
        "summary_tags": str(boost_dir / "summary_tags.txt"),
        "transcript": str(boost_dir / "transcript.md"),
    }
    _write_json(boost_dir / "boost_manifest.json", manifest)
    print(f"Wrote {boost_dir / 'boost_manifest.json'}")
    print(f"Wrote {boost_dir / 'transcript.md'}")
    print(f"Wrote {report_path}")
    print(f"Wrote {html_path}")
    return 0 if final_result and not errors else 1
