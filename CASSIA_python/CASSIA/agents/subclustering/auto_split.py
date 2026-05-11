"""
Experimental auto-split orchestration for CASSIA subclustering.

This module keeps the default ``runCASSIA_subclusters`` behavior unchanged and
exposes an opt-in wrapper that:

1. Builds compact marker summaries for all input subclusters.
2. Uses an LLM to group clusters into smaller biologically coherent batches.
3. Validates that every input cluster appears exactly once.
4. Runs normal CASSIA subclustering per group.
5. Merges outputs and retries missing rows as single-cluster jobs.

The feature is intentionally conservative: if the LLM split plan is invalid,
it falls back to deterministic marker-module grouping.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

try:
    from CASSIA.core.llm_utils import call_llm
except ImportError:
    try:
        from ...core.llm_utils import call_llm
    except ImportError:
        from llm_utils import call_llm


DEFAULT_MAX_GROUP_SIZE = 5


MARKER_MODULES: Dict[str, Dict[str, Any]] = {
    "cd8_exhaustion_cytotoxic": {
        "markers": {
            "CD8A", "CD8B", "GZMK", "GZMA", "GZMB", "GZMH", "PRF1",
            "NKG7", "GNLY", "CCL5", "IFNG", "PDCD1", "LAG3", "HAVCR2",
            "TOX", "CXCL13", "ENTPD1", "TIGIT", "TCF7", "SLAMF6",
        },
        "hint": "CD8 T-cell cytotoxic, effector-memory, interferon-stimulated, and exhaustion continuum",
        "reference_hint": "CD8 T cell",
    },
    "cd4_helper_treg": {
        "markers": {
            "CD4", "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "FOXP3",
            "IL2RA", "CTLA4", "TIGIT", "IKZF2", "CCR8", "TNFRSF18",
            "TNFRSF4", "CXCR5", "BCL6", "ICOS", "IL21", "CD40LG",
            "TOX2", "TBX21", "CXCR3", "RORC", "IL17A", "IL17F", "CCR6",
            "GATA3", "IL4", "IL5", "IL13",
        },
        "hint": "CD4 T-cell helper, Treg, Tfh, Th17/Th1/Th2, and checkpoint states",
        "reference_hint": "CD4 T cell",
    },
    "innate_like_t": {
        "markers": {
            "KLRB1", "SLC4A10", "TRAV1-2", "ZBTB16", "TRDC", "TRGC1",
            "TRGC2", "TRDV1", "TRDV2", "KLRC1", "KLRC2", "KLRC3",
            "KLRD1", "KLRF1", "NCAM1",
        },
        "hint": "innate-like T-cell, MAIT, gamma-delta, NKT-like, and NK-like lymphocyte programs",
        "reference_hint": "innate-like T cell",
    },
    "state_module": {
        "markers": {
            "MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", "PCLAF",
            "CENPF", "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFITM1",
            "IFITM3", "MX1", "OAS1", "OAS2", "RSAD2", "HSPA1A",
            "HSPA1B", "HSPA6", "DNAJB1", "FOS", "JUN", "CD69", "ITGAE",
            "CXCR6",
        },
        "hint": "T-cell state modules such as proliferation, interferon stimulation, tissue residency, or stress",
        "reference_hint": "T cell",
    },
    "myeloid_macrophage": {
        "markers": {
            "LYZ", "LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC",
            "APOE", "APOC1", "SPP1", "TREM2", "FOLR2", "SELENOP",
            "SLC40A1", "IL1B", "TNF", "CXCL9", "CXCL10",
        },
        "hint": "myeloid, monocyte, macrophage, or TAM programs",
        "reference_hint": "macrophage",
    },
    "epithelial_contaminant": {
        "markers": {
            "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CLDN7",
            "CLDN8", "CDX2", "AGR2", "LGALS4",
        },
        "hint": "epithelial contaminant or non-parent epithelial program",
        "reference_hint": "epithelial contaminant",
    },
    "stromal_neural_contaminant": {
        "markers": {
            "COL1A1", "COL1A2", "COL3A1", "SPARC", "DCN", "LUM",
            "PDGFRA", "PLP1", "PMP22", "SOX2", "HAND2", "NRXN1",
            "CHL1", "LAMA4",
        },
        "hint": "stromal, neural, glial, Schwann-like, or other non-parent contaminant program",
        "reference_hint": "non-T contaminant",
    },
}


def runCASSIA_subclusters_auto_split(
    marker,
    major_cluster_info: str,
    output_name: str,
    *,
    model: Optional[str] = None,
    temperature: Optional[float] = None,
    provider: str = "openrouter",
    n_genes: int = 50,
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    additional_context: Optional[str] = None,
    use_reference: bool = True,
    reference_provider: Optional[str] = None,
    reference_model: Optional[str] = None,
    reference_cell_type_hint: Optional[str] = None,
    reference_depth: str = "detailed",
    reference_max_content_length: int = 5000,
    reference_max_context_length: int = 12000,
    split_provider: Optional[str] = None,
    split_model: Optional[str] = None,
    split_temperature: float = 0,
    max_group_size: int = DEFAULT_MAX_GROUP_SIZE,
    force_split: bool = True,
    retry_missing: bool = True,
) -> pd.DataFrame:
    """Run experimental LLM-assisted auto-split subclustering.

    Args mirror ``runCASSIA_subclusters`` where practical. The existing
    subclustering implementation is called once per planned group. The final
    merged CSV is written to ``{output_name}.csv``. Group-level CSV/HTML files
    are also retained with ``_<group_id>`` suffixes for inspection.
    """
    from .subclustering import _prepare_subcluster_marker_dataframe, runCASSIA_subclusters

    marker_df = _prepare_subcluster_marker_dataframe(marker, n_genes=n_genes)
    marker_df["cluster"] = marker_df["cluster"].astype(str)
    cluster_ids = marker_df["cluster"].tolist()
    if not force_split and len(cluster_ids) <= max_group_size:
        runCASSIA_subclusters(
            marker=marker,
            major_cluster_info=major_cluster_info,
            output_name=output_name,
            model=model,
            temperature=temperature,
            provider=provider,
            n_genes=n_genes,
            tissue=tissue,
            species=species,
            additional_context=additional_context,
            use_reference=use_reference,
            reference_provider=reference_provider,
            reference_model=reference_model,
            reference_cell_type_hint=reference_cell_type_hint,
            reference_depth=reference_depth,
            reference_max_content_length=reference_max_content_length,
            reference_max_context_length=reference_max_context_length,
        )
        return _read_result_csv(output_name)

    plan = _plan_split(
        marker_df=marker_df,
        major_cluster_info=major_cluster_info,
        provider=split_provider or reference_provider or provider,
        model=split_model or reference_model or model,
        temperature=split_temperature,
        max_group_size=max_group_size,
    )
    output_base = Path(output_name)
    plan_path = output_base.with_name(output_base.name + "_auto_split_plan.json")
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")

    results: List[pd.DataFrame] = []
    group_summaries: List[Dict[str, Any]] = []
    for group in plan["groups"]:
        result = _run_group(
            runCASSIA_subclusters=runCASSIA_subclusters,
            marker_df=marker_df,
            group=group,
            output_name=output_name,
            major_cluster_info=major_cluster_info,
            model=model,
            temperature=temperature,
            provider=provider,
            n_genes=n_genes,
            tissue=tissue,
            species=species,
            additional_context=additional_context,
            use_reference=use_reference,
            reference_provider=reference_provider,
            reference_model=reference_model,
            reference_cell_type_hint=reference_cell_type_hint,
            reference_depth=reference_depth,
            reference_max_content_length=reference_max_content_length,
            reference_max_context_length=reference_max_context_length,
        )
        results.append(result)
        group_summaries.append({
            "group_id": group.get("group_id"),
            "cluster_ids": group.get("cluster_ids", []),
            "output_rows": result["Result ID"].astype(str).tolist() if "Result ID" in result.columns else [],
        })

    merged = _merge_results(results)
    expected = set(cluster_ids)
    observed = set(merged["Result ID"].astype(str)) if "Result ID" in merged.columns else set()
    missing = sorted(expected - observed)

    retry_summaries: List[Dict[str, Any]] = []
    if retry_missing and missing:
        for missing_id in missing:
            retry_group = {
                "group_id": f"retry_missing_{_safe_id(missing_id)}",
                "cluster_ids": [missing_id],
                "group_hint": f"retry missing cluster {missing_id} from {major_cluster_info}",
                "reference_cell_type_hint": reference_cell_type_hint or major_cluster_info,
                "reason": "Single-cluster retry because the grouped run missed this Result ID.",
            }
            result = _run_group(
                runCASSIA_subclusters=runCASSIA_subclusters,
                marker_df=marker_df,
                group=retry_group,
                output_name=output_name,
                major_cluster_info=major_cluster_info,
                model=model,
                temperature=temperature,
                provider=provider,
                n_genes=n_genes,
                tissue=tissue,
                species=species,
                additional_context=additional_context,
                use_reference=use_reference,
                reference_provider=reference_provider,
                reference_model=reference_model,
                reference_cell_type_hint=reference_cell_type_hint,
                reference_depth=reference_depth,
                reference_max_content_length=reference_max_content_length,
                reference_max_context_length=reference_max_context_length,
            )
            results.append(result)
            retry_summaries.append({
                "cluster_id": missing_id,
                "output_rows": result["Result ID"].astype(str).tolist() if "Result ID" in result.columns else [],
            })
        merged = _merge_results(results)

    final_missing = sorted(expected - set(merged["Result ID"].astype(str))) if "Result ID" in merged.columns else sorted(expected)
    merged.to_csv(f"{output_name}.csv", index=False)

    summary = {
        "mode": "auto_split",
        "input_cluster_ids": cluster_ids,
        "planner_status": plan.get("planner_status"),
        "planner_validation": plan.get("planner_validation"),
        "groups": group_summaries,
        "missing_after_group_runs": missing,
        "retry_missing": retry_missing,
        "retry_summaries": retry_summaries,
        "final_missing": final_missing,
        "final_output": f"{output_name}.csv",
    }
    summary_path = output_base.with_name(output_base.name + "_auto_split_summary.json")
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return merged


def _plan_split(
    marker_df: pd.DataFrame,
    major_cluster_info: str,
    provider: str,
    model: Optional[str],
    temperature: float,
    max_group_size: int,
) -> Dict[str, Any]:
    prompt = _build_split_prompt(marker_df, major_cluster_info, max_group_size)
    try:
        response = call_llm(
            prompt=prompt,
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=1400,
        )
        plan = _parse_json_response(response)
        ok, reason = _validate_split_plan(plan, marker_df, max_group_size)
        if ok:
            plan["planner_status"] = "llm_valid"
            plan["planner_validation"] = reason
            return plan
        fallback = _fallback_split_plan(marker_df, max_group_size)
        fallback["planner_status"] = "llm_invalid_fallback_used"
        fallback["planner_validation"] = reason
        fallback["llm_plan"] = plan
        return fallback
    except Exception as exc:
        fallback = _fallback_split_plan(marker_df, max_group_size)
        fallback["planner_status"] = "llm_failed_fallback_used"
        fallback["planner_validation"] = str(exc)
        return fallback


def _build_split_prompt(marker_df: pd.DataFrame, major_cluster_info: str, max_group_size: int) -> str:
    lines = []
    for _, row in marker_df.iterrows():
        cluster_id = str(row["cluster"])
        markers = _split_markers(row["markers"])[:30]
        tags = _module_scores(markers)[:3]
        tag_text = ", ".join(f"{name}:{score}" for name, score in tags if score > 0) or "no strong module"
        lines.append(f"- Cluster {cluster_id}: markers={', '.join(markers)}; rough_modules={tag_text}")

    return f"""You are planning an experimental CASSIA subclustering run.

You are NOT annotating final cell types. Your only task is to split input
subclusters into smaller biologically coherent groups for downstream
annotation.

Parent cluster context:
{major_cluster_info}

Rules:
- Every cluster ID must appear exactly once.
- Each group must contain at most {max_group_size} clusters.
- Keep confusing states together when they need direct comparison.
- Separate clear non-parent contaminants when useful.
- Use short, stable snake_case group IDs.
- Return JSON only.

Cluster marker summaries:
{chr(10).join(lines)}

Return this JSON schema:
{{
  "strategy": "short explanation",
  "groups": [
    {{
      "group_id": "short_snake_case_id",
      "cluster_ids": ["cluster_id_1", "cluster_id_2"],
      "group_hint": "biological context for downstream annotation",
      "reference_cell_type_hint": "CD8 T cell|CD4 T cell|macrophage|T cell|non-parent contaminant",
      "reason": "why these clusters should be compared together"
    }}
  ]
}}"""


def _fallback_split_plan(marker_df: pd.DataFrame, max_group_size: int) -> Dict[str, Any]:
    assignments: Dict[str, List[str]] = {}
    for _, row in marker_df.iterrows():
        cluster_id = str(row["cluster"])
        markers = _split_markers(row["markers"])[:40]
        best = _best_module(markers)
        assignments.setdefault(best, []).append(cluster_id)

    groups = []
    for module_name, ids in assignments.items():
        module = MARKER_MODULES.get(module_name, {})
        for idx, chunk in enumerate(_chunks(ids, max_group_size)):
            suffix = f"_{idx + 1}" if len(ids) > max_group_size else ""
            groups.append({
                "group_id": f"{module_name}{suffix}",
                "cluster_ids": chunk,
                "group_hint": module.get("hint", "mixed marker-program group"),
                "reference_cell_type_hint": module.get("reference_hint", "cell subtype"),
                "reason": "Deterministic fallback grouping by marker-module overlap.",
            })
    return {
        "strategy": "deterministic marker-module fallback grouping",
        "groups": groups,
    }


def _run_group(
    *,
    runCASSIA_subclusters,
    marker_df: pd.DataFrame,
    group: Dict[str, Any],
    output_name: str,
    major_cluster_info: str,
    model: Optional[str],
    temperature: Optional[float],
    provider: str,
    n_genes: int,
    tissue: Optional[str],
    species: Optional[str],
    additional_context: Optional[str],
    use_reference: bool,
    reference_provider: Optional[str],
    reference_model: Optional[str],
    reference_cell_type_hint: Optional[str],
    reference_depth: str,
    reference_max_content_length: int,
    reference_max_context_length: int,
) -> pd.DataFrame:
    ids = [str(item) for item in group.get("cluster_ids", [])]
    group_marker_df = marker_df[marker_df["cluster"].astype(str).isin(ids)].copy()
    group_id = _safe_id(str(group.get("group_id", "group")))
    group_output = f"{output_name}_{group_id}"
    group_context = f"{major_cluster_info}\nAuto-split subgroup context: {group.get('group_hint', '')}"
    group_reference_hint = (
        group.get("reference_cell_type_hint")
        or reference_cell_type_hint
        or major_cluster_info
    )

    runCASSIA_subclusters(
        marker=group_marker_df,
        major_cluster_info=group_context,
        output_name=group_output,
        model=model,
        temperature=temperature,
        provider=provider,
        n_genes=n_genes,
        tissue=tissue,
        species=species,
        additional_context=additional_context,
        use_reference=use_reference,
        reference_provider=reference_provider,
        reference_model=reference_model,
        reference_cell_type_hint=group_reference_hint,
        reference_depth=reference_depth,
        reference_max_content_length=reference_max_content_length,
        reference_max_context_length=reference_max_context_length,
    )
    return _read_result_csv(group_output)


def _read_result_csv(output_name: str) -> pd.DataFrame:
    result = pd.read_csv(f"{output_name}.csv")
    if "Result ID" in result.columns:
        result["Result ID"] = result["Result ID"].astype(str)
    return result


def _merge_results(results: Sequence[pd.DataFrame]) -> pd.DataFrame:
    frames = [frame for frame in results if frame is not None and not frame.empty]
    if not frames:
        return pd.DataFrame()
    merged = pd.concat(frames, ignore_index=True)
    if "Result ID" in merged.columns:
        merged["Result ID"] = merged["Result ID"].astype(str)
        merged = merged.drop_duplicates(subset=["Result ID"], keep="last")
    return merged


def _validate_split_plan(plan: Dict[str, Any], marker_df: pd.DataFrame, max_group_size: int) -> Tuple[bool, str]:
    groups = plan.get("groups") if isinstance(plan, dict) else None
    if not isinstance(groups, list) or not groups:
        return False, "groups missing or empty"
    expected = set(marker_df["cluster"].astype(str).tolist())
    observed: List[str] = []
    for group in groups:
        ids = group.get("cluster_ids") if isinstance(group, dict) else None
        if not isinstance(ids, list) or not ids:
            return False, f"group {group.get('group_id') if isinstance(group, dict) else '?'} has no cluster_ids"
        if len(ids) > max_group_size:
            return False, f"group {group.get('group_id')} has more than {max_group_size} clusters"
        observed.extend(str(item) for item in ids)
    observed_set = set(observed)
    if observed_set != expected:
        missing = sorted(expected - observed_set)
        extra = sorted(observed_set - expected)
        return False, f"split mismatch missing={missing} extra={extra}"
    if len(observed) != len(observed_set):
        duplicates = sorted({item for item in observed if observed.count(item) > 1})
        return False, f"duplicate cluster ids: {duplicates}"
    return True, "ok"


def _parse_json_response(response: str) -> Dict[str, Any]:
    match = re.search(r"\{[\s\S]*\}", str(response or ""))
    if not match:
        return {}
    try:
        return json.loads(match.group())
    except json.JSONDecodeError:
        return {}


def _split_markers(value: Any) -> List[str]:
    if isinstance(value, list):
        raw = value
    else:
        raw = re.split(r"[,;|\n]+", str(value or ""))
    markers = []
    seen = set()
    for marker in raw:
        clean = str(marker).strip().strip("'\"`")
        if not clean:
            continue
        key = clean.upper()
        if key in seen:
            continue
        seen.add(key)
        markers.append(clean)
    return markers


def _module_scores(markers: Iterable[str]) -> List[Tuple[str, int]]:
    marker_set = {str(marker).upper() for marker in markers}
    scores = [
        (name, len(marker_set & set(module["markers"])))
        for name, module in MARKER_MODULES.items()
    ]
    return sorted(scores, key=lambda item: (-item[1], item[0]))


def _best_module(markers: Iterable[str]) -> str:
    scores = _module_scores(markers)
    if not scores or scores[0][1] == 0:
        return "mixed_other"
    return scores[0][0]


def _chunks(items: Sequence[str], size: int) -> List[List[str]]:
    size = max(1, int(size or DEFAULT_MAX_GROUP_SIZE))
    return [list(items[index:index + size]) for index in range(0, len(items), size)]


def _safe_id(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_").lower()
    return safe or "group"
