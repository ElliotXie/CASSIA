"""Marker-set computation for parent-cluster subclustering workflows."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


DEFAULT_SUBCLUSTER_MARKER_MODES = ("local", "global_relabel", "target_vs_background")


def _require_scanpy():
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError(
            "compute_subcluster_marker_sets requires scanpy for differential marker ranking. "
            "Install it with: pip install scanpy"
        ) from exc
    return sc


def _load_anndata(data: Any):
    if hasattr(data, "obs") and hasattr(data, "X") and hasattr(data, "var_names"):
        return data

    if isinstance(data, (str, Path)):
        path = Path(data)
        if not path.exists():
            raise FileNotFoundError(f"AnnData file not found: {path}")
        if path.suffix.lower() != ".h5ad":
            raise ValueError("Path input must point to a .h5ad file")
        sc = _require_scanpy()
        return sc.read_h5ad(path)

    raise TypeError("data must be an AnnData object or a path to a .h5ad file")


def _validate_modes(modes: Sequence[str]) -> Tuple[str, ...]:
    valid = set(DEFAULT_SUBCLUSTER_MARKER_MODES)
    cleaned: List[str] = []
    for mode in modes:
        mode = str(mode).strip()
        if mode not in valid:
            raise ValueError(f"Unknown marker mode '{mode}'. Expected one of {sorted(valid)}")
        if mode not in cleaned:
            cleaned.append(mode)
    if not cleaned:
        raise ValueError("At least one marker mode must be requested")
    return tuple(cleaned)


def _obs_labels(adata, column: str) -> pd.Series:
    if column not in adata.obs:
        raise ValueError(f"Column '{column}' was not found in adata.obs")
    return adata.obs[column].astype(str)


def _ordered_subclusters(adata, parent_mask: np.ndarray, subcluster_col: str) -> List[str]:
    if subcluster_col not in adata.obs:
        raise ValueError(f"Column '{subcluster_col}' was not found in adata.obs")

    series = adata.obs.loc[parent_mask, subcluster_col]
    if hasattr(series.dtype, "categories"):
        observed = set(series.dropna().astype(str))
        return [str(value) for value in series.dtype.categories if str(value) in observed]

    return [str(value) for value in pd.unique(series.dropna().astype(str))]


def _matrix_and_var_names(adata, use_raw: Optional[bool], layer: Optional[str]):
    if use_raw and layer:
        raise ValueError("use_raw and layer cannot both be set")
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True was requested, but adata.raw is not available")
        return adata.raw.X, list(adata.raw.var_names)
    if layer:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' was not found in adata.layers")
        return adata.layers[layer], list(adata.var_names)
    return adata.X, list(adata.var_names)


def _column_vector(matrix: Any, index: int) -> np.ndarray:
    values = matrix[:, index]
    if hasattr(values, "toarray"):
        values = values.toarray()
    elif hasattr(values, "A"):
        values = values.A
    values = np.asarray(values).reshape(-1)
    return values.astype(float, copy=False)


def _expression_stats(
    adata,
    genes: Sequence[str],
    group_mask: np.ndarray,
    reference_mask: np.ndarray,
    use_raw: Optional[bool],
    layer: Optional[str],
    expression_threshold: float,
) -> Dict[str, Dict[str, float]]:
    matrix, var_names = _matrix_and_var_names(adata, use_raw=use_raw, layer=layer)
    gene_to_index = {str(gene): idx for idx, gene in enumerate(var_names)}
    stats: Dict[str, Dict[str, float]] = {}

    for gene in genes:
        idx = gene_to_index.get(str(gene))
        if idx is None:
            continue
        expr = _column_vector(matrix, idx)
        group_expr = expr[group_mask]
        ref_expr = expr[reference_mask]
        stats[str(gene)] = {
            "pct.1": float(np.mean(group_expr > expression_threshold)) if len(group_expr) else np.nan,
            "pct.2": float(np.mean(ref_expr > expression_threshold)) if len(ref_expr) else np.nan,
            "mean.1": float(np.mean(group_expr)) if len(group_expr) else np.nan,
            "mean.2": float(np.mean(ref_expr)) if len(ref_expr) else np.nan,
        }
    return stats


def _normalize_rank_df(
    df: pd.DataFrame,
    mode: str,
    subcluster: str,
    group_label: str,
    reference_label: str,
    comparison: str,
    n_genes: int,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()

    column_map = {
        "names": "gene",
        "scores": "score",
        "logfoldchanges": "avg_log2FC",
        "pvals": "p_val",
        "pvals_adj": "p_val_adj",
    }
    work = df.rename(columns=column_map).head(n_genes).copy()
    if "gene" not in work.columns:
        raise ValueError("Scanpy rank_genes_groups output did not include gene names")

    work.insert(0, "mode", mode)
    work.insert(1, "subcluster", str(subcluster))
    work.insert(2, "comparison", comparison)
    work.insert(3, "rank", range(1, len(work) + 1))
    work["group_label"] = str(group_label)
    work["reference_label"] = str(reference_label)
    return work


def _apply_filters(
    df: pd.DataFrame,
    min_log2fc: Optional[float],
    max_p_val_adj: Optional[float],
    n_genes: int,
) -> pd.DataFrame:
    work = df.copy()
    if min_log2fc is not None and "avg_log2FC" in work.columns:
        work = work[pd.to_numeric(work["avg_log2FC"], errors="coerce") >= float(min_log2fc)]
    if max_p_val_adj is not None and "p_val_adj" in work.columns:
        work = work[pd.to_numeric(work["p_val_adj"], errors="coerce") <= float(max_p_val_adj)]
    work = work.groupby(["mode", "subcluster"], sort=False, observed=True).head(n_genes).copy()
    work["rank"] = work.groupby(["mode", "subcluster"], sort=False, observed=True).cumcount() + 1
    return work.reset_index(drop=True)


def _rank_groups(
    adata,
    groupby: str,
    groups: Sequence[str],
    reference: str,
    mode: str,
    group_to_subcluster: Mapping[str, str],
    comparison_by_group: Mapping[str, str],
    n_genes: int,
    method: str,
    use_raw: Optional[bool],
    layer: Optional[str],
    rankby_abs: bool,
    tie_correct: bool,
    expression_threshold: float,
    key_added: str,
) -> pd.DataFrame:
    sc = _require_scanpy()
    groups = [str(group) for group in groups]
    if not groups:
        return pd.DataFrame()

    adata.obs[groupby] = adata.obs[groupby].astype(str)
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups,
        reference=reference,
        n_genes=n_genes,
        method=method,
        use_raw=use_raw,
        layer=layer,
        rankby_abs=rankby_abs,
        tie_correct=tie_correct,
        pts=False,
        key_added=key_added,
    )

    rows = []
    labels = adata.obs[groupby].astype(str)
    for group in groups:
        group_df = sc.get.rank_genes_groups_df(adata, group=group, key=key_added)
        subcluster = group_to_subcluster[group]
        comparison = comparison_by_group[group]
        reference_label = "rest" if reference == "rest" else reference
        normalized = _normalize_rank_df(
            group_df,
            mode=mode,
            subcluster=subcluster,
            group_label=group,
            reference_label=reference_label,
            comparison=comparison,
            n_genes=n_genes,
        )
        if normalized.empty:
            continue

        group_mask = (labels == group).to_numpy()
        reference_mask = (labels != group).to_numpy() if reference == "rest" else (labels == reference).to_numpy()
        stats = _expression_stats(
            adata,
            genes=normalized["gene"].astype(str).tolist(),
            group_mask=group_mask,
            reference_mask=reference_mask,
            use_raw=use_raw,
            layer=layer,
            expression_threshold=expression_threshold,
        )
        for stat_col in ("pct.1", "pct.2", "mean.1", "mean.2"):
            normalized[stat_col] = normalized["gene"].astype(str).map(
                lambda gene: stats.get(gene, {}).get(stat_col, np.nan)
            )
        rows.append(normalized)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def _top_gene_text(df: pd.DataFrame, subcluster: str, mode: str, n_genes: int) -> str:
    if df.empty:
        return ""
    mask = (df["subcluster"].astype(str) == str(subcluster)) & (df["mode"].astype(str) == mode)
    genes = df.loc[mask].sort_values("rank")["gene"].astype(str).head(n_genes).tolist()
    return ", ".join(dict.fromkeys(genes))


def _build_combined_marker_table(
    mode_tables: Mapping[str, pd.DataFrame],
    subclusters: Sequence[str],
    n_genes: int,
) -> pd.DataFrame:
    all_rows = [df for df in mode_tables.values() if df is not None and not df.empty]
    long_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()

    rows = []
    for subcluster in subclusters:
        local = _top_gene_text(long_df, subcluster, "local", n_genes)
        global_relabel = _top_gene_text(long_df, subcluster, "global_relabel", n_genes)
        target_vs_background = _top_gene_text(long_df, subcluster, "target_vs_background", n_genes)
        sections = []
        if local:
            sections.append(f"local markers: {local}")
        if global_relabel:
            sections.append(f"global relabel markers: {global_relabel}")
        if target_vs_background:
            sections.append(f"target-vs-background markers: {target_vs_background}")
        rows.append({
            "cluster": str(subcluster),
            "markers": " | ".join(sections),
            "local_markers": local,
            "global_relabel_markers": global_relabel,
            "target_vs_background_markers": target_vs_background,
        })
    return pd.DataFrame(rows)


def _write_outputs(output_dir: Path, results: Dict[str, Any]) -> Dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: Dict[str, str] = {}
    for key in (
        "local_markers",
        "global_relabel_markers",
        "target_vs_background_markers",
        "combined_markers",
    ):
        value = results.get(key)
        if isinstance(value, pd.DataFrame):
            path = output_dir / f"{key}.csv"
            value.to_csv(path, index=False)
            written[key] = str(path)

    metadata_path = output_dir / "marker_set_metadata.json"
    metadata_path.write_text(
        json.dumps(results.get("metadata", {}), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    written["metadata"] = str(metadata_path)
    return written


def compute_subcluster_marker_sets(
    data: Any,
    parent_col: str,
    parent_label: str,
    subcluster_col: str,
    modes: Sequence[str] = DEFAULT_SUBCLUSTER_MARKER_MODES,
    n_genes: int = 50,
    method: str = "wilcoxon",
    output_dir: Optional[Any] = None,
    use_raw: Optional[bool] = None,
    layer: Optional[str] = None,
    rankby_abs: bool = False,
    tie_correct: bool = False,
    expression_threshold: float = 0.0,
    min_log2fc: Optional[float] = 0.25,
    max_p_val_adj: Optional[float] = None,
    key_prefix: str = "cassia_subcluster_marker_sets",
) -> Dict[str, Any]:
    """
    Compute three complementary marker views for subcluster annotation.

    Parameters
    ----------
    data
        AnnData object or path to a .h5ad file. The object must already contain
        broad parent labels and subcluster labels in ``adata.obs``.
    parent_col
        Observation column containing broad labels, e.g. ``"broad_annotation"``.
    parent_label
        Parent population to analyze, e.g. ``"Immune cell"``.
    subcluster_col
        Observation column containing subcluster IDs for cells inside the parent
        population.
    modes
        Marker modes to compute:
        ``local`` compares subclusters within the parent subset;
        ``global_relabel`` places all parent subclusters back into the full
        object and compares each subcluster to the rest;
        ``target_vs_background`` compares each target subcluster to cells
        outside the parent population, excluding sibling subclusters.
    n_genes
        Number of ranked marker genes per subcluster per mode.
    method
        Scanpy ``rank_genes_groups`` method, e.g. ``"wilcoxon"`` or ``"t-test"``.
    output_dir
        Optional directory where CSV files and metadata JSON are written.
    min_log2fc
        Optional positive log2 fold-change filter. Defaults to 0.25 to keep
        upregulated markers and drop genes that only rank because of ties.

    Returns
    -------
    dict
        Contains one DataFrame per marker mode, a ``combined_markers`` DataFrame
        suitable for CASSIA subclustering input, and metadata.
    """
    if n_genes <= 0:
        raise ValueError("n_genes must be a positive integer")

    modes = _validate_modes(modes)
    adata = _load_anndata(data)
    parent_labels = _obs_labels(adata, parent_col)
    parent_mask = (parent_labels == str(parent_label)).to_numpy()
    if not parent_mask.any():
        raise ValueError(f"No cells found where {parent_col} == {parent_label!r}")

    subclusters = _ordered_subclusters(adata, parent_mask, subcluster_col)
    if not subclusters:
        raise ValueError(f"No subcluster labels found in '{subcluster_col}' for parent label {parent_label!r}")

    mode_tables: Dict[str, pd.DataFrame] = {}

    if "local" in modes:
        local = adata[parent_mask].copy()
        local_group_col = "_cassia_local_subcluster"
        local.obs[local_group_col] = local.obs[subcluster_col].astype(str).values
        local_groups = [subcluster for subcluster in subclusters if (local.obs[local_group_col].astype(str) == subcluster).sum() > 0]
        if len(local_groups) >= 2:
            mode_tables["local"] = _rank_groups(
                local,
                groupby=local_group_col,
                groups=local_groups,
                reference="rest",
                mode="local",
                group_to_subcluster={group: group for group in local_groups},
                comparison_by_group={group: f"{group} vs sibling subclusters inside {parent_label}" for group in local_groups},
                n_genes=n_genes,
                method=method,
                use_raw=use_raw,
                layer=layer,
                rankby_abs=rankby_abs,
                tie_correct=tie_correct,
                expression_threshold=expression_threshold,
                key_added=f"{key_prefix}_local",
            )
        else:
            mode_tables["local"] = pd.DataFrame()

    if "global_relabel" in modes:
        global_data = adata.copy()
        global_group_col = "_cassia_global_relabel"
        global_labels = parent_labels.astype(str).copy()
        subcluster_labels = adata.obs[subcluster_col].astype(str)
        target_groups = []
        for subcluster in subclusters:
            target_label = f"{parent_label}::{subcluster}"
            target_mask = parent_mask & (subcluster_labels == str(subcluster)).to_numpy()
            global_labels.loc[target_mask] = target_label
            target_groups.append(target_label)
        global_data.obs[global_group_col] = global_labels.values
        mode_tables["global_relabel"] = _rank_groups(
            global_data,
            groupby=global_group_col,
            groups=target_groups,
            reference="rest",
            mode="global_relabel",
            group_to_subcluster={f"{parent_label}::{subcluster}": subcluster for subcluster in subclusters},
            comparison_by_group={
                f"{parent_label}::{subcluster}": f"{subcluster} relabeled into full object vs rest"
                for subcluster in subclusters
            },
            n_genes=n_genes,
            method=method,
            use_raw=use_raw,
            layer=layer,
            rankby_abs=rankby_abs,
            tie_correct=tie_correct,
            expression_threshold=expression_threshold,
            key_added=f"{key_prefix}_global_relabel",
        )

    if "target_vs_background" in modes:
        target_rows = []
        subcluster_labels = adata.obs[subcluster_col].astype(str)
        outside_parent_mask = ~parent_mask
        if not outside_parent_mask.any():
            mode_tables["target_vs_background"] = pd.DataFrame()
        else:
            for subcluster in subclusters:
                target_mask = parent_mask & (subcluster_labels == str(subcluster)).to_numpy()
                subset_mask = target_mask | outside_parent_mask
                target_data = adata[subset_mask].copy()
                target_group_col = "_cassia_target_vs_background"
                labels = np.where(target_mask[subset_mask], "target", "outside_parent_background")
                target_data.obs[target_group_col] = labels
                table = _rank_groups(
                    target_data,
                    groupby=target_group_col,
                    groups=["target"],
                    reference="outside_parent_background",
                    mode="target_vs_background",
                    group_to_subcluster={"target": subcluster},
                    comparison_by_group={
                        "target": f"{subcluster} vs cells outside {parent_label}; sibling subclusters excluded"
                    },
                    n_genes=n_genes,
                    method=method,
                    use_raw=use_raw,
                    layer=layer,
                    rankby_abs=rankby_abs,
                    tie_correct=tie_correct,
                    expression_threshold=expression_threshold,
                    key_added=f"{key_prefix}_target_{str(subcluster).replace(' ', '_')}",
                )
                target_rows.append(table)
            mode_tables["target_vs_background"] = (
                pd.concat(target_rows, ignore_index=True) if target_rows else pd.DataFrame()
            )

    for mode, table in list(mode_tables.items()):
        if not table.empty:
            mode_tables[mode] = _apply_filters(table, min_log2fc=min_log2fc, max_p_val_adj=max_p_val_adj, n_genes=n_genes)

    combined_markers = _build_combined_marker_table(mode_tables, subclusters=subclusters, n_genes=n_genes)

    metadata = {
        "parent_col": parent_col,
        "parent_label": str(parent_label),
        "subcluster_col": subcluster_col,
        "subclusters": subclusters,
        "modes": list(modes),
        "n_genes": n_genes,
        "method": method,
        "use_raw": use_raw,
        "layer": layer,
        "target_vs_background_reference": "cells outside parent_label; sibling subclusters excluded",
        "n_parent_cells": int(parent_mask.sum()),
        "n_total_cells": int(adata.n_obs),
    }

    results: Dict[str, Any] = {
        "local_markers": mode_tables.get("local", pd.DataFrame()),
        "global_relabel_markers": mode_tables.get("global_relabel", pd.DataFrame()),
        "target_vs_background_markers": mode_tables.get("target_vs_background", pd.DataFrame()),
        "combined_markers": combined_markers,
        "metadata": metadata,
    }

    if output_dir is not None:
        results["files"] = _write_outputs(Path(output_dir), results)

    return results
