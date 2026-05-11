"""Input validation for CASSIA CLI marker tables."""

from __future__ import annotations

import contextlib
import io
import json
import re
import shlex
from pathlib import Path
from statistics import median
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from CASSIA.core.marker_utils import split_markers

from .runner import load_marker_clusters


CLUSTER_COLUMN_CANDIDATES = [
    "cluster",
    "group",
    "Cluster",
    "Cluster ID",
    "celltype",
    "cell_type",
    "Broad.cell.type",
    "Broad Cell Type",
    "seurat_clusters",
    "leiden",
    "louvain",
]

GENE_COLUMN_CANDIDATES = [
    "gene",
    "genes",
    "names",
    "Gene",
    "Gene Symbol",
    "marker",
    "markers",
    "Top.Markers",
    "Top Markers",
    "Marker List",
]

SEURAT_REQUIRED_COLUMNS = ["cluster", "gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2"]
SCANPY_REQUIRED_COLUMNS = ["group", "names", "logfoldchanges", "pvals_adj"]
VALID_RANKING_METHODS = ["avg_log2FC", "p_val_adj", "pct_diff", "Score"]


def _column_map(columns: Sequence[str]) -> Dict[str, str]:
    return {str(column).casefold(): str(column) for column in columns}


def _find_column(columns: Sequence[str], candidates: Sequence[str], override: Optional[str] = None) -> Optional[str]:
    mapping = _column_map(columns)
    if override:
        return mapping.get(str(override).casefold())
    for candidate in candidates:
        found = mapping.get(candidate.casefold())
        if found is not None:
            return found
    return None


def _looks_like_marker_list(series: pd.Series) -> bool:
    sample = series.dropna().astype(str).head(20)
    if sample.empty:
        return False
    return bool(sample.str.contains(r"[,;]|\s+[A-Za-z0-9_.-]+\s+[A-Za-z0-9_.-]+").any())


def _looks_ensembl(value: str) -> bool:
    return bool(re.match(r"ENS[A-Z]*G\d+(?:\.\d+)?$", value.strip(), flags=re.IGNORECASE))


def _detect_format(df: pd.DataFrame, cluster_column: Optional[str], gene_column: Optional[str]) -> str:
    columns = list(df.columns)
    if {"group", "names", "logfoldchanges"}.issubset(set(columns)):
        return "scanpy_long"
    if {"cluster", "gene", "avg_log2FC"}.issubset(set(columns)):
        return "seurat_long"
    if len(columns) == 2:
        return "preformatted_marker_list"
    if gene_column and gene_column in df.columns and _looks_like_marker_list(df[gene_column]):
        return "preformatted_marker_list"
    marker_col = _find_column(columns, ["markers", "Top.Markers", "Top Markers", "Marker List", "genes"])
    if marker_col and _looks_like_marker_list(df[marker_col]):
        return "preformatted_marker_list"
    if cluster_column and gene_column and cluster_column in df.columns and gene_column in df.columns:
        return "generic_long"
    if _find_column(columns, CLUSTER_COLUMN_CANDIDATES) and _find_column(columns, ["gene", "names", "Gene"]):
        return "generic_long"
    return "unknown"


def _detect_columns(
    df: pd.DataFrame,
    detected_format: str,
    cluster_column: Optional[str],
    gene_column: Optional[str],
) -> Tuple[Optional[str], Optional[str]]:
    columns = list(df.columns)
    if cluster_column:
        cluster_col = _find_column(columns, [], cluster_column)
    elif detected_format == "scanpy_long":
        cluster_col = _find_column(columns, ["group"])
    else:
        cluster_col = _find_column(columns, CLUSTER_COLUMN_CANDIDATES) or (columns[0] if columns else None)

    if gene_column:
        gene_col = _find_column(columns, [], gene_column)
    elif detected_format == "scanpy_long":
        gene_col = _find_column(columns, ["names"])
    elif detected_format in {"seurat_long", "generic_long"}:
        gene_col = _find_column(columns, ["gene", "names", "Gene"])
    else:
        gene_col = _find_column(columns, GENE_COLUMN_CANDIDATES) or (columns[1] if len(columns) > 1 else None)
    return cluster_col, gene_col


def _required_columns_for_format(detected_format: str, ranking_method: str) -> List[str]:
    if detected_format == "scanpy_long":
        required = list(SCANPY_REQUIRED_COLUMNS)
        if ranking_method == "pct_diff":
            required.extend(["pct.1", "pct.2"])
        if ranking_method == "Score":
            required.append("Score")
        return required
    if detected_format in {"seurat_long", "generic_long"}:
        required = list(SEURAT_REQUIRED_COLUMNS)
        if ranking_method == "Score":
            required.append("Score")
        return required
    return []


def _summarize_marker_counts(clusters: Sequence[Any]) -> Dict[str, Any]:
    counts = [len(cluster.markers) for cluster in clusters]
    if not counts:
        return {
            "cluster_count": 0,
            "total_markers": 0,
            "min_markers": 0,
            "median_markers": 0,
            "max_markers": 0,
            "low_marker_clusters": [],
        }
    low = [str(cluster.cluster_id) for cluster in clusters if len(cluster.markers) < 5]
    return {
        "cluster_count": len(clusters),
        "total_markers": sum(counts),
        "min_markers": min(counts),
        "median_markers": median(counts),
        "max_markers": max(counts),
        "low_marker_clusters": low[:20],
    }


def _recommended_command(path: Path, cluster_col: Optional[str], gene_col: Optional[str], args: Any) -> str:
    parts = [
        "cassia",
        "annotate",
        "-i",
        str(path),
        "--backend",
        getattr(args, "backend", None) or "codex-cli",
        "--tissue",
        getattr(args, "tissue", None) or "TISSUE",
        "--species",
        getattr(args, "species", None) or "human",
    ]
    if cluster_col:
        parts.extend(["--celltype-column", cluster_col])
    if gene_col:
        parts.extend(["--gene-column", gene_col])
    parts.extend(["--n-genes", str(getattr(args, "n_genes", 50))])
    parts.extend(["--out", "runs/cassia_run"])
    return " ".join(shlex.quote(str(part)) for part in parts)


def validate_marker_input(args: Any) -> Dict[str, Any]:
    """Validate a marker CSV and return structured diagnostics."""
    input_path = Path(args.input).expanduser()
    errors: List[str] = []
    warnings: List[str] = []
    info: List[str] = []
    payload: Dict[str, Any] = {
        "input": str(input_path),
        "status": "error",
        "format": "unknown",
        "row_count": 0,
        "column_count": 0,
        "columns": [],
        "detected_columns": {
            "cluster_column": None,
            "gene_column": None,
            "ranking_method": args.ranking_method,
        },
        "summary": {},
        "errors": errors,
        "warnings": warnings,
        "info": info,
        "recommended_command": "",
    }

    if args.ranking_method not in VALID_RANKING_METHODS:
        errors.append(f"Invalid ranking method '{args.ranking_method}'. Expected one of: {', '.join(VALID_RANKING_METHODS)}")
        return payload
    if args.n_genes <= 0:
        errors.append("--n-genes must be greater than 0")
        return payload
    if not input_path.exists():
        errors.append(f"Input file does not exist: {input_path}")
        return payload
    if not input_path.is_file():
        errors.append(f"Input path is not a file: {input_path}")
        return payload

    try:
        df = pd.read_csv(input_path)
    except Exception as exc:
        errors.append(f"Could not read CSV: {exc}")
        return payload

    unnamed_cols = [column for column in df.columns if str(column).startswith("Unnamed:")]
    if unnamed_cols:
        df = df.drop(columns=unnamed_cols)
        info.append(f"Ignored index-like columns: {', '.join(map(str, unnamed_cols))}")

    payload["row_count"] = int(len(df))
    payload["column_count"] = int(len(df.columns))
    payload["columns"] = [str(column) for column in df.columns]

    if df.empty:
        errors.append("Input CSV has no rows")
        return payload
    if len(df.columns) < 2:
        errors.append("Input CSV must have at least two usable columns")
        return payload

    override_cluster = args.celltype_column
    override_gene = args.gene_column
    missing_overrides = [
        name
        for name, value in (("celltype column", override_cluster), ("gene column", override_gene))
        if value and _find_column(df.columns, [], value) is None
    ]
    if missing_overrides:
        errors.append(
            "Missing requested {items}. Available columns: {columns}".format(
                items=", ".join(missing_overrides),
                columns=", ".join(map(str, df.columns)),
            )
        )

    detected_format = _detect_format(df, override_cluster, override_gene)
    cluster_col, gene_col = _detect_columns(df, detected_format, override_cluster, override_gene)
    payload["format"] = detected_format
    payload["detected_columns"]["cluster_column"] = cluster_col
    payload["detected_columns"]["gene_column"] = gene_col

    if detected_format == "unknown":
        errors.append(
            "Could not determine whether this is a marker-list table or a long differential-expression table. "
            "Use --celltype-column and --gene-column if the columns are non-standard."
        )
    if not cluster_col:
        errors.append("Could not detect cluster/cell type column")
    if not gene_col:
        errors.append("Could not detect gene/marker column")

    if detected_format in {"seurat_long", "scanpy_long", "generic_long"}:
        required = _required_columns_for_format(detected_format, args.ranking_method)
        missing = [column for column in required if column not in df.columns]
        if missing:
            errors.append(
                "Long differential-expression tables need columns for CASSIA ranking/filtering. "
                f"Missing: {', '.join(missing)}"
            )
    elif detected_format == "preformatted_marker_list" and gene_col:
        empty_marker_rows = int(df[gene_col].isna().sum())
        if empty_marker_rows:
            warnings.append(f"{empty_marker_rows} row(s) have empty marker lists")

    if cluster_col and cluster_col in df.columns:
        empty_clusters = int(df[cluster_col].isna().sum())
        if empty_clusters:
            warnings.append(f"{empty_clusters} row(s) have empty cluster IDs")
        duplicates = int(df[cluster_col].astype(str).duplicated().sum())
        if duplicates and detected_format == "preformatted_marker_list":
            warnings.append(f"{duplicates} duplicated cluster ID row(s) found in marker-list format")

    if gene_col and gene_col in df.columns:
        sample_genes: List[str] = []
        for value in df[gene_col].dropna().astype(str).head(50):
            sample_genes.extend(split_markers(value)[:5])
        if sample_genes:
            ensembl_fraction = sum(1 for gene in sample_genes if _looks_ensembl(gene)) / len(sample_genes)
            if ensembl_fraction >= 0.5:
                warnings.append("Most sampled markers look like Ensembl gene IDs; CASSIA prompts work best with gene symbols")

    clusters = []
    if not errors:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                clusters = load_marker_clusters(
                    input_path,
                    n_genes=args.n_genes,
                    celltype_column=cluster_col,
                    gene_column=gene_col,
                    ranking_method=args.ranking_method,
                    ascending=args.ascending,
                    limit=args.limit,
                )
        except Exception as exc:
            errors.append(f"CASSIA could not build marker lists from this input: {exc}")

    if clusters:
        summary = _summarize_marker_counts(clusters)
        payload["summary"] = summary
        if summary["low_marker_clusters"]:
            shown = ", ".join(summary["low_marker_clusters"][:5])
            warnings.append(f"{len(summary['low_marker_clusters'])} cluster(s) have fewer than 5 markers after filtering: {shown}")
        payload["recommended_command"] = _recommended_command(input_path, cluster_col, gene_col, args)
    else:
        payload["summary"] = _summarize_marker_counts([])

    if errors:
        payload["status"] = "error"
    elif warnings:
        payload["status"] = "warning"
    else:
        payload["status"] = "ok"
    return payload


def _print_human_report(payload: Dict[str, Any]) -> None:
    status_label = {
        "ok": "OK",
        "warning": "WARNING",
        "error": "ERROR",
    }.get(payload["status"], payload["status"].upper())
    print(f"CASSIA input validation: {payload['input']}")
    print(f"Status: {status_label}")
    print("")
    print(f"Detected format: {payload['format']}")
    print(f"Rows: {payload['row_count']}  Columns: {payload['column_count']}")
    detected = payload["detected_columns"]
    print(f"Cluster column: {detected.get('cluster_column') or 'not detected'}")
    print(f"Gene/marker column: {detected.get('gene_column') or 'not detected'}")
    print(f"Ranking method: {detected.get('ranking_method')}")

    summary = payload.get("summary") or {}
    if summary.get("cluster_count"):
        print("")
        print(f"Clusters: {summary['cluster_count']}")
        print(
            "Markers per cluster: min={min_markers}, median={median_markers}, max={max_markers}".format(
                **summary
            )
        )
        print(f"Total markers prepared: {summary['total_markers']}")

    if payload.get("errors"):
        print("")
        print("Errors:")
        for error in payload["errors"]:
            print(f"  - {error}")
    if payload.get("warnings"):
        print("")
        print("Warnings:")
        for warning in payload["warnings"]:
            print(f"  - {warning}")
    if payload.get("info"):
        print("")
        print("Info:")
        for item in payload["info"]:
            print(f"  - {item}")

    if payload.get("recommended_command"):
        print("")
        print("Suggested next command:")
        print(f"  {payload['recommended_command']}")


def run_validate(args: Any) -> int:
    """Validate marker input from CLI arguments."""
    payload = validate_marker_input(args)
    if args.json:
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        _print_human_report(payload)
    if payload["status"] == "error":
        return 1
    if payload["status"] == "warning" and args.strict:
        return 1
    return 0
