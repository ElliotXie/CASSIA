"""
Build the core runCASSIA_batch benchmark case manifest.

The marker source is CASSIA_example/Benchmark/marker100.xlsx. Historical
baseline correctness is optionally joined from Supplementary_Data 4.xlsx when
that file is available locally.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import openpyxl
import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = ROOT / "CASSIA_example" / "Benchmark" / "marker100.xlsx"
DEFAULT_BASELINE_SOURCE = Path(
    "/Volumes/Untitled/bioprojectbackup/cassia/ultimatesubmission/"
    "Final_submission/Supplement/Supplementary_Data 4.xlsx"
)
DEFAULT_OUTPUT = SUITE_DIR / "cases" / "core_50.csv"
DEFAULT_INVENTORY = SUITE_DIR / "cases" / "source_inventory.csv"


CASE_FIELDS = [
    "case_id",
    "dataset",
    "tissue",
    "species",
    "expected_cell_type",
    "broad_cell_type",
    "difficulty",
    "case_type",
    "expected_terms",
    "min_score",
    "selection_reason",
    "source_file",
    "source_row",
    "source_cassia_correctness",
    "source_gptcelltype_4o_correctness",
    "source_gptcelltype_4_correctness",
    "source_singleR_correctness",
    "source_celltypist_correctness",
    "marker_list",
]


# 50 selected rows from marker100.xlsx:
#   GTEx lung: all 15 high-quality cases.
#   Azimuth kidney: 22 representative canonical, subtype, and hard renal cases.
#   TS large intestine: 8 non-eye cases, including rare epithelial states.
#   HCL fetal skin: 5 diverse fetal skin cases with historical baseline rows.
SELECTED_CASES: List[Tuple[str, str, str]] = [
    ("GTEx", "lung", "Endothelial cell (lymphatic)"),
    ("GTEx", "lung", "Endothelial cell (vascular)"),
    ("GTEx", "lung", "Epithelial cell (alveolar type I)"),
    ("GTEx", "lung", "Epithelial cell (alveolar type II)"),
    ("GTEx", "lung", "Epithelial cell (basal)"),
    ("GTEx", "lung", "Epithelial cell (ciliated)"),
    ("GTEx", "lung", "Epithelial cell (club)"),
    ("GTEx", "lung", "Fibroblast"),
    ("GTEx", "lung", "Immune (alveolar macrophage)"),
    ("GTEx", "lung", "Immune (B cell)"),
    ("GTEx", "lung", "Immune (DC/macrophage)"),
    ("GTEx", "lung", "Immune (mast cell)"),
    ("GTEx", "lung", "Immune (NK cell)"),
    ("GTEx", "lung", "Immune (T cell)"),
    ("GTEx", "lung", "Pericyte/SMC"),
    ("Azimuth", "kidney", "B cell"),
    ("Azimuth", "kidney", "CD4 T cell"),
    ("Azimuth", "kidney", "CD8 T cell"),
    ("Azimuth", "kidney", "NK cell"),
    ("Azimuth", "kidney", "NKT cell"),
    ("Azimuth", "kidney", "Mast cell"),
    ("Azimuth", "kidney", "Neutrophil"),
    ("Azimuth", "kidney", "MNP-a/classical monocyte derived"),
    ("Azimuth", "kidney", "MNP-b/non-classical monocyte derived"),
    ("Azimuth", "kidney", "MNP-c/dendritic cell"),
    ("Azimuth", "kidney", "MNP-d/Tissue macrophage"),
    ("Azimuth", "kidney", "Plasmacytoid dendritic cell"),
    ("Azimuth", "kidney", "Fibroblast"),
    ("Azimuth", "kidney", "Myofibroblast"),
    ("Azimuth", "kidney", "Podocyte"),
    ("Azimuth", "kidney", "Proximal tubule"),
    ("Azimuth", "kidney", "Thick ascending limb of Loop of Henle"),
    ("Azimuth", "kidney", "Intercalated cell"),
    ("Azimuth", "kidney", "Principal cell"),
    ("Azimuth", "kidney", "Connecting tubule"),
    ("Azimuth", "kidney", "Glomerular endothelium"),
    ("Azimuth", "kidney", "Peritubular capillary endothelium"),
    ("TS", "large intestine", "gut endothelial cell"),
    ("TS", "large intestine", "large intestine goblet cell"),
    ("TS", "large intestine", "intestinal enteroendocrine cell"),
    ("TS", "large intestine", "intestinal tuft cell"),
    ("TS", "large intestine", "transit amplifying cell of large intestine"),
    ("TS", "large intestine", "intestinal crypt stem cell"),
    ("TS", "large intestine", "plasma cell"),
    ("TS", "large intestine", "paneth cell of epithelium of large intestine"),
    ("HCL", "Fetal skin", "Dermis fibroblast"),
    ("HCL", "Fetal skin", "Endothelial cell"),
    ("HCL", "Fetal skin", "Fibroblast_EFEMP1 high"),
    ("HCL", "Fetal skin", "Macrophage"),
    ("HCL", "Fetal skin", "Melanocyte_MLANA high"),
]


HARD_CASES = {
    ("GTEx", "lung", "Immune (DC/macrophage)"),
    ("GTEx", "lung", "Pericyte/SMC"),
    ("Azimuth", "kidney", "CD8 T cell"),
    ("Azimuth", "kidney", "NKT cell"),
    ("Azimuth", "kidney", "Myofibroblast"),
    ("Azimuth", "kidney", "Connecting tubule"),
    ("Azimuth", "kidney", "Glomerular endothelium"),
    ("Azimuth", "kidney", "Peritubular capillary endothelium"),
    ("TS", "large intestine", "intestinal tuft cell"),
    ("TS", "large intestine", "transit amplifying cell of large intestine"),
    ("TS", "large intestine", "intestinal crypt stem cell"),
    ("TS", "large intestine", "paneth cell of epithelium of large intestine"),
}


EASY_CASES = {
    ("GTEx", "lung", "Epithelial cell (basal)"),
    ("GTEx", "lung", "Epithelial cell (ciliated)"),
    ("GTEx", "lung", "Fibroblast"),
    ("GTEx", "lung", "Immune (alveolar macrophage)"),
    ("GTEx", "lung", "Immune (B cell)"),
    ("GTEx", "lung", "Immune (mast cell)"),
    ("GTEx", "lung", "Immune (NK cell)"),
    ("GTEx", "lung", "Immune (T cell)"),
    ("Azimuth", "kidney", "B cell"),
    ("Azimuth", "kidney", "CD4 T cell"),
    ("Azimuth", "kidney", "NK cell"),
    ("Azimuth", "kidney", "Mast cell"),
    ("Azimuth", "kidney", "Neutrophil"),
    ("Azimuth", "kidney", "Fibroblast"),
    ("Azimuth", "kidney", "Podocyte"),
    ("Azimuth", "kidney", "Proximal tubule"),
    ("TS", "large intestine", "plasma cell"),
    ("HCL", "Fetal skin", "Macrophage"),
    ("HCL", "Fetal skin", "Melanocyte_MLANA high"),
}


BASELINE_COLUMNS = [
    "CASSIA_correctness",
    "gptcelltype_4o_correctness",
    "gptcelltype_4_correctness",
    "singleR_correctness",
    "celltypist_correctness",
]


def normalize(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def slugify(value: object) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", str(value or "").lower()).strip("_")
    return slug or "case"


def source_key(dataset: object, tissue: object, cell_type: object) -> Tuple[str, str, str]:
    return normalize(dataset), normalize(tissue), normalize(cell_type)


def case_id(dataset: str, tissue: str, expected_cell_type: str) -> str:
    return f"{slugify(dataset)}_{slugify(tissue)}_{slugify(expected_cell_type)}"


def broad_cell_type(label: str) -> str:
    lower = label.lower()
    if "alveolar type i" in lower:
        return "Alveolar type I epithelial cell"
    if "alveolar type ii" in lower:
        return "Alveolar type II epithelial cell"
    if "basal" in lower:
        return "Basal epithelial cell"
    if "ciliated" in lower:
        return "Ciliated epithelial cell"
    if "club" in lower:
        return "Club epithelial cell"
    if "epithelial" in lower or "enterocyte" in lower or "goblet" in lower or "paneth" in lower:
        if "goblet" in lower:
            return "Goblet cell"
        if "enteroendocrine" in lower:
            return "Enteroendocrine cell"
        if "tuft" in lower:
            return "Tuft cell"
        if "crypt stem" in lower:
            return "Intestinal crypt stem cell"
        if "transit" in lower:
            return "Transit-amplifying epithelial cell"
        if "paneth" in lower:
            return "Paneth cell"
        return "Epithelial cell"
    if "lymphatic" in lower:
        return "Lymphatic endothelial cell"
    if "endothelial" in lower or "endothelium" in lower:
        return "Endothelial cell"
    if "fibroblast" in lower:
        return "Fibroblast"
    if "pericyte" in lower or "smc" in lower or "myofibroblast" in lower:
        return "Perivascular stromal cell"
    if "macrophage" in lower or "mnp" in lower or "monocyte" in lower:
        return "Mononuclear phagocyte"
    if "dendritic" in lower:
        return "Dendritic cell"
    if "plasmacytoid" in lower:
        return "Plasmacytoid dendritic cell"
    if "mast" in lower:
        return "Mast cell"
    if "neutrophil" in lower:
        return "Neutrophil"
    if "nk" in lower:
        return "NK/NKT cell" if "nkt" in lower else "NK cell"
    if "cd4" in lower:
        return "CD4 T cell"
    if "cd8" in lower:
        return "CD8 T cell"
    if "t cell" in lower or "t-cell" in lower:
        return "T cell"
    if "b cell" in lower:
        return "B cell"
    if "podocyte" in lower:
        return "Podocyte"
    if "proximal tubule" in lower:
        return "Proximal tubule cell"
    if "thick ascending" in lower:
        return "Thick ascending limb cell"
    if "intercalated" in lower:
        return "Intercalated cell"
    if "principal" in lower:
        return "Principal cell"
    if "connecting tubule" in lower:
        return "Connecting tubule cell"
    if "melanocyte" in lower:
        return "Melanocyte"
    return label


def expected_terms(label: str) -> str:
    lower = label.lower()
    if "lymphatic" in lower:
        return "lymphatic;endothelial"
    if "vascular" in lower:
        return "vascular;endothelial"
    if "glomerular" in lower:
        return "glomerular;endothelial"
    if "peritubular" in lower:
        return "peritubular/capillary;endothelial"
    if "endothelial" in lower or "endothelium" in lower:
        return "endothelial"
    if "alveolar type i" in lower:
        return "alveolar type i/AT1/type I;epithelial"
    if "alveolar type ii" in lower:
        return "alveolar type ii/AT2/type II;epithelial"
    if "basal" in lower:
        return "basal;epithelial"
    if "ciliated" in lower:
        return "ciliated;epithelial"
    if "club" in lower:
        return "club/Clara;epithelial"
    if "dc/macrophage" in lower:
        return "dendritic/macrophage/myeloid;antigen-presenting/APC"
    if "alveolar macrophage" in lower:
        return "alveolar;macrophage"
    if "macrophage" in lower:
        return "macrophage/myeloid"
    if "mnp-a" in lower:
        return "monocyte/macrophage/myeloid;classical"
    if "mnp-b" in lower:
        return "monocyte/macrophage/myeloid;non-classical/nonclassical/FCGR3A"
    if "mnp-c" in lower:
        return "dendritic;myeloid/antigen-presenting/APC"
    if "mnp-d" in lower:
        return "macrophage;tissue-resident/resident"
    if "plasmacytoid" in lower:
        return "plasmacytoid/pDC;dendritic"
    if "mast" in lower:
        return "mast"
    if "neutrophil" in lower:
        return "neutrophil"
    if "nkt" in lower:
        return "NKT/NK T/natural killer T;T cell/NK"
    if "nk" in lower:
        return "NK/natural killer"
    if "cd4" in lower:
        return "CD4;T cell/T-cells/T cells"
    if "cd8" in lower:
        return "CD8/cytotoxic;T cell/T-cells/T cells"
    if "b cell" in lower:
        return "B cell/B cells"
    if "fibroblast" in lower:
        return "fibroblast"
    if "myofibroblast" in lower:
        return "myofibroblast;smooth muscle/pericyte/fibroblast"
    if "pericyte" in lower or "smc" in lower:
        return "pericyte/smooth muscle/SMC"
    if "podocyte" in lower:
        return "podocyte"
    if "proximal tubule" in lower:
        return "proximal tubule"
    if "thick ascending" in lower:
        return "thick ascending/TAL;loop of henle"
    if "intercalated" in lower:
        return "intercalated;collecting duct"
    if "principal" in lower:
        return "principal;collecting duct"
    if "connecting tubule" in lower:
        return "connecting tubule/CNT;renal tubular/tubule"
    if "goblet" in lower:
        return "goblet;intestinal/large intestine/colon"
    if "enteroendocrine" in lower:
        return "enteroendocrine;intestinal"
    if "tuft" in lower:
        return "tuft;intestinal"
    if "transit amplifying" in lower:
        return "transit amplifying/transit-amplifying;progenitor/proliferating"
    if "crypt stem" in lower:
        return "intestinal stem/crypt stem;stem/LGR5"
    if "paneth" in lower:
        return "paneth"
    if "plasma" in lower:
        return "plasma cell/B cell"
    if "melanocyte" in lower:
        return "melanocyte;MLANA"
    return broad_cell_type(label)


def case_type(key: Tuple[str, str, str]) -> str:
    label = key[2].lower()
    if any(token in label for token in ["tuft", "crypt stem", "transit", "paneth"]):
        return "rare_or_state"
    if key in HARD_CASES:
        return "mixed_or_ambiguous"
    if any(token in label for token in ["subtype", "type ", "mnp", "plasmacytoid"]):
        return "subtype_resolution"
    return "canonical_lineage"


def difficulty(key: Tuple[str, str, str]) -> str:
    if key in HARD_CASES:
        return "hard"
    if key in EASY_CASES:
        return "easy"
    return "medium"


def min_score_for_terms(terms: str, diff: str) -> int:
    return 1


def read_marker_source(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Marker source not found: {path}")
    df = pd.read_excel(path)
    required = {"True Cell Type", "Marker List", "Tissue", "Species", "Dataset"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Marker source is missing columns: {', '.join(missing)}")
    df = df.copy()
    df.insert(0, "source_row", range(2, len(df) + 2))
    return df


def read_baseline_source(path: Optional[Path]) -> Dict[Tuple[str, str, str], Dict[str, object]]:
    if not path or not path.exists():
        return {}

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    if "Annotation" not in workbook.sheetnames:
        return {}
    worksheet = workbook["Annotation"]
    headers = [cell.value for cell in next(worksheet.iter_rows(min_row=1, max_row=1))]
    index = {header: i for i, header in enumerate(headers)}
    lookup: Dict[Tuple[str, str, str], Dict[str, object]] = {}
    for values in worksheet.iter_rows(min_row=2, values_only=True):
        dataset = values[index["dataset"]]
        tissue = values[index["Tissue"]]
        cell_type = values[index["True.Cell.Type"]]
        row = {
            column: values[index[column]]
            for column in BASELINE_COLUMNS
            if column in index
        }
        lookup.setdefault(source_key(dataset, tissue, cell_type), row)
    return lookup


def build_case_rows(source_path: Path, baseline_path: Optional[Path]) -> List[Dict[str, object]]:
    source = read_marker_source(source_path)
    lookup: Dict[Tuple[str, str, str], pd.Series] = {}
    for _idx, row in source.iterrows():
        lookup.setdefault(source_key(row["Dataset"], row["Tissue"], row["True Cell Type"]), row)
    baseline_lookup = read_baseline_source(baseline_path)

    rows: List[Dict[str, object]] = []
    seen_ids = set()
    missing = []
    for dataset, tissue, label in SELECTED_CASES:
        key = source_key(dataset, tissue, label)
        row = lookup.get(key)
        if row is None:
            missing.append(f"{dataset} / {tissue} / {label}")
            continue

        original_key = (dataset, tissue, label)
        terms = expected_terms(label)
        diff = difficulty(original_key)
        cid = case_id(dataset, tissue, label)
        if cid in seen_ids:
            raise RuntimeError(f"Duplicate case_id: {cid}")
        seen_ids.add(cid)
        baseline = baseline_lookup.get(key, {})
        rows.append({
            "case_id": cid,
            "dataset": dataset,
            "tissue": str(row["Tissue"]).lower(),
            "species": str(row["Species"]).lower(),
            "expected_cell_type": label,
            "broad_cell_type": broad_cell_type(label),
            "difficulty": diff,
            "case_type": case_type(original_key),
            "expected_terms": terms,
            "min_score": min_score_for_terms(terms, diff),
            "selection_reason": selection_reason(original_key),
            "source_file": str(source_path),
            "source_row": int(row["source_row"]),
            "source_cassia_correctness": baseline.get("CASSIA_correctness"),
            "source_gptcelltype_4o_correctness": baseline.get("gptcelltype_4o_correctness"),
            "source_gptcelltype_4_correctness": baseline.get("gptcelltype_4_correctness"),
            "source_singleR_correctness": baseline.get("singleR_correctness"),
            "source_celltypist_correctness": baseline.get("celltypist_correctness"),
            "marker_list": row["Marker List"],
        })

    if missing:
        raise RuntimeError("Selected labels missing from marker source:\n- " + "\n- ".join(missing))
    if len(rows) != 50:
        raise RuntimeError(f"Expected 50 selected cases, built {len(rows)}")
    return rows


def selection_reason(key: Tuple[str, str, str]) -> str:
    dataset, _tissue, label = key
    if dataset == "GTEx":
        return "Selected as high-quality GTEx lung benchmark coverage."
    if dataset == "Azimuth":
        if key in HARD_CASES:
            return "Selected as a hard Azimuth kidney subtype/ambiguity case."
        return "Selected as high-quality Azimuth kidney coverage."
    if dataset == "TS":
        return "Selected as non-eye Tabula Sapiens large-intestine coverage."
    if dataset == "HCL":
        return "Selected as diverse HCL fetal-skin coverage with marker-list input."
    return f"Selected from marker100.xlsx: {label}"


def build_inventory(source_path: Path, selected_rows: List[Dict[str, object]], baseline_path: Optional[Path]) -> pd.DataFrame:
    source = read_marker_source(source_path)
    baseline_lookup = read_baseline_source(baseline_path)
    selected_ids = {row["case_id"] for row in selected_rows}

    rows = []
    for _idx, row in source.iterrows():
        dataset = str(row["Dataset"])
        tissue = str(row["Tissue"])
        label = str(row["True Cell Type"])
        cid = case_id(dataset, tissue, label)
        baseline = baseline_lookup.get(source_key(dataset, tissue, label), {})
        rows.append({
            "selected_for_core_50": cid in selected_ids,
            "case_id": cid,
            "dataset": dataset,
            "tissue": tissue,
            "species": row["Species"],
            "expected_cell_type": label,
            "broad_cell_type": broad_cell_type(label),
            "source_row": int(row["source_row"]),
            "source_cassia_correctness": baseline.get("CASSIA_correctness"),
            "source_gptcelltype_4o_correctness": baseline.get("gptcelltype_4o_correctness"),
            "source_gptcelltype_4_correctness": baseline.get("gptcelltype_4_correctness"),
            "marker_list": row["Marker List"],
        })
    return pd.DataFrame(rows)


def write_cases(rows: List[Dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CASE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_summary(rows: List[Dict[str, object]]) -> None:
    df = pd.DataFrame(rows)
    print("By dataset:")
    print(df["dataset"].value_counts().sort_index().to_string())
    print("By difficulty:")
    print(df["difficulty"].value_counts().sort_index().to_string())
    print("Historical baselines on selected rows:")
    for column in [
        "source_cassia_correctness",
        "source_gptcelltype_4o_correctness",
        "source_gptcelltype_4_correctness",
    ]:
        values = pd.to_numeric(df[column], errors="coerce").dropna()
        if values.empty:
            print(f"- {column}: no matched rows")
        else:
            print(
                f"- {column}: n={len(values)}, weighted={values.sum():g}/{len(values)} "
                f"({values.mean():.3f}), full_correct={(values == 1).sum()}/{len(values)}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Build core runCASSIA_batch benchmark cases")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--baseline-source", type=Path, default=DEFAULT_BASELINE_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--inventory-output", type=Path, default=DEFAULT_INVENTORY)
    args = parser.parse_args()

    baseline_source = args.baseline_source if args.baseline_source.exists() else None
    rows = build_case_rows(args.source, baseline_source)
    write_cases(rows, args.output)
    inventory = build_inventory(args.source, rows, baseline_source)
    args.inventory_output.parent.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(args.inventory_output, index=False)

    print(f"Wrote {len(rows)} cases to {args.output}")
    print(f"Wrote source inventory to {args.inventory_output}")
    print_summary(rows)
    if baseline_source is None:
        print("Historical baseline source was not found; baseline columns were left blank.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
