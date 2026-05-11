"""
Build literature-derived macrophage subtype benchmark case files.

The generated CSVs use marker genes and subtype labels from local downloaded
supplementary tables, so benchmark inputs are traceable to papers rather than
hand-crafted examples.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from urllib.request import urlretrieve
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


SUITE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
PAPERS = ROOT / "CASSIA_python" / "CASSIA" / "agents" / "reference_agent" / "macrophage_test" / "papers" / "downloads"
CASES_DIR = SUITE_DIR / "cases"
OUT_DIR = CASES_DIR / "literature"
HELDOUT_DIR = CASES_DIR / "heldout"
SOURCE_DATA_DIR = SUITE_DIR / "source_data"
TOP_MARKERS_PER_CASE = 30

COULTON_XLSX = (
    PAPERS
    / "Coulton_2024_Using_a_pan-cancer_atlas_to"
    / "supplementary"
    / "41467_2024_49885_MOESM4_ESM"
    / "41467_2024_49885_MOESM4_ESM.xlsx"
)
WANG_XLSX = (
    PAPERS
    / "Wang_2023_An_immune_cell_atlas_reveals"
    / "supplementary"
    / "mmc3"
    / "mmc3.xlsx"
)
LI_CSV = (
    PAPERS
    / "Li_2023_Single-cell_characterization_of_macrophages_in"
    / "supplementary"
    / "12276_2023_1115_MOESM3_ESM"
    / "12276_2023_1115_MOESM3_ESM.csv"
)
LI_2024_MARKER_PDF = (
    PAPERS
    / "Li_2024_PanCancer_Myeloid_ICB"
    / "supplementary"
    / "41467_2024_50478_MOESM5_ESM"
    / "41467_2024_50478_MOESM5_ESM.pdf"
)
LI_2024_MARKER_PDF_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-50478-8/"
    "MediaObjects/41467_2024_50478_MOESM5_ESM.pdf"
)
QI_2022_XLSX = (
    SOURCE_DATA_DIR
    / "heldout"
    / "qi_2022_crc"
    / "41467_2022_29366_MOESM5_ESM.xlsx"
)
QI_2022_XLSX_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-29366-6/"
    "MediaObjects/41467_2022_29366_MOESM5_ESM.xlsx"
)

LI_EXPECTED_TERMS = {
    "M-C1": "IL1B/TNF inflammatory TAM/inflammatory macrophage;CCL3/CCL4/TNF/IL1B;macrophage/TAM",
    "M-C2": "inflammatory activated macrophage/immediate-early inflammatory macrophage;CXCL8/CCL3/IL1B/TNFAIP3/NFKBIA;macrophage/TAM",
    "M-C3": "C1QC/C3 complement antigen-presenting TAM/antigen-presenting macrophage;C1QA/C1QB/C1QC/HLA/CD74;macrophage/TAM",
    "M-C4": "metabolic/proliferation/low inflammatory/S100A6/FABP5;uveal melanoma macrophage/macrophage",
}

LI_2024_EXPECTED_TERMS = {
    "Macro_FOLR2-APOE+": "APOE/TREM2;lipid/lipid-phagolysosomal/lipid-associated;macrophage/TAM",
    "Macro_FOLR2+APOE-": "FOLR2/SELENOP/SLC40A1;resident-like/tissue-resident/iron-handling;macrophage/TAM",
    "Macro_FOLR2+APOE+": "GPNMB/CCL18/CXCL9;inflammatory/immunosuppressive/antigen-presenting/lipid-phagolysosomal;macrophage/TAM",
    "Macro_IER3": "IER3/TNF/CCL3/CCL4/JUN;immediate-early/inflammatory/stress-activated;macrophage/TAM",
    "Macro_IFI27": "IFI27/interferon/IFN;APOE/C1QA/C1QB/C1QC/complement;macrophage/TAM",
    "Macro_ISG15": "ISG15/IFIT1/IFIT/type-I interferon;interferon/IFN/ISG;macrophage/TAM",
    "Macro_LYVE1": "FOLR2/SELENOP/SLC40A1/F13A1;resident-like/tissue-resident/iron-handling;macrophage/TAM",
    "Macro_NLRP3": "NLRP3/IL1B/EREG;inflammasome/inflammatory/IL-1;macrophage/TAM",
    "Macro_OLFML3": "CXCL9/CXCL10/CXCL11/GBP;IFNG/IFN-gamma/M1-like/interferon-gamma;macrophage/TAM",
}

LI_2024_MACROPHAGE_LABELS = [
    "Macro_FOLR2-APOE+",
    "Macro_FOLR2+APOE-",
    "Macro_FOLR2+APOE+",
    "Macro_IER3",
    "Macro_IFI27",
    "Macro_ISG15",
    "Macro_LYVE1",
    "Macro_NLRP3",
    "Macro_OLFML3",
]

QI_2022_EXPECTED_TERMS = {
    "THBS1+ Macrophage": (
        "STAB1/CD163/FOLR2/SELENOP/MERTK;"
        "resident-like/immunoregulatory/phagocytic;"
        "macrophage/TAM"
    ),
    "MARCO+ Macrophage": (
        "MARCO/SPP1;"
        "SPP1/TREM2/lipid/inflammatory/ECM;"
        "macrophage/TAM"
    ),
    "VCAN+ Monocyte": (
        "VCAN/FCN1;"
        "IL1B/EREG/inflammatory;"
        "monocyte/macrophage"
    ),
    "Proliferating Myeloid cells": (
        "MKI67/TOP2A;"
        "proliferating/cycling/cell cycle;"
        "myeloid/macrophage"
    ),
}

QI_2022_HELDOUT_LABELS = [
    "THBS1+ Macrophage",
    "MARCO+ Macrophage",
    "VCAN+ Monocyte",
    "Proliferating Myeloid cells",
]


COULTON_EXPECTED_TERMS = {
    "0_AlvMac": "AlvMac/alveolar;macrophage/TAM",
    "1_MetM2Mac": "MetM2Mac/FOLR2/resident;SELENOP/SLC40A1;macrophage/TAM",
    "2_C3Mac": "C3Mac/C3;complement;macrophage/TAM",
    "3_ICIMac1": "ICIMac;SPP1/TREM2;macrophage/TAM",
    "4_ICIMac2": "ICIMac/TREM2/lipid;APOE;macrophage/TAM",
    "5_StressMac": "StressMac/stress;heat/HSP;macrophage/TAM",
    "6_SPP1AREGMac": "SPP1AREGMac/SPP1/AREG/EREG;inflammatory/angiogenic;macrophage/TAM",
    "7_IFNMac": "IFNMac/interferon/IFN;CCL2/CCL8;macrophage/TAM",
    "8_IFNGMac": "IFNGMac/IFNG/IFN-gamma;CXCL9/CXCL10;macrophage/TAM",
    "9_AngioMac": "AngioMac/angiogenic/angiogenesis;AREG/EREG;macrophage/TAM",
    "10_InflamMac": "InflamMac/inflammatory;IL1B/TNF;macrophage/TAM",
    "11_MetalloMac": "MetalloMac/metallothionein;MT1/MT2/metal;macrophage/TAM",
    "12_MBMMac": "MBMMac;KCNMA1/LRMDA;macrophage/TAM",
    "13_CalciumMac": "CalciumMac;S100A6/S100A10;macrophage/TAM",
    "14_ProliMac": "ProliMac/proliferating;MKI67/cell cycle;macrophage/TAM",
    "15_LYZMac": "LYZMac;LYZ;macrophage/TAM",
    "16_ECMHomeoMac": "ECMHomeoMac/ECM/matrix;MMP9/TIMP1/SPP1;macrophage/TAM",
    "17_IFNMac3": "IFNMac/interferon/IFN;ISG15/IFIT;macrophage/TAM",
    "18_ECMMac": "ECMMac/ECM/matrix;COL1A1/COL1A2/collagen;macrophage/TAM",
    "19_ClassMono": "ClassMono/classical monocyte;S100A8/S100A9/FCN1;monocyte",
    "21_HemeMac": "HemeMac/heme/iron;HMOX1/SLC40A1;macrophage/TAM",
    "22_IFNMac4": "IFNMac/interferon/IFN;IFITM2/IFITM3;macrophage/TAM",
}

COULTON_COVERED_LABELS = [
    "1_MetM2Mac",
    "2_C3Mac",
    "4_ICIMac2",
    "5_StressMac",
    "6_SPP1AREGMac",
    "8_IFNGMac",
    "9_AngioMac",
    "10_InflamMac",
    "11_MetalloMac",
    "17_IFNMac3",
    "18_ECMMac",
    "21_HemeMac",
]


def _cluster_num(label: str) -> int:
    return int(str(label).split("_", 1)[0])


def _join_markers(markers: Iterable[str]) -> str:
    return ", ".join(str(marker).strip() for marker in markers if str(marker).strip())


def _parse_scientific_number(value: str) -> float:
    return float(str(value).replace(",", "."))


def _top_positive_markers(
    df: pd.DataFrame,
    cluster_value,
    cluster_col: str,
    gene_col: str,
    fc_col: str,
    n: int = TOP_MARKERS_PER_CASE,
    order_col: Optional[str] = None,
) -> List[str]:
    subset = df[(df[cluster_col] == cluster_value) & (pd.to_numeric(df[fc_col], errors="coerce") > 0)].copy()
    if order_col and order_col in subset.columns:
        subset = subset.sort_values([order_col, fc_col], ascending=[True, False])
    elif "p_val" in subset.columns:
        subset = subset.sort_values(["p_val", fc_col], ascending=[True, False])
    else:
        subset = subset.sort_values(fc_col, ascending=False)
    return list(subset[gene_col].astype(str).head(n))


def build_coulton_cases() -> Dict[str, pd.DataFrame]:
    markers = pd.read_excel(COULTON_XLSX, sheet_name="S. Data 6", header=1)
    composition = pd.read_excel(COULTON_XLSX, sheet_name="S. Data 4", header=1)
    labels = sorted(composition["short.label"].dropna().astype(str).unique(), key=_cluster_num)

    rows = []
    case_index = 1
    for label in labels:
        if label in {"20_TDoub", "23_NA"}:
            continue
        genes = _top_positive_markers(
            markers,
            cluster_value=_cluster_num(label),
            cluster_col="cluster",
            gene_col="gene",
            fc_col="avg_log2FC",
            order_col="order",
        )
        if not genes:
            continue
        rows.append({
            "case_id": f"coulton_2024_case_{case_index:02d}",
            "source": "Coulton 2024 Nat Commun",
            "doi": "10.1038/s41467-024-49885-8",
            "paper_subtype": label,
            "markers": _join_markers(genes),
            "expected_terms": COULTON_EXPECTED_TERMS.get(label, f"{label};macrophage/TAM"),
            "min_score": 2,
            "source_file": str(COULTON_XLSX.relative_to(ROOT)),
            "source_sheet": "S. Data 6",
            "evidence_note": f"Top {TOP_MARKERS_PER_CASE} positive differentially expressed markers from Supplementary Data 6.",
        })
        case_index += 1

    all_cases = pd.DataFrame(rows)
    covered_cases = all_cases[all_cases["paper_subtype"].isin(COULTON_COVERED_LABELS)].reset_index(drop=True)
    return {
        "coulton_2024_pan_cancer_tam_all.csv": all_cases.reset_index(drop=True),
        "coulton_2024_pan_cancer_tam_covered.csv": covered_cases,
    }


def build_wang_cases() -> pd.DataFrame:
    df = pd.read_excel(WANG_XLSX, sheet_name="Macrophage", header=1)
    rows = []
    for case_index, subtype in enumerate(sorted(df["Subtype"].dropna().astype(str).unique()), start=1):
        genes = _top_positive_markers(
            df,
            cluster_value=subtype,
            cluster_col="Subtype",
            gene_col="Gene",
            fc_col="avg_log2FC",
            n=TOP_MARKERS_PER_CASE,
        )
        rows.append({
            "case_id": f"wang_2023_case_{case_index:02d}",
            "source": "Wang 2023 Cell",
            "doi": "10.1016/j.cell.2023.08.019",
            "paper_subtype": subtype,
            "markers": _join_markers(genes),
            "expected_terms": f"{subtype};macrophage",
            "min_score": 1,
            "source_file": str(WANG_XLSX.relative_to(ROOT)),
            "source_sheet": "Macrophage",
            "evidence_note": f"Top {TOP_MARKERS_PER_CASE} positive macrophage subtype DEGs from Supplementary Table mmc3.",
        })
    return pd.DataFrame(rows)


def build_li_cases() -> pd.DataFrame:
    df = pd.read_csv(LI_CSV)
    rows = []
    for case_index, cluster in enumerate(sorted(df["cluster"].dropna().astype(str).unique()), start=1):
        subset = df[
            (df["cluster"].astype(str) == cluster)
            & (pd.to_numeric(df["avg_log2FC"], errors="coerce") > 0)
        ].copy()
        genes = list(
            subset.sort_values("avg_log2FC", ascending=False)["gene"]
            .astype(str)
            .head(TOP_MARKERS_PER_CASE)
        )
        rows.append({
            "case_id": f"li_2023_uveal_melanoma_case_{case_index:02d}",
            "source": "Li 2023 Exp Mol Med",
            "doi": "10.1038/s12276-023-01115-9",
            "paper_subtype": cluster,
            "markers": _join_markers(genes),
            "expected_terms": LI_EXPECTED_TERMS.get(cluster, f"{cluster};uveal melanoma;macrophage"),
            "min_score": 2 if cluster != "M-C4" else 1,
            "source_file": str(LI_CSV.relative_to(ROOT)),
            "source_sheet": "MOESM3",
            "evidence_note": (
                f"Top {TOP_MARKERS_PER_CASE} positive macrophage cluster DEGs from Supplementary Data MOESM3, "
                "ranked by avg_log2FC. Expected terms map paper clusters to the current "
                "human cancer macrophage consensus where possible; M-C4 is UM-specific "
                "and weakly covered by the current consensus layer."
            ),
        })
    return pd.DataFrame(rows)


def _li_2024_marker_rows() -> pd.DataFrame:
    if not LI_2024_MARKER_PDF.exists():
        LI_2024_MARKER_PDF.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(LI_2024_MARKER_PDF_URL, LI_2024_MARKER_PDF)

    pdftotext = shutil.which("pdftotext")
    if pdftotext is None:
        raise RuntimeError(
            "pdftotext is required to parse Li 2024 Supplementary Data 2 PDF. "
            "Install poppler or provide a pre-extracted marker table."
        )

    completed = subprocess.run(
        [pdftotext, "-layout", str(LI_2024_MARKER_PDF), "-"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    row_pattern = re.compile(
        r"^\s*(?:\d+)?"
        r"(?P<cluster>[A-Za-z][A-Za-z0-9_+\-()]+)\s+"
        r"(?P<gene>[A-Za-z0-9.\-]+)\s+"
        r"(?P<p_val>[\d,]+E[+\-]\d+)\s+"
        r"(?P<avg_log2FC>[+\-]?[\d,]+E[+\-]\d+)\s+"
        r"(?P<pct_1>[\d,]+)\s+"
        r"(?P<pct_2>[\d,]+)\s+"
        r"(?P<p_val_adj>[\d,]+E[+\-]\d+)"
    )
    rows = []
    for line in completed.stdout.splitlines():
        match = row_pattern.match(line)
        if not match:
            continue
        parsed = match.groupdict()
        rows.append({
            "cluster": parsed["cluster"],
            "gene": parsed["gene"],
            "p_val": _parse_scientific_number(parsed["p_val"]),
            "avg_log2FC": _parse_scientific_number(parsed["avg_log2FC"]),
            "pct_1": _parse_scientific_number(parsed["pct_1"]),
            "pct_2": _parse_scientific_number(parsed["pct_2"]),
            "p_val_adj": _parse_scientific_number(parsed["p_val_adj"]),
        })
    return pd.DataFrame(rows)


def build_li_2024_cases() -> pd.DataFrame:
    df = _li_2024_marker_rows()
    rows = []
    for case_index, cluster in enumerate(LI_2024_MACROPHAGE_LABELS, start=1):
        subset = df[(df["cluster"] == cluster) & (df["avg_log2FC"] > 0)]
        genes = list(subset["gene"].astype(str).head(TOP_MARKERS_PER_CASE))
        rows.append({
            "case_id": f"li_2024_pan_cancer_icb_case_{case_index:02d}",
            "source": "Li 2024 Nat Commun",
            "doi": "10.1038/s41467-024-50478-8",
            "paper_subtype": cluster,
            "markers": _join_markers(genes),
            "expected_terms": LI_2024_EXPECTED_TERMS[cluster],
            "min_score": 2,
            "source_file": str(LI_2024_MARKER_PDF.relative_to(ROOT)),
            "source_sheet": "Supplementary Data 2 PDF",
            "evidence_note": (
                f"Top {TOP_MARKERS_PER_CASE} positive macrophage subtype marker genes parsed from Supplementary Data 2. "
                "The paper reports a pan-cancer ICB myeloid atlas across eight cancer types and "
                "identifies macrophage/monocyte subtypes from 47,750 myeloid cells."
            ),
        })
    return pd.DataFrame(rows)


def build_qi_2022_crc_heldout_cases() -> pd.DataFrame:
    if not QI_2022_XLSX.exists():
        QI_2022_XLSX.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(QI_2022_XLSX_URL, QI_2022_XLSX)

    df = pd.read_excel(QI_2022_XLSX, sheet_name="Myeloid Cells")
    rows = []
    for case_index, cluster in enumerate(QI_2022_HELDOUT_LABELS, start=1):
        genes = _top_positive_markers(
            df,
            cluster_value=cluster,
            cluster_col="cluster",
            gene_col="gene",
            fc_col="avg_logFC",
            n=TOP_MARKERS_PER_CASE,
        )
        rows.append({
            "case_id": f"qi_2022_crc_heldout_case_{case_index:02d}",
            "source": "Qi 2022 Nat Commun",
            "doi": "10.1038/s41467-022-29366-6",
            "paper_subtype": cluster,
            "markers": _join_markers(genes),
            "expected_terms": QI_2022_EXPECTED_TERMS[cluster],
            "min_score": 2,
            "source_file": str(QI_2022_XLSX.relative_to(ROOT)),
            "source_sheet": "Myeloid Cells",
            "evidence_note": (
                f"Held-out colorectal cancer myeloid subtype case. Top {TOP_MARKERS_PER_CASE} "
                "positive marker genes were taken from Supplementary Data 2, Myeloid Cells. "
                "This paper is intentionally not included in the macrophage reference brain."
            ),
        })
    return pd.DataFrame(rows)


def build_inventory(case_files: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    return pd.DataFrame([
        {
            "paper": "Coulton 2024, pan-cancer TAM atlas",
            "doi": "10.1038/s41467-024-49885-8",
            "local_status": "PDF and Supplementary Data 1-7 available locally",
            "usable_for_benchmark": "yes - primary",
            "reason": "Author-defined TAM clusters with differential marker table; best match to current tumor macrophage reference docs.",
            "case_file": "cases/literature/coulton_2024_pan_cancer_tam_covered.csv",
            "num_cases": len(case_files["coulton_2024_pan_cancer_tam_covered.csv"]),
        },
        {
            "paper": "Coulton 2024, pan-cancer TAM atlas",
            "doi": "10.1038/s41467-024-49885-8",
            "local_status": "PDF and Supplementary Data 1-7 available locally",
            "usable_for_benchmark": "yes - extended",
            "reason": "Broader cluster set includes monocyte-like and less-covered TAM labels; useful after reference brain is expanded.",
            "case_file": "cases/literature/coulton_2024_pan_cancer_tam_all.csv",
            "num_cases": len(case_files["coulton_2024_pan_cancer_tam_all.csv"]),
        },
        {
            "paper": "Wang 2023, prenatal macrophage atlas",
            "doi": "10.1016/j.cell.2023.08.019",
            "local_status": "Supplementary XLSX files available locally",
            "usable_for_benchmark": "yes - future domain",
            "reason": "Fifteen developmental macrophage subtypes with DEGs; not a tumor macrophage benchmark, so current reference coverage is partial.",
            "case_file": "cases/literature/wang_2023_prenatal_macrophage.csv",
            "num_cases": len(case_files["wang_2023_prenatal_macrophage.csv"]),
        },
        {
            "paper": "Li 2023, uveal melanoma macrophages",
            "doi": "10.1038/s12276-023-01115-9",
            "local_status": "PDF and two CSV supplements available locally",
            "usable_for_benchmark": "partial",
            "reason": "Four macrophage clusters have DEGs, but labels are M-C1 to M-C4; needs paper-derived semantic mapping before strict subtype scoring.",
            "case_file": "cases/literature/li_2023_uveal_melanoma_macrophage.csv",
            "num_cases": len(case_files["li_2023_uveal_melanoma_macrophage.csv"]),
        },
        {
            "paper": "Li 2024, pan-cancer ICB myeloid atlas",
            "doi": "10.1038/s41467-024-50478-8",
            "local_status": "Supplementary Data 2 PDF marker table downloaded locally",
            "usable_for_benchmark": "yes - independent pan-cancer macrophage",
            "reason": "Nine author-defined Macro_* clusters have positive marker genes and map well to the human cancer macrophage consensus layer.",
            "case_file": "cases/literature/li_2024_pan_cancer_icb_macrophage.csv",
            "num_cases": len(case_files["li_2024_pan_cancer_icb_macrophage.csv"]),
        },
        {
            "paper": "Qi 2022, colorectal cancer myeloid atlas",
            "doi": "10.1038/s41467-022-29366-6",
            "local_status": "Supplementary Data 2 XLSX downloaded on demand to Benchmark/reference_subtype/source_data",
            "usable_for_benchmark": "yes - heldout",
            "reason": "Author-defined myeloid/macrophage subtypes with DEG table; intentionally not used as reference-brain evidence.",
            "case_file": "cases/heldout/qi_2022_crc_macrophage.csv",
            "num_cases": len(case_files["qi_2022_crc_macrophage.csv"]),
        },
        {
            "paper": "Cheng 2021, pan-cancer tumor-infiltrating myeloid atlas",
            "doi": "10.1016/j.cell.2021.01.010",
            "local_status": "PDF and Supplementary Tables mmc1-mmc6 available locally",
            "usable_for_benchmark": "reference evidence",
            "reason": "High-impact source for angiogenic TAM facts and signatures, but local supplements do not provide a clean author-labeled macrophage subtype marker benchmark table.",
            "case_file": "",
            "num_cases": 0,
        },
        {
            "paper": "Mulder 2021, MNP-VERSE",
            "doi": "10.1016/j.immuni.2021.07.007",
            "local_status": "PDF and Supplementary Tables mmc2-mmc8 available locally",
            "usable_for_benchmark": "partial",
            "reason": "Strong cross-tissue MoMac signatures; cluster IDs need semantic label mapping before strict subtype scoring.",
            "case_file": "",
            "num_cases": 0,
        },
        {
            "paper": "Ochocka 2021, glioma-associated brain macrophages",
            "doi": "10.1038/s41467-021-21407-w",
            "local_status": "PDF and Supplementary Data 3 available locally",
            "usable_for_benchmark": "partial",
            "reason": "Top DEGs by MG cluster are available; needs microglia/macrophage label mapping and a CNS-specific reference branch.",
            "case_file": "",
            "num_cases": 0,
        },
    ])


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    HELDOUT_DIR.mkdir(parents=True, exist_ok=True)
    case_files = build_coulton_cases()
    case_files["wang_2023_prenatal_macrophage.csv"] = build_wang_cases()
    case_files["li_2023_uveal_melanoma_macrophage.csv"] = build_li_cases()
    case_files["li_2024_pan_cancer_icb_macrophage.csv"] = build_li_2024_cases()
    case_files["qi_2022_crc_macrophage.csv"] = build_qi_2022_crc_heldout_cases()

    for filename, cases in case_files.items():
        outdir = HELDOUT_DIR if filename == "qi_2022_crc_macrophage.csv" else OUT_DIR
        cases.to_csv(outdir / filename, index=False)

    inventory = build_inventory(case_files)
    inventory.to_csv(OUT_DIR / "literature_case_inventory.csv", index=False)

    print(f"Wrote literature benchmark cases to {OUT_DIR}")
    for filename, cases in case_files.items():
        print(f"  {filename}: {len(cases)} cases")
    print(f"  literature_case_inventory.csv: {len(inventory)} papers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
