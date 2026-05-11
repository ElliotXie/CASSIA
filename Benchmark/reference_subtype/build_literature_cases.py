"""
Build literature-derived macrophage subtype benchmark case files.

The generated CSVs use marker genes and subtype labels from local downloaded
supplementary tables, so benchmark inputs are traceable to papers rather than
hand-crafted examples.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PAPERS = ROOT / "CASSIA_python" / "CASSIA" / "agents" / "reference_agent" / "macrophage_test" / "papers" / "downloads"
OUT_DIR = Path(__file__).resolve().parent / "literature_cases"

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


COULTON_EXPECTED_TERMS = {
    "0_AlvMac": "AlvMac/alveolar;macrophage/TAM",
    "1_MetM2Mac": "MetM2Mac/FOLR2/resident;SELENOP/SLC40A1;macrophage/TAM",
    "2_C3Mac": "C3Mac/C3;complement;macrophage/TAM",
    "3_ICIMac1": "ICIMac;SPP1/TREM2;macrophage/TAM",
    "4_ICIMac2": "ICIMac/TREM2/lipid;APOE;macrophage/TAM",
    "5_StressMac": "StressMac/stress;heat/HSP;macrophage/TAM",
    "6_SPP1AREGMac": "SPP1AREGMac/SPP1/AREG/EREG;inflammatory/angiogenic;macrophage/TAM",
    "7_IFNMac": "IFNMac/interferon/IFN;CCL2/CCL8;macrophage/TAM",
    "8_IFNGMac": "IFNGMac/IFNG/IFN-gamma/IFN-γ;CXCL9/CXCL10;macrophage/TAM",
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


def _top_positive_markers(
    df: pd.DataFrame,
    cluster_value,
    cluster_col: str,
    gene_col: str,
    fc_col: str,
    n: int = 10,
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
            "evidence_note": "Top positive differentially expressed markers from Supplementary Data 6.",
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
            n=10,
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
            "evidence_note": "Top positive macrophage subtype DEGs from Supplementary Table mmc3.",
        })
    return pd.DataFrame(rows)


def build_li_cases() -> pd.DataFrame:
    df = pd.read_csv(LI_CSV)
    rows = []
    for case_index, cluster in enumerate(sorted(df["cluster"].dropna().astype(str).unique()), start=1):
        genes = _top_positive_markers(
            df,
            cluster_value=cluster,
            cluster_col="cluster",
            gene_col="gene",
            fc_col="avg_log2FC",
            n=10,
        )
        rows.append({
            "case_id": f"li_2023_uveal_melanoma_case_{case_index:02d}",
            "source": "Li 2023 Exp Mol Med",
            "doi": "10.1038/s12276-023-01115-9",
            "paper_subtype": cluster,
            "markers": _join_markers(genes),
            "expected_terms": f"{cluster};uveal melanoma;macrophage",
            "min_score": 1,
            "source_file": str(LI_CSV.relative_to(ROOT)),
            "source_sheet": "MOESM3",
            "evidence_note": "Top positive macrophage cluster DEGs from Supplementary Data MOESM3.",
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
            "case_file": "literature_cases/coulton_2024_pan_cancer_tam_covered.csv",
            "num_cases": len(case_files["coulton_2024_pan_cancer_tam_covered.csv"]),
        },
        {
            "paper": "Coulton 2024, pan-cancer TAM atlas",
            "doi": "10.1038/s41467-024-49885-8",
            "local_status": "PDF and Supplementary Data 1-7 available locally",
            "usable_for_benchmark": "yes - extended",
            "reason": "Broader cluster set includes monocyte-like and less-covered TAM labels; useful after reference brain is expanded.",
            "case_file": "literature_cases/coulton_2024_pan_cancer_tam_all.csv",
            "num_cases": len(case_files["coulton_2024_pan_cancer_tam_all.csv"]),
        },
        {
            "paper": "Wang 2023, prenatal macrophage atlas",
            "doi": "10.1016/j.cell.2023.08.019",
            "local_status": "Supplementary XLSX files available locally",
            "usable_for_benchmark": "yes - future domain",
            "reason": "Fifteen developmental macrophage subtypes with DEGs; not a tumor macrophage benchmark, so current reference coverage is partial.",
            "case_file": "literature_cases/wang_2023_prenatal_macrophage.csv",
            "num_cases": len(case_files["wang_2023_prenatal_macrophage.csv"]),
        },
        {
            "paper": "Li 2023, uveal melanoma macrophages",
            "doi": "10.1038/s12276-023-01115-9",
            "local_status": "PDF and two CSV supplements available locally",
            "usable_for_benchmark": "partial",
            "reason": "Four macrophage clusters have DEGs, but labels are M-C1 to M-C4; needs paper-derived semantic mapping before strict subtype scoring.",
            "case_file": "literature_cases/li_2023_uveal_melanoma_macrophage.csv",
            "num_cases": len(case_files["li_2023_uveal_melanoma_macrophage.csv"]),
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
    case_files = build_coulton_cases()
    case_files["wang_2023_prenatal_macrophage.csv"] = build_wang_cases()
    case_files["li_2023_uveal_melanoma_macrophage.csv"] = build_li_cases()

    for filename, cases in case_files.items():
        cases.to_csv(OUT_DIR / filename, index=False)

    inventory = build_inventory(case_files)
    inventory.to_csv(OUT_DIR / "literature_case_inventory.csv", index=False)

    print(f"Wrote literature benchmark cases to {OUT_DIR}")
    for filename, cases in case_files.items():
        print(f"  {filename}: {len(cases)} cases")
    print(f"  literature_case_inventory.csv: {len(inventory)} papers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
