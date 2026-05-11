---
id: macrophage_overview
category: myeloid
cell_types:
  - Macrophage
  - Monocyte-derived macrophage
  - Tissue-resident macrophage
  - Tumor-associated macrophage
trigger_markers:
  - CD68
  - CSF1R
  - C1QA
  - C1QB
  - C1QC
  - APOE
  - APOC1
  - LYZ
  - TYROBP
  - FCER1G
  - LST1
exclusion_markers:
  - CD3D
  - CD3E
  - MS4A1
  - CD79A
  - EPCAM
  - KRT19
  - PECAM1
sources:
  - Cheng et al. Cell 2021. https://doi.org/10.1016/j.cell.2021.01.010
  - Mulder et al. Immunity 2021. https://doi.org/10.1016/j.immuni.2021.07.007
  - Coulton et al. Nat Commun 2024. https://doi.org/10.1038/s41467-024-49885-8
---

# Macrophage Subtype Annotation Guide

## Overview

Macrophage subtype calls should be made as state-aware annotations, not as
fixed lineages. In scRNA-seq, macrophages often separate by tissue residency,
monocyte recruitment, interferon stimulation, inflammatory cytokines, lipid
handling, complement/phagocytosis, matrix remodeling, angiogenic programs,
stress, proliferation, and tumor conditioning.

For subtype-level CASSIA work, use the broad macrophage identity first, then
choose the most specific state supported by the top markers and tissue context.
Avoid forcing an M1/M2 label when the marker set supports a named scRNA-seq
state such as SPP1/AREG TAM, IFNG/IFN macrophage, inflammatory macrophage,
resident-like FOLR2/C1QC macrophage, lipid/TREM2 macrophage, heme macrophage,
or ECM-remodeling macrophage.

## Canonical Macrophage Identity

- Pan macrophage / mononuclear phagocyte markers: `CSF1R`, `CD68`, `LYZ`,
  `TYROBP`, `FCER1G`, `LST1`, `AIF1`, `ITGAM`, `CD14`, `MS4A7`.
- Phagocytic/complement macrophage markers: `C1QA`, `C1QB`, `C1QC`, `APOE`,
  `APOC1`, `CTSB`, `CTSZ`, `LIPA`, `GPNMB`, `ACP5`.
- Antigen-presenting macrophage markers: `HLA-DRA`, `HLA-DRB1`, `HLA-DPA1`,
  `HLA-DPB1`, `CD74`.
- Monocyte-like markers: `FCN1`, `S100A8`, `S100A9`, `VCAN`, `CCR2`, `IL1B`.

## Subtype Decision Points

### Monocyte-like vs macrophage-like

Call monocyte-like or recruited monocyte-derived macrophage when `FCN1`,
`VCAN`, `S100A8`, `S100A9`, `LST1`, `CCR2`, and inflammatory chemokines
dominate. Call macrophage-like when `APOE`, `APOC1`, `C1QA/B/C`, `RNASE1`,
`FOLR2`, `LIPA`, `CTSB`, `GPNMB`, `CD63`, `CD81`, `PLTP`, or `ACP5` dominate.

### Tissue-resident-like macrophage

Resident-like macrophages tend to show `FOLR2`, `SELENOP`, `SLC40A1`,
`STAB1`, `LYVE1`, `MRC1`, `CD163`, `C1QA/B/C`, `APOE`, and `RNASE1`, with
organ-specific additions such as `FABP4/MARCO/PPARG` in alveolar macrophages,
`TMEM119/P2RY12/SALL1` in microglia, and `CLEC4F/VSIG4/MARCO/TIMD4` in
Kupffer cells.

### Tumor-associated macrophage states

Tumor-associated macrophages are better described by programs than by M1/M2:
SPP1/AREG inflammatory angiogenic, APOE/TREM2 lipid-associated, IFNG/IFN
response, metallothionein, heme/iron handling, ECM-remodeling, inflammatory
IL1B/TNF, and resident-like FOLR2/C1QC states. Use `tam_pan_cancer.md` for
these calls.

### Inflammatory / interferon states

Use inflammatory or IFN labels when `CXCL9`, `CXCL10`, `GBP1`, `STAT1`,
`ISG15`, `IFIT1`, `IFITM1/3`, `MX1`, `IL1B`, `TNF`, `CXCL8`, `CCL3`, `CCL4`,
or `NLRP3` dominate. Use `inflammatory_interferon.md` for state separation.

## Common Pitfalls

- `SPP1` alone is not enough to call a tumor-associated macrophage; check for
  macrophage background markers and co-programs such as `AREG`, `TREM2`,
  `MMP9`, `VEGFA`, `LIPA`, `APOE`, or inflammatory chemokines.
- `C1QA/B/C` can mark resident-like, phagocytic, or complement-rich TAMs and
  is not automatically anti-inflammatory.
- `FOLR2`, `SELENOP`, and `SLC40A1` usually support resident-like or iron-
  handling macrophages, especially when paired with `CD163`, `MRC1`, `LYVE1`,
  or `STAB1`.
- Avoid M1/M2-only subtype outputs. If M1/M2 is useful, mention it as a broad
  bias after giving the more precise scRNA-seq state.
- Tumor, inflamed tissue, and dissociation stress can induce overlapping
  cytokine, chemokine, heat-shock, and interferon programs.

## Recommended Output Language

Use concise labels like:

- SPP1/AREG inflammatory angiogenic TAM
- IFNG/CXCL9 interferon-activated macrophage
- FOLR2/SELENOP resident-like macrophage
- APOE/TREM2 lipid-associated TAM
- C1QC complement/phagocytic macrophage
- IL1B/TNF inflammatory macrophage
- HMOX1/SLC40A1 heme-handling macrophage
- COL1A1/SPARC ECM-remodeling TAM
