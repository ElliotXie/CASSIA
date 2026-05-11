---
id: macrophage_resident_like
category: myeloid
cell_types:
  - Tissue-resident macrophage
  - Resident-like macrophage
  - Macrophage
  - Tumor-associated macrophage
trigger_markers:
  - C1QA
  - C1QB
  - C1QC
  - APOE
  - APOC1
  - RNASE1
  - FOLR2
  - SELENOP
  - SLC40A1
  - STAB1
  - LYVE1
  - MRC1
  - CD163
  - MARCO
  - FABP4
  - PPARG
  - TMEM119
  - P2RY12
  - CLEC4F
  - VSIG4
exclusion_markers:
  - CD3D
  - MS4A1
  - EPCAM
sources:
  - Mulder et al. Immunity 2021. https://doi.org/10.1016/j.immuni.2021.07.007
  - Cheng et al. Cell 2021. https://doi.org/10.1016/j.cell.2021.01.010
  - Coulton et al. Nat Commun 2024. https://doi.org/10.1038/s41467-024-49885-8
---

# Resident-like Macrophage Reference

## Overview

Resident-like macrophage annotation should combine conserved macrophage markers
with tissue context. Across datasets, macrophage identity is often supported by
`APOE`, `APOC1`, `C1QB`, `C1QC`, `RNASE1`, `ACP5`, `GPNMB`, `PLD3`, `CTSB`,
`PLTP`, `DAB2`, `CD63`, `CD81`, and `LIPA`. Resident-like or tissue-adapted
states add genes such as `FOLR2`, `SELENOP`, `SLC40A1`, `STAB1`, `LYVE1`,
`MRC1`, `CD163`, and organ-specific markers.

## Conserved Macrophage Signature

Mulder et al. report a cross-tissue macrophage signature including:

`APOE`, `APOC1`, `C1QB`, `C1QC`, `RNASE1`, `ACP5`, `GPNMB`, `PLD3`, `CTSB`,
`PLTP`, `CD9`, `PRDX1`, `CTSZ`, `DAB2`, `CD63`, `CD81`, `LIPA`, `GLUL`,
`SLCO2B1`, `CREG1`, `LGALS3`, `LAMP1`.

Use this signature to distinguish macrophages from classical monocytes and
from dendritic cells.

## Resident-like / Tissue-adapted Programs

### FOLR2/SELENOP/SLC40A1 macrophage

**Core markers:** `FOLR2`, `SELENOP`, `SLC40A1`, `F13A1`, `STAB1`, `RNASE1`,
`PLTP`, `LGMN`, `DAB2`, `MS4A6A`, `MS4A4A`, `CD163`, `MRC1`.

**Interpretation:** Resident-like, tissue-adapted, iron-handling macrophage.
In the Coulton atlas, this is close to `1_MetM2Mac`; in tumor contexts it can
be called FOLR2/SELENOP resident-like TAM if tumor-associated.

### C1QC complement/phagocytic macrophage

**Core markers:** `C1QA`, `C1QB`, `C1QC`, `APOE`, `APOC1`, `PLD4`, `CD74`,
`HLA-DPA1`, `HLA-DPB1`, `HLA-DRB5`, `CX3CR1`, `GPR34`.

**Interpretation:** Complement-rich phagocytic macrophage, often antigen-
presenting. In tumor references, this can be C1QC+ TAM; outside tumor it may
be a resident or resident-like macrophage state.

### LYVE1/MRC1/CD163 macrophage

**Core markers:** `LYVE1`, `MRC1`, `CD163`, `FOLR2`, `STAB1`, `SELENOP`,
`SLC40A1`, `C1QA`, `C1QB`, `C1QC`.

**Interpretation:** Perivascular/interstitial resident-like macrophage. Use
tissue labels when the organ is known, such as interstitial lung macrophage,
perivascular macrophage, or LYVE1+ resident-like TAM.

## Organ-specific Anchors

- Lung alveolar macrophage: `FABP4`, `MARCO`, `PPARG`, `MCEMP1`, `MRC1`,
  `INHBA`, `RETN` in human; `Siglecf` in mouse.
- Lung interstitial macrophage: `LYVE1`, `FOLR2`, `MRC1`, `CD163`, `SEPP1`
  or `SELENOP`.
- Brain microglia: `TMEM119`, `P2RY12`, `SALL1`, `CX3CR1`, `CSF1R`.
- Border-associated macrophage: `LYVE1`, `MRC1`, `CD163`, `FOLR2`, `MSR1`.
- Liver Kupffer cell: `CLEC4F`, `VSIG4`, `MARCO`, `TIMD4`, `CD5L`.
- Peritoneal macrophage: `GATA6`, `ICAM2`, `TIMD4`.
- Osteoclast / osteoclast-like macrophage: `CTSK`, `ACP5`, `MMP9`.

## Monocyte-vs-Resident Calls

Call a cluster monocyte-like when `FCN1`, `VCAN`, `S100A8`, `S100A9`,
`S100A12`, `CCR2`, and `CSF3R` dominate. Call resident-like macrophage when
the top markers shift toward `APOE/APOC1/C1QA/B/C/RNASE1/FOLR2/SELENOP`.
Mixed clusters can be described as "monocyte-derived macrophages transitioning
toward resident-like macrophage" if both programs are strong.

## Common Pitfalls

- `C1QA/B/C` is a macrophage/complement signal, not a plasma cell signal.
- `FOLR2` and `SELENOP` are more informative for resident-like macrophages
  than generic M2 labels.
- Alveolar macrophage markers are tissue-specific; do not call AlvMac in a
  non-lung dataset unless there is a plausible lung or metastatic-lung context.
- In tumor data, resident-like TAMs may coexist with SPP1/TREM2/IFN programs.
