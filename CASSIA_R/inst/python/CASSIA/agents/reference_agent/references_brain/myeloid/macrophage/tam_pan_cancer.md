---
id: macrophage_tam_pan_cancer
category: myeloid
cell_types:
  - Tumor-associated macrophage
  - TAM
  - Macrophage
trigger_markers:
  - SPP1
  - AREG
  - TREM2
  - APOE
  - APOC1
  - GPNMB
  - LIPA
  - MMP9
  - VEGFA
  - CXCL9
  - CXCL10
  - FOLR2
  - SLC40A1
  - HMOX1
  - COL1A1
  - COL1A2
  - SPARC
exclusion_markers:
  - EPCAM
  - KRT8
  - KRT18
  - KRT19
sources:
  - Cheng et al. Cell 2021. https://doi.org/10.1016/j.cell.2021.01.010
  - Coulton et al. Nat Commun 2024. https://doi.org/10.1038/s41467-024-49885-8
---

# Pan-Cancer TAM Subtype Reference

## Overview

Pan-cancer TAM references are most useful when macrophages are being
subclustered inside tumor or inflamed tissue. Cheng et al. profiled
tumor-infiltrating myeloid cells across 15 cancer types and emphasized
SPP1+, C1QC+, ISG15+, FN1+, LYVE1+, NLRP3+, and INHBA+ macrophage/TAM
patterns. Coulton et al. built a larger pan-cancer TAM atlas across 17 tumor
types and named recurrent macrophage clusters such as `SPP1AREGMac`,
`IFNGMac`, `MetalloMac`, `HemeMac`, `ECMMac`, `InflamMac`, `AngioMac`,
`MetM2Mac`, `AlvMac`, and `ICIMac`.

Use this document to convert vague "TAM/M2 macrophage" outputs into specific
subtype labels tied to marker programs.

## Coulton 2024 Author Cluster Labels

These label-to-marker mappings come from Coulton et al. Nat Commun 2024
Supplementary Data 6. When a query marker set strongly matches one of these
programs, preserve the author label in the annotation reason or subtype label
so downstream users can trace the call back to the paper.

- `1_MetM2Mac`: `SELENOP`, `SLC40A1`, `F13A1`, `RNASE1`, `FOLR2`, `STAB1`,
  `LGMN`, `DAB2`.
- `2_C3Mac`: `C3`, `PLD4`, `RGS1`, `HLA-DPA1`, `CX3CR1`, `CD74`,
  `HLA-DPB1`, `HLA-DRB5`.
- `4_ICIMac2`: `APOE`, `APOC1`, `CCL18`, `GPNMB`, `CTSD`, `LIPA`,
  `LGMN`, `PLA2G7`, `TREM2`.
- `5_StressMac`: `HSPA6`, `HSPA1B`, `HSPA1A`, `DNAJB1`, `HSPB1`, `BAG3`,
  `HSPH1`, `HSP90AA1`.
- `6_SPP1AREGMac`: `CCL20`, `CXCL3`, `IL1B`, `CXCL2`, `CXCL8`, `EREG`,
  `G0S2`, `TIMP1`, `CXCL1`.
- `8_IFNGMac`: `CXCL9`, `CXCL10`, `GBP1`, `MMP9`, `GBP5`, `WARS1`,
  `STAT1`, `SLAMF7`.
- `9_AngioMac`: `AREG`, `THBS1`, `EREG`, `NAMPT`, `IL1B`, `GPR183`,
  `NR4A3`.
- `10_InflamMac`: `CCL3L3`, `CCL4L2`, `CXCL8`, `IL1B`, `TNF`, `CCL4`,
  `CCL3`, `CXCL2`.
- `11_MetalloMac`: `MT1G`, `MT1X`, `MT2A`, `MT1E`, `MT1H`, `MT1F`, `MT1M`,
  `MIF`, `SPP1`.
- `17_IFNMac3`: `ISG15`, `CXCL10`, `IFIT1`, `IFIT2`, `CCL8`, `IFIT3`,
  `IFITM1`, `MX1`, `RSAD2`.
- `18_ECMMac`: `COL1A2`, `COL1A1`, `COL3A1`, `IGFBP7`, `SPARC`, `MGP`,
  `LUM`, `DCN`, `POSTN`.
- `21_HemeMac`: `HMOX1`, `SLC40A1`, `CCL18`, `CD163`, `LGMN`, `CTSB`,
  `CTSL`.

## High-Confidence TAM Programs

### SPP1/AREG inflammatory angiogenic TAM

**Core markers:** `SPP1`, `AREG`, `EREG`, `CXCL3`, `CXCL2`, `CXCL8`, `CXCL1`,
`CCL20`, `IL1B`, `TIMP1`, `SOD2`, `IL1RN`, `MMP9`, `VEGFA`.

**Interpretation:** This is a tumor-conditioned inflammatory and angiogenic
macrophage program. In the Coulton atlas, `6_SPP1AREGMac` is marked by
epithelial growth factor ligands and inflammatory chemokines. Cheng et al.
reported SPP1+ TAMs as angiogenesis-associated and often tumor-enriched. Use
labels such as "SPP1/AREG inflammatory angiogenic TAM" or "SPP1+ pro-
angiogenic TAM" rather than generic M2.

**Distinguish from:** Monocyte-like inflammation if `FCN1`, `VCAN`, `S100A8`,
and `S100A9` dominate; ECM-remodeling TAM if collagen genes dominate.

### APOE/TREM2 lipid-associated TAM

**Core markers:** `TREM2`, `APOE`, `APOC1`, `GPNMB`, `LIPA`, `CTSD`, `LGMN`,
`PLA2G7`, `ACP5`, `PLD3`, `PSAP`.

**Interpretation:** Lipid/phagolysosomal TAM state, often immunosuppressive or
tumor-conditioned. In Coulton, `4_ICIMac2` is strongly APOE/APOC1/TREM2/GPNMB.
Use "APOE/TREM2 lipid-associated TAM" when lipid handling and lysosomal genes
are present with macrophage markers.

**Distinguish from:** Resident-like FOLR2 macrophage when `FOLR2`, `SELENOP`,
`SLC40A1`, and `STAB1` dominate without strong `TREM2/GPNMB/LIPA`.

### FOLR2/SELENOP resident-like TAM

**Core markers:** `FOLR2`, `SELENOP`, `SLC40A1`, `F13A1`, `STAB1`, `RNASE1`,
`CD163`, `MRC1`, `LYVE1`, `C1QA`, `C1QB`, `C1QC`, `APOE`.

**Interpretation:** Resident-like or tissue-adapted TAM/macrophage program.
Coulton `1_MetM2Mac` is enriched for `FOLR2`, `SELENOP`, `SLC40A1`, `STAB1`;
Cheng discusses LYVE1+ resident tissue macrophages and C1QC+ TAMs with weaker
connectivity to CD14+ monocytes.

**Distinguish from:** Complement/phagocytic C1QC macrophages when `C1QA/B/C`
and MHC genes dominate; heme macrophages when `HMOX1` and `SLC40A1` dominate.

### C1QC complement/phagocytic TAM

**Core markers:** `C1QA`, `C1QB`, `C1QC`, `APOE`, `APOC1`, `CD74`, `HLA-DRA`,
`HLA-DPA1`, `HLA-DPB1`, `CX3CR1`, `PLD4`, `GPR34`.

**Interpretation:** Complement-rich, phagocytic, antigen-presenting macrophage.
Cheng described C1QC+ TAMs as a major tumor-enriched macrophage pattern and
reported higher phagocytosis scores than SPP1+ TAMs.

**Distinguish from:** B/plasma cells if immunoglobulin genes dominate; DCs if
`CLEC9A`, `XCR1`, `CD1C`, `FCER1A`, or `LILRA4` dominate.

### IFNG / CXCL9 TAM

**Core markers:** `CXCL9`, `CXCL10`, `GBP1`, `GBP5`, `STAT1`, `WARS1`,
`SLAMF7`, `LGALS2`, `VAMP5`, `CALHM6`, `MMP9`.

**Interpretation:** IFN-gamma-activated TAM/macrophage, often linked to T cell
interaction and antigen presentation. Coulton `8_IFNGMac` is the clearest
program. Use "IFNG/CXCL9 interferon-activated TAM" when `CXCL9/CXCL10/GBP`
genes dominate.

### Metallothionein macrophage

**Core markers:** `MT1G`, `MT1X`, `MT2A`, `MT1E`, `MT1H`, `MT1F`, `MT1M`,
`MIF`, `SPP1`, `LDHA`, `LGALS1`.

**Interpretation:** Metal-ion/stress-associated TAM state. Coulton
`11_MetalloMac` is a gold-standard signature. Do not label this as a generic
SPP1 TAM when metallothioneins dominate.

### Heme / iron-handling macrophage

**Core markers:** `HMOX1`, `SLC40A1`, `HAMP`, `CD163`, `CCL18`, `LGMN`,
`STAB1`, `SELENOP`, `CTSB`, `CTSL`.

**Interpretation:** Heme-processing or iron-handling macrophage. Coulton
`21_HemeMac` is a gold-standard signature. This may overlap with resident-like
macrophage signatures through `SLC40A1`, `CD163`, and `STAB1`.

### ECM-remodeling TAM

**Core markers:** `COL1A1`, `COL1A2`, `COL3A1`, `SPARC`, `COL6A1`, `COL6A2`,
`COL6A3`, `LUM`, `DCN`, `POSTN`, `BGN`, `CALD1`.

**Interpretation:** Matrix-remodeling TAM-like cluster from the Coulton atlas,
not a fibroblast by default if macrophage markers are present. However, because
many ECM genes are stromal, check for `PTPRC`, `LST1`, `TYROBP`, `C1QA/B/C`,
`APOE`, or `CD68` before calling it macrophage.

## Minimal Marker-to-Subtype Rules

- `SPP1 + AREG/EREG + CXCL1/2/3/8` -> SPP1/AREG inflammatory angiogenic TAM.
- `TREM2 + APOE/APOC1 + GPNMB/LIPA/CTSD` -> APOE/TREM2 lipid-associated TAM.
- `FOLR2 + SELENOP/SLC40A1/STAB1 + CD163/MRC1/LYVE1` -> resident-like TAM.
- `C1QA/B/C + HLA-DRA/CD74 + APOE` -> C1QC complement/phagocytic TAM.
- `CXCL9/CXCL10 + GBP1/GBP5 + STAT1` -> IFNG/CXCL9 TAM.
- `MT1G/MT1X/MT2A` dominant -> metallothionein macrophage.
- `HMOX1 + SLC40A1 + CD163/HAMP` -> heme/iron-handling macrophage.
- `COL1A1/COL1A2/SPARC` plus macrophage markers -> ECM-remodeling TAM.

## Common Pitfalls

- `COL1A1/COL1A2/SPARC` without macrophage markers is more likely fibroblast
  or stromal contamination than an ECM TAM.
- `SPP1` can also appear in epithelial cells, osteoclasts, fibroblasts, and
  some dendritic/myeloid states; require macrophage context.
- Cheng et al. showed SPP1 and C1QC can be mutually exclusive in some cancers
  but co-expressed in others. Do not assume they always define separate cells.
- "M2-like TAM" is too vague for subtype output when FOLR2, TREM2, SPP1, C1QC,
  heme, IFN, or ECM programs can be named directly.
