---
id: macrophage_inflammatory_interferon
category: myeloid
cell_types:
  - Inflammatory macrophage
  - Interferon-stimulated macrophage
  - IFNG macrophage
  - Tumor-associated macrophage
  - Macrophage
trigger_markers:
  - CXCL9
  - CXCL10
  - GBP1
  - GBP5
  - STAT1
  - ISG15
  - IFIT1
  - IFITM1
  - IFITM3
  - MX1
  - IL1B
  - TNF
  - CXCL8
  - CCL3
  - CCL4
  - NLRP3
  - HSPA1A
  - HSPA1B
  - HSPA6
  - DNAJB1
exclusion_markers:
  - CD3D
  - CD79A
  - EPCAM
sources:
  - Cheng et al. Cell 2021. https://doi.org/10.1016/j.cell.2021.01.010
  - Coulton et al. Nat Commun 2024. https://doi.org/10.1038/s41467-024-49885-8
---

# Inflammatory and Interferon Macrophage Reference

## Overview

Inflammatory and interferon macrophage states are common in tumor, infection,
autoimmune, and dissociation-stressed samples. They should be separated into
IFN-gamma/chemokine, type-I IFN/ISG, IL1B/TNF inflammatory, inflammasome-like,
and heat-shock/stress programs when possible.

## IFNG/CXCL9 Interferon-activated Macrophage

**Core markers:** `CXCL9`, `CXCL10`, `GBP1`, `GBP5`, `STAT1`, `WARS1`,
`VAMP5`, `SLAMF7`, `LGALS2`, `CALHM6`, `MMP9`.

**Interpretation:** IFN-gamma-activated macrophage/TAM, often associated with
T cell interaction, antigen presentation, and chemokine recruitment. In the
Coulton atlas, `8_IFNGMac` is marked by this program.

**Suggested label:** IFNG/CXCL9 interferon-activated macrophage or IFNG TAM.

## Type-I IFN / ISG Macrophage

**Core markers:** `ISG15`, `IFIT1`, `IFIT2`, `IFIT3`, `IFITM1`, `IFITM2`,
`IFITM3`, `MX1`, `OAS1`, `OAS2`, `RSAD2`, `CXCL10`, `TNFSF10`, `CASP1`,
`CASP4`.

**Interpretation:** Type-I IFN or antiviral response macrophage. Cheng et al.
reported ISG15+ TAMs with interferon-induced proteins and death/pyroptosis
regulators, with an M1-like bias. Coulton includes multiple IFN macrophage
clusters, including `IFNMac` and gold-standard `IFNMac3/IFNMac4` signatures.

**Suggested label:** ISG15+ interferon-stimulated macrophage/TAM.

## IL1B/TNF Inflammatory Macrophage

**Core markers:** `IL1B`, `TNF`, `CXCL8`, `CXCL1`, `CXCL2`, `CXCL3`, `CCL3`,
`CCL4`, `CCL3L3`, `CCL4L2`, `DUSP2`, `IER3`, `BCL2A1`, `NFKBIA`.

**Interpretation:** Acute inflammatory macrophage or monocyte-derived
macrophage state. Coulton `10_InflamMac` is marked by inflammatory cytokines
and chemokines. If `FCN1`, `VCAN`, `S100A8`, and `S100A9` dominate, call it
inflammatory monocyte-like macrophage rather than mature resident macrophage.

**Suggested label:** IL1B/TNF inflammatory macrophage.

## NLRP3 / Inflammasome-like Macrophage

**Core markers:** `NLRP3`, `IL1B`, `CASP1`, `CASP4`, `CXCL8`, `TNF`,
`BCL2A1`, `SOD2`, `NFKBIA`.

**Interpretation:** Inflammasome-related inflammatory macrophage state.
Use this call when `NLRP3` and inflammasome/caspase genes are present with
macrophage markers.

## Heat-shock / Stress Macrophage

**Core markers:** `HSPA6`, `HSPA1B`, `HSPA1A`, `DNAJB1`, `HSPB1`, `BAG3`,
`HSPH1`, `HSP90AA1`, `ZFAND2A`, `HSPD1`, `HSPE1`, `IER5`.

**Interpretation:** Heat-shock/stress macrophage state. Coulton `5_StressMac`
is a gold-standard signature. Treat this as a state label and consider
dissociation or sample-processing stress if it appears broadly across cell
types.

## Practical Distinctions

- `CXCL9/CXCL10/GBP1/STAT1` -> IFNG/CXCL9 macrophage.
- `ISG15/IFIT/MX1/OAS/RSAD2` -> type-I IFN or ISG15 macrophage.
- `IL1B/TNF/CXCL8/CCL3/CCL4` -> inflammatory macrophage.
- `NLRP3/CASP1/CASP4/IL1B` -> inflammasome-like macrophage.
- `HSPA1A/HSPA1B/HSPA6/DNAJB1` -> heat-shock/stress macrophage.

## Common Pitfalls

- IFN markers can occur in many immune cells. Require macrophage markers such
  as `LST1`, `TYROBP`, `FCER1G`, `C1QA/B/C`, `APOE`, or `CD68`.
- Inflammatory monocytes and macrophages overlap. `FCN1/VCAN/S100A8/S100A9`
  point toward monocyte-like cells; `APOE/C1QA/B/C/RNASE1/FOLR2` point toward
  macrophages.
- Heat-shock programs may reflect technical stress rather than a stable
  biological subtype.
