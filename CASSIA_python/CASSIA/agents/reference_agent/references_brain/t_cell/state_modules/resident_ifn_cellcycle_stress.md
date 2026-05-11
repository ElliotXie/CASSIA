---
id: t_cell_state_modules_resident_ifn_cellcycle_stress
category: t_cell
parent: immune_t_cell_consensus
cell_types:
  - Tissue-resident memory T cell
  - Interferon-stimulated T cell
  - Proliferating T cell
  - Stress-response T cell
condition:
  - human immune tissue
  - inflamed tissue
  - tumor microenvironment
  - peripheral blood
trigger_markers:
  - CD69
  - ITGAE
  - CXCR6
  - ITGA1
  - ZNF683
  - XCL1
  - XCL2
  - ISG15
  - IFIT1
  - IFIT2
  - IFIT3
  - IFITM1
  - IFITM3
  - MX1
  - OAS1
  - OAS2
  - RSAD2
  - STAT1
  - IRF7
  - MKI67
  - TOP2A
  - STMN1
  - TYMS
  - UBE2C
  - PCLAF
  - CENPF
  - HSPA1A
  - HSPA1B
  - HSPA6
  - DNAJB1
  - HSP90AA1
  - HSPB1
  - BAG3
  - FOS
  - JUN
  - DUSP1
exclusion_markers:
  - MS4A1
  - CD79A
  - LST1
  - FCER1G
  - C1QA
  - C1QB
  - C1QC
  - EPCAM
sources:
  - Szabo et al. Nat Commun 2019. https://doi.org/10.1038/s41467-019-12464-3
  - Andreatta et al. Nat Commun 2021. https://doi.org/10.1038/s41467-021-23324-4
  - Zheng et al. Science 2021. https://doi.org/10.1126/science.abe6474
  - Chu et al. Nat Med 2023. https://doi.org/10.1038/s41591-023-02371-y
---

# T Cell State Modules: Residency, IFN, Cell Cycle, and Stress

## Scope

Use this document as a state-module reference, not as a primary lineage
ontology. These programs modify a T-cell subtype call and should usually be
reported as suffixes or cautions: tissue-resident memory, interferon-stimulated,
proliferating/cycling, and heat-shock or stress-response.

Example output:

`TOX/CXCL13 exhausted CD8 T cell with MKI67+ proliferating state`

Avoid output:

`MKI67+ cell` or `stress cell` without T-cell identity and subtype context.

## State Module Catalog

### CD69/ITGAE/CXCR6 tissue-resident memory module

**Core markers:** `CD69`, `ITGAE`, `CXCR6`, `ITGA1`, `ZNF683`, `XCL1`,
`XCL2`, with reduced circulation markers such as `SELL` and `S1PR1`.

**Interpretation:** Tissue-resident memory program. This module can occur in
CD8, CD4, MAIT, gamma-delta, and exhausted T-cell states.

**Recommended use:** Add `tissue-resident memory` or `TRM-like` to the best
underlying T-cell subtype.

**Pitfalls:** `CD69` alone is an early activation marker. Require `ITGAE`,
`CXCR6`, `ITGA1`, or `ZNF683` for stronger residency evidence.

### ISG15/IFIT interferon-stimulated module

**Core markers:** `ISG15`, `IFIT1`, `IFIT2`, `IFIT3`, `IFITM1`, `IFITM3`,
`MX1`, `OAS1`, `OAS2`, `RSAD2`, `STAT1`, `IRF7`, `IFI6`.

**Interpretation:** Type-I interferon or antiviral-response T-cell state. This
can reflect infection, inflammation, tumor nucleic-acid sensing, treatment, or
sample-wide interferon exposure.

**Recommended use:** Add `ISG15+ interferon-stimulated` to the best T-cell
subtype.

**Pitfalls:** IFN signatures are not lineage-specific and can occur in many
lineages. Confirm T-cell identity before applying this module.

### MKI67/TOP2A proliferating module

**Core markers:** `MKI67`, `TOP2A`, `STMN1`, `TYMS`, `UBE2C`, `PCLAF`,
`CENPF`, `HMGB2`, `TUBA1B`, `TUBB`.

**Interpretation:** Cell-cycle or proliferating state. This may occur in
expanding Tregs, exhausted T cells, cytotoxic T cells, or activated helper T
cells depending on the accompanying markers.

**Recommended use:** Add `MKI67+ proliferating` to the underlying subtype.

**Pitfalls:** Cell-cycle genes are not lineage markers. Do not assign
`proliferating T cell` unless CD3/TCR genes support T-cell identity.

### HSPA/DNAJB1 stress-response module

**Core markers:** `HSPA1A`, `HSPA1B`, `HSPA6`, `DNAJB1`, `HSP90AA1`,
`HSPB1`, `BAG3`, `FOS`, `JUN`, `DUSP1`, `IER2`.

**Interpretation:** Heat-shock or stress-response T-cell state. Chu et al.
describe a pan-cancer T-cell stress-response state, but the same genes can
also reflect dissociation, sample handling, or broad activation stress.

**Recommended use:** Add `HSPA/DNAJB1 stress-response` as a caution or suffix
when T-cell identity is otherwise clear.

**Pitfalls:** If stress genes appear across many unrelated lineages, prefer a
technical-stress caution rather than a biological subtype.

## Practical Distinctions

- `CD69/ITGAE/CXCR6/ZNF683` -> tissue-resident memory module.
- `ISG15/IFIT/MX1/OAS/RSAD2` -> interferon-stimulated module.
- `MKI67/TOP2A/STMN1/TYMS/UBE2C` -> proliferating module.
- `HSPA1A/HSPA1B/HSPA6/DNAJB1/FOS/JUN` -> stress-response module.

## Recommended Output Language

- CD69/ITGAE tissue-resident memory T cell
- ISG15+ interferon-stimulated T cell
- MKI67+ proliferating T cell
- HSPA/DNAJB1 stress-response T cell
- exhausted CD8 T cell with tissue-resident memory module
- effector Treg with MKI67+ proliferating state
