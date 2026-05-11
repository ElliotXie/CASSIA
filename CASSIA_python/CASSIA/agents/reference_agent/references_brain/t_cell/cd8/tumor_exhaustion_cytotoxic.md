---
id: cd8_t_cell_tumor_exhaustion_cytotoxic
category: t_cell
parent: immune_t_cell_consensus
cell_types:
  - CD8+ T cell
  - Cytotoxic CD8 T cell
  - Exhausted CD8 T cell
  - Dysfunctional CD8 T cell
  - Tumor-infiltrating CD8 T cell
condition:
  - human tumor
  - chronic infection
  - inflamed tissue
  - peripheral blood
trigger_markers:
  - CD8A
  - CD8B
  - GZMK
  - GZMA
  - GZMB
  - GZMH
  - PRF1
  - NKG7
  - GNLY
  - CCL5
  - IFNG
  - PDCD1
  - LAG3
  - HAVCR2
  - TOX
  - CXCL13
  - ENTPD1
  - TIGIT
  - CTLA4
  - TCF7
  - SLAMF6
  - XCL1
  - XCL2
  - CD69
  - ITGAE
  - CXCR6
exclusion_markers:
  - CD4
  - FOXP3
  - IL2RA
  - MS4A1
  - CD79A
  - LST1
  - FCER1G
  - NCAM1
sources:
  - Zheng et al. Cell 2017. https://doi.org/10.1016/j.cell.2017.05.035
  - Guo et al. Nat Med 2018. https://doi.org/10.1038/s41591-018-0045-3
  - Andreatta et al. Nat Commun 2021. https://doi.org/10.1038/s41467-021-23324-4
  - Zheng et al. Science 2021. https://doi.org/10.1126/science.abe6474
  - Chu et al. Nat Med 2023. https://doi.org/10.1038/s41591-023-02371-y
---

# CD8 T Cell Tumor, Exhaustion, and Cytotoxic Reference

## Scope

Use this document for CD8 T-cell subclustering, especially tumor-infiltrating
or chronically stimulated T cells. It separates naive/central-memory,
GZMK/transitional, cytotoxic effector, precursor-exhausted, terminal exhausted,
and tissue-resident CD8 programs.

Do not use this document to label NK cells. Cytotoxic genes are shared between
CD8 T cells and NK cells, so require T-cell identity markers such as `CD3D`,
`CD3E`, `CD3G`, `TRAC`, or `TRBC1/2`.

## CD8 Identity Gate

**Required support:** `CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1`, or `TRBC2`.

**CD8 support:** `CD8A` and especially `CD8B`. `CD8A` alone is weaker because
it can be detected in NK-like and dendritic contexts.

**NK caution:** If `NCAM1`, `KLRF1`, `TYROBP`, `FCER1G`, `FCGR3A`, and low
TCR/CD3 genes dominate, prefer NK or NK-like lymphocyte instead of CD8 T cell.

## Consensus CD8 State Catalog

### CCR7/TCF7 naive-like or central-memory CD8 T cell

**Core markers:** `CCR7`, `SELL`, `TCF7`, `LEF1`, `IL7R`, `LTB`, `MAL`,
`NOSIP`, `PIK3IP1`.

**Interpretation:** Naive-like or central-memory CD8 program. In blood and
lymphoid tissue this is usually conventional memory/naive biology. In tumor
samples, evaluate whether the same `TCF7/LEF1` block is accompanied by
checkpoint markers that suggest precursor exhaustion.

**Suggested label:** `CCR7/TCF7 central-memory CD8 T cell`.

### GZMK/CCL5 transitional or effector-memory CD8 T cell

**Core markers:** `GZMK`, `CCL5`, `CXCR3`, `DUSP2`, `CST7`, `IL7R`, `GZMA`,
with variable `NKG7`.

**Interpretation:** Effector-memory, transitional, or pre-effector CD8 state.
Tumor T-cell atlases often place GZMK+ cells between memory-like and more
cytotoxic or dysfunctional programs.

**Suggested label:** `GZMK+ transitional CD8 T cell`.

**Pitfalls:** Do not call terminal cytotoxic effector unless `GZMB`, `GZMH`,
`PRF1`, `GNLY`, and `FGFBP2` are strong. Do not call exhausted unless
checkpoint and `TOX/CXCL13/ENTPD1` evidence is present.

### GNLY/PRF1 cytotoxic effector CD8 T cell

**Core markers:** `NKG7`, `GNLY`, `PRF1`, `GZMB`, `GZMH`, `GZMA`, `CTSW`,
`FGFBP2`, `KLRD1`, `KLRG1`, `CCL5`, `IFNG`.

**Interpretation:** Active cytotoxic effector CD8 program, common in antiviral,
inflammatory, and tumor immune responses.

**Suggested label:** `GNLY/PRF1 cytotoxic CD8 T cell`.

**Pitfalls:** NK cells share `NKG7/GNLY/PRF1/GZMB`. Require T-cell identity
and consider NK contamination when TCR/CD3 genes are weak.

### TCF7/PDCD1 precursor-exhausted CD8 T cell

**Core markers:** `TCF7`, `LEF1`, `IL7R`, `SLAMF6`, `PDCD1`, `TOX`, `XCL1`,
`XCL2`, `CXCR5`, with lower terminal checkpoint and cytotoxic genes.

**Interpretation:** Stem-like or precursor-exhausted CD8 state, most relevant
in tumor or chronic antigen exposure. Andreatta et al. support quiescent or
precursor-like dysfunctional states that can coexist with terminal exhausted
programs.

**Suggested label:** `TCF7+ precursor-exhausted CD8 T cell`.

**Pitfalls:** In healthy blood, the same `TCF7/LEF1/IL7R` block usually means
naive or central-memory CD8 T cell. Use tumor/chronic-stimulation context and
checkpoint support before calling precursor exhaustion.

### TOX/CXCL13 exhausted or dysfunctional CD8 T cell

**Core markers:** `PDCD1`, `LAG3`, `HAVCR2`, `TIGIT`, `CTLA4`, `TOX`,
`CXCL13`, `ENTPD1`, `LAYN`, `TNFRSF9`, `GZMB`, `PRF1`.

**Interpretation:** Chronic antigen-experienced exhausted or dysfunctional CD8
T-cell program. In tumor references, `CXCL13`, `ENTPD1`, `LAYN`, and `TOX`
support tumor-reactive exhausted TIL interpretation.

**Suggested label:** `TOX/CXCL13 exhausted CD8 T cell`.

**Pitfalls:** `PDCD1` alone is activation-compatible. Require a multi-marker
checkpoint program plus `TOX`, `CXCL13`, `ENTPD1`, `LAYN`, or tumor context.
`CTLA4` and `TIGIT` also occur in Tregs.

### CD69/ITGAE tissue-resident CD8 T cell

**Core markers:** `CD69`, `ITGAE`, `CXCR6`, `ITGA1`, `ZNF683`, `XCL1`,
`XCL2`, with reduced circulation markers such as `SELL` and `S1PR1`.

**Interpretation:** Tissue-resident memory CD8 state. Use this as a state
modifier when cytotoxic, memory, or exhausted markers are also present.

**Suggested label:** `CD69/ITGAE tissue-resident memory CD8 T cell`.

**Pitfalls:** `CD69` alone is an activation marker. Prefer a residency call
when `ITGAE`, `CXCR6`, `ITGA1`, or `ZNF683` co-occur.

## Practical Distinctions

- `CCR7/SELL/TCF7/LEF1` -> naive-like or central-memory CD8 unless tumor checkpoint markers suggest precursor exhaustion.
- `GZMK/CCL5/CXCR3` -> transitional or effector-memory CD8.
- `GNLY/PRF1/GZMB/GZMH/FGFBP2` -> cytotoxic effector CD8, with NK caution.
- `TCF7/PDCD1/SLAMF6/TOX` -> precursor-exhausted CD8 in tumor or chronic antigen context.
- `PDCD1/LAG3/HAVCR2/TOX/CXCL13/ENTPD1` -> exhausted or dysfunctional CD8.
- `CD69/ITGAE/CXCR6` -> tissue-resident modifier.

## Recommended Output Language

- CCR7/TCF7 central-memory CD8 T cell
- GZMK+ transitional CD8 T cell
- GNLY/PRF1 cytotoxic CD8 T cell
- TCF7+ precursor-exhausted CD8 T cell
- TOX/CXCL13 exhausted CD8 T cell
- CD69/ITGAE tissue-resident memory CD8 T cell
