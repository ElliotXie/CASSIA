---
id: immune_t_cell_consensus
category: t_cell
cell_types:
  - T cell
  - Immune T cell
  - CD4+ T cell
  - CD8+ T cell
  - Tumor-infiltrating T cell
  - Regulatory T cell
  - T follicular helper cell
  - Tissue-resident memory T cell
  - MAIT cell
  - Gamma-delta T cell
condition:
  - human immune tissue
  - inflamed tissue
  - tumor microenvironment
  - peripheral blood
trigger_markers:
  - CD3D
  - CD3E
  - CD3G
  - TRAC
  - TRBC1
  - TRBC2
  - CD4
  - CD8A
  - CD8B
  - IL7R
  - CCR7
  - SELL
  - TCF7
  - LEF1
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
  - TIGIT
  - CTLA4
  - FOXP3
  - IL2RA
  - CCR8
  - CXCR5
  - BCL6
  - ICOS
  - CD69
  - ITGAE
  - CXCR6
  - TRDC
  - TRGC1
  - TRGC2
  - KLRB1
  - SLC4A10
  - TRAV1-2
  - ZBTB16
  - MKI67
  - TOP2A
  - HSPA1A
  - HSPA1B
exclusion_markers:
  - MS4A1
  - CD79A
  - CD79B
  - MZB1
  - JCHAIN
  - LST1
  - FCER1G
  - C1QA
  - C1QB
  - C1QC
  - EPCAM
  - KRT8
  - KRT18
  - KRT19
  - PECAM1
  - COL1A1
sources:
  - Szabo et al. Nat Commun 2019. https://doi.org/10.1038/s41467-019-12464-3
  - Zheng et al. Cell 2017. https://doi.org/10.1016/j.cell.2017.05.035
  - Guo et al. Nat Med 2018. https://doi.org/10.1038/s41591-018-0045-3
  - Andreatta et al. Nat Commun 2021. https://doi.org/10.1038/s41467-021-23324-4
  - Zheng et al. Science 2021. https://doi.org/10.1126/science.abe6474
  - Chu et al. Nat Med 2023. https://doi.org/10.1038/s41591-023-02371-y
  - Zemmour et al. Nat Immunol 2018. https://doi.org/10.1038/s41590-018-0051-0
  - Garner et al. Nat Immunol 2023. https://doi.org/10.1038/s41590-023-01575-1
---

# Immune T Cell Consensus Reference

## Scope

This document is a consensus reference for human immune T-cell annotation and
subclustering in scRNA-seq data. It is designed for mixed T-cell populations in
blood, inflamed tissue, and tumor microenvironments. Use it to choose
state-aware T-cell labels rather than broad labels such as "T cell",
"activated T cell", or "cytotoxic lymphocyte" when marker evidence supports a
more specific call.

Good output pattern:

`CXCL13/PDCD1 exhausted CD8 T cell (tumor-infiltrating T-cell atlas-like)`

Bad output pattern:

`Exhausted` with no CD8/T-cell identity, marker program, or context.

## Evidence Model

- **Consensus name**: preferred CASSIA subtype wording.
- **Core markers**: marker program that should drive the call.
- **Context**: blood, tissue, tumor, activation, or technical-state cautions.
- **Paper aliases**: author-defined states can be mentioned as evidence, not
  as universal cell ontology.
- **Pitfalls**: common overcalls and contaminating lineage programs.

Szabo et al. Nat Commun 2019 provides a healthy human T-cell reference across
lung, lymph node, bone marrow, blood, resting states, and activation states.
Zheng et al. Cell 2017 and Guo et al. Nat Med 2018 provide early single-cell
tumor-infiltrating T-cell maps in HCC and NSCLC. Andreatta et al. Nat Commun
2021 and Zheng et al. Science 2021 support cross-study or pan-cancer T-cell
reference states. Chu et al. Nat Med 2023 adds a large pan-cancer T-cell atlas
with stress-response, Tfh, Treg, and proliferative heterogeneity. Zemmour et
al. Nat Immunol 2018 supports regulatory T-cell heterogeneity. Garner et al.
Nat Immunol 2023 supports MAIT-cell transcriptional, tissue, functional, and
clonal diversity.

## Identity Gate

Before subtype annotation, confirm that the cluster is truly T lineage:

- Strong T-cell identity: `CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1`, `TRBC2`.
- CD4 lineage support: `CD4`, `IL7R`, `CCR7`, `LTB`, `ANXA1`.
- CD8 lineage support: `CD8A`, `CD8B`, `GZMK`, `CCL5`, `NKG7`, `PRF1`.
- Do not call a T-cell subtype from `CD4`, `CD8A`, `NKG7`, or cytotoxic genes
  alone. Monocytes can express `CD4`; NK cells share `NKG7`, `GNLY`, `PRF1`,
  `GZMB`; dendritic cells may show low `CD8A`.

## Consensus State Catalog

### CCR7/TCF7 naive-like or central-memory T cell

**Core markers:** `CCR7`, `SELL`, `TCF7`, `LEF1`, `IL7R`, `LTB`, `MAL`,
`NOSIP`, `PIK3IP1`.

**Interpretation:** Naive-like or central-memory T cell program. Use CD4 or
CD8 as a prefix when lineage markers are clear, for example
`CCR7/TCF7 naive-like CD4 T cell` or `CCR7/TCF7 central-memory CD8 T cell`.

**Pitfalls:** `TCF7` can also mark precursor exhausted T cells in tumors. If
`PDCD1`, `TOX`, `XCL1`, `CXCL13`, or `ENTPD1` are present, consider a
precursor-exhausted or dysfunctional tumor T-cell state instead of naive.

### GZMK/CCL5 effector-memory or pre-effector T cell

**Core markers:** `GZMK`, `CCL5`, `CXCR3`, `DUSP2`, `CST7`, `IL7R`, `CD44`,
with variable `NKG7` and low to moderate cytotoxic genes.

**Interpretation:** Effector-memory, transitional, or pre-effector T-cell
state. In tumor TIL references this can resemble pre-dysfunctional or
pre-exhausted CD8 states when it occurs with weak checkpoint expression.

**Suggested label:** `GZMK+ effector-memory CD8 T cell` or
`GZMK+ transitional T cell`.

**Pitfalls:** Do not call terminal cytotoxic effector if `GZMB`, `PRF1`,
`GNLY`, and `FGFBP2` are weak. Do not call exhausted unless the checkpoint and
TOX/CXCL13 program is convincing.

### GNLY/PRF1/GZMB cytotoxic effector T cell

**Core markers:** `NKG7`, `GNLY`, `PRF1`, `GZMB`, `GZMH`, `GZMA`, `CTSW`,
`FGFBP2`, `KLRD1`, `KLRG1`, `CCL5`.

**Interpretation:** Cytotoxic effector T-cell program, usually CD8 when
`CD8A/CD8B` and CD3/TCR genes are present. This program can be anti-viral,
anti-tumor, or inflammatory depending on tissue context.

**Suggested label:** `GNLY/PRF1 cytotoxic CD8 T cell`.

**Pitfalls:** NK cells share this program. Require T-cell identity
(`CD3D/E/G`, `TRAC/TRBC`) and check NK lineage markers such as `NCAM1`, `KLRF1`,
`FCGR3A`, `TYROBP`, and `FCER1G`.

### TOX/CXCL13 exhausted or dysfunctional CD8 T cell

**Core markers:** `PDCD1`, `LAG3`, `HAVCR2`, `TIGIT`, `CTLA4`, `TOX`,
`CXCL13`, `ENTPD1`, `LAYN`, `TNFRSF9`, `GZMB`, `PRF1`.

**Interpretation:** Chronic antigen-experienced dysfunctional or exhausted
T-cell program, common in tumor and chronic infection. In tumor references,
`CXCL13`, `ENTPD1`, `LAYN`, and `TOX` strengthen tumor-reactive exhausted TIL
interpretation.

**Suggested label:** `TOX/CXCL13 exhausted CD8 T cell` or
`PDCD1/LAG3 dysfunctional CD8 T cell`.

**Pitfalls:** `PDCD1` alone is activation-compatible and is not enough for
exhaustion. Require a multi-marker checkpoint program plus `TOX`, `CXCL13`,
`ENTPD1`, or tumor context. `CTLA4` and `TIGIT` also occur in Tregs.

### TCF7/PDCD1 precursor-exhausted or stem-like exhausted T cell

**Core markers:** `TCF7`, `LEF1`, `IL7R`, `SLAMF6`, `PDCD1`, `TOX`, `XCL1`,
`CXCR5`, with lower `HAVCR2` and lower terminal cytotoxic genes than terminal
exhausted cells.

**Interpretation:** Precursor-exhausted, stem-like exhausted, or quiescent
dysfunctional tumor T-cell state. Andreatta et al. describe rare quiescent
precursor states that can coexist with terminal exhausted CD8 TILs.

**Suggested label:** `TCF7+ precursor-exhausted CD8 T cell`.

**Pitfalls:** In non-tumor blood samples, `TCF7/LEF1/IL7R` usually supports
naive or central memory rather than precursor exhaustion. Context matters.

### CD69/ITGAE/CXCR6 tissue-resident memory T cell

**Core markers:** `CD69`, `ITGAE`, `CXCR6`, `ITGA1`, `ZNF683`, `XCL1`,
`XCL2`, with reduced circulation markers such as `SELL` and `S1PR1`.

**Interpretation:** Tissue-resident memory T-cell program. It can occur in
CD8, CD4, MAIT, or other T-cell subsets and should usually be used as a state
suffix if another subtype is clear.

**Suggested label:** `CD69/ITGAE tissue-resident memory CD8 T cell`.

**Pitfalls:** `CD69` is also an early activation marker. Require a residency
program such as `ITGAE`, `CXCR6`, `ITGA1`, or `ZNF683` when possible.

### FOXP3/CTLA4 regulatory T cell

**Core markers:** `FOXP3`, `IL2RA`, `CTLA4`, `TIGIT`, `IKZF2`, `TNFRSF18`,
`TNFRSF4`, `ENTPD1`, `CCR8`, `LAYN`, `TNFRSF9`, `RTKN2`.

**Interpretation:** Regulatory T-cell program. In tumors, an effector Treg
state is supported by `CCR8`, `TNFRSF18`, `TNFRSF4`, `TNFRSF9`, `LAYN`,
`ENTPD1`, and high checkpoint genes.

**Suggested label:** `FOXP3/CTLA4 regulatory T cell` or
`CCR8+ effector Treg`.

**Pitfalls:** Activated conventional CD4 T cells can transiently express
`FOXP3` or `IL2RA`. Stable Treg calls should combine `FOXP3`, `IL2RA`,
`CTLA4/TIGIT`, and low `IL7R` when possible.

### CXCR5/BCL6 T follicular helper T cell

**Core markers:** `CXCR5`, `BCL6`, `ICOS`, `PDCD1`, `TOX2`, `IL21`,
`CD40LG`, `SH2D1A`, `MAF`, `BTLA`.

**Interpretation:** B-cell helper or follicular helper T-cell program. This
state is most plausible in lymph node, tertiary lymphoid structure, inflamed
tissue, or tumor samples with B-cell-rich neighborhoods.

**Suggested label:** `CXCR5/BCL6 T follicular helper cell`.

**Pitfalls:** `PDCD1` is shared with exhausted T cells. Require `CXCR5`,
`BCL6`, `ICOS`, `IL21`, or `CD40LG` and usually CD4 lineage support.

### Th1, Th17, and Th2 helper T-cell programs

**Th1 markers:** `TBX21`, `IFNG`, `CXCR3`, `CCR5`, `STAT1`, `CCL5`.

**Th17 markers:** `RORC`, `IL17A`, `IL17F`, `CCR6`, `IL23R`, `KLRB1`,
`RORA`.

**Th2 markers:** `GATA3`, `IL4`, `IL5`, `IL13`, `CCR4`, `PTGDR2`.

**Interpretation:** CD4 helper-polarization programs. Cytokine transcripts
are often sparse in droplet scRNA-seq, so transcription factors and chemokine
receptors can be more useful than cytokines alone.

**Pitfalls:** `KLRB1/CCR6/RORC` can also occur in MAIT-like cells. Require CD4
or conventional helper T-cell context before calling Th17.

### ISG15/IFIT interferon-stimulated T cell

**Core markers:** `ISG15`, `IFIT1`, `IFIT2`, `IFIT3`, `IFITM1`, `IFITM3`,
`MX1`, `OAS1`, `OAS2`, `RSAD2`, `STAT1`, `IRF7`, `IFI6`.

**Interpretation:** Type-I interferon or antiviral-response T-cell state. Use
as a state modifier, for example `ISG15+ interferon-stimulated CD4 T cell`.

**Pitfalls:** IFN signatures appear in many lineages and often reflect a
sample-wide response. Confirm T-cell identity and avoid treating IFN response
as a lineage.

### MKI67/TOP2A proliferating T cell

**Core markers:** `MKI67`, `TOP2A`, `STMN1`, `TYMS`, `UBE2C`, `PCLAF`,
`CENPF`, `HMGB2`.

**Interpretation:** Cycling/proliferating T-cell state. Use as a suffix to the
best lineage or subtype call, such as `proliferating exhausted CD8 T cell` or
`proliferating Treg`.

**Pitfalls:** Cell-cycle genes are not lineage markers. Do not assign
`proliferating T cell` unless CD3/TCR genes support T-cell identity.

### HSPA/DNAJB1 stress-response T cell

**Core markers:** `HSPA1A`, `HSPA1B`, `HSPA6`, `DNAJB1`, `HSP90AA1`,
`HSPB1`, `BAG3`, `FOS`, `JUN`, `DUSP1`.

**Interpretation:** Heat-shock or stress-response T-cell state. Chu et al.
describe a pan-cancer T-cell stress-response state associated with
immunotherapy resistance signals, but the same genes can also reflect
dissociation, sample handling, or general activation stress.

**Suggested label:** `HSPA/DNAJB1 stress-response T cell`, usually as a suffix
or caution rather than a stable subtype.

**Pitfalls:** Stress genes are not T-cell specific. If this program appears in
many unrelated lineages, prefer a technical-stress caution.

### TRDC/TRGC gamma-delta T cell

**Core markers:** `TRDC`, `TRGC1`, `TRGC2`, `TRDV1`, `TRDV2`, `TRGV9`,
`TRGV10`, with pan-T markers `CD3D/E/G` and `TRAC` often lower or absent
depending on capture and annotation.

**Interpretation:** Gamma-delta T-cell lineage. In tissues and tumors, these
cells can show cytotoxic, innate-like, tissue-resident, or IL17-like states.

**Suggested label:** `TRDC+ gamma-delta T cell`.

**Pitfalls:** Use a multi-gene TCR-gamma/delta signature when possible. A
single TCR constant gene can be noisy in droplet data.

### KLRB1/SLC4A10 MAIT cell

**Core markers:** `KLRB1`, `SLC4A10`, `TRAV1-2`, `ZBTB16`, `RORA`, `RORC`,
`IL18RAP`, `CCR6`, `CXCR6`, `DPP4`, `NCR3`, with variable cytotoxic markers.

**Interpretation:** Mucosal-associated invariant T-cell program. Garner et al.
support tissue-localization, activation-state, and clonotype-linked diversity
in human MAIT cells, with core MAIT markers including `KLRB1` and `SLC4A10`.

**Suggested label:** `KLRB1/SLC4A10 MAIT cell`.

**Pitfalls:** `KLRB1` is not specific alone. Require `SLC4A10`, `TRAV1-2`, or
other MAIT-supporting genes. MAIT cells can resemble cytotoxic CD8 or Th17-like
states depending on activation and tissue.

## Practical Distinctions

- `CD3D/E/G` + `TRAC/TRBC` -> T lineage; without these, check NK, myeloid, B,
  endothelial, epithelial, and stromal contamination.
- `CCR7/SELL/TCF7/LEF1` -> naive-like or central memory; in tumors with
  `PDCD1/TOX/XCL1`, consider precursor exhaustion.
- `GZMK/CCL5/CXCR3` -> effector-memory or transitional T cell.
- `GNLY/PRF1/GZMB/NKG7` -> cytotoxic effector, but distinguish CD8 T from NK.
- `PDCD1/LAG3/HAVCR2/TOX/CXCL13/ENTPD1` -> exhausted or dysfunctional TIL.
- `FOXP3/IL2RA/CTLA4/TIGIT/CCR8` -> Treg, especially effector Treg in tumors.
- `CXCR5/BCL6/ICOS/IL21/CD40LG` -> Tfh or B-helper CD4 T cell.
- `CD69/ITGAE/CXCR6/ZNF683` -> tissue-resident memory program.
- `TRDC/TRGC1/TRGC2` -> gamma-delta T cell.
- `KLRB1/SLC4A10/TRAV1-2/ZBTB16` -> MAIT cell.
- `MKI67/TOP2A` and `HSPA1A/HSPA1B` should usually be suffixes or cautions,
  not standalone T-cell subtypes.

## Recommended Output Language

Use concise marker-program labels such as:

- CCR7/TCF7 naive-like CD4 T cell
- CCR7/TCF7 central-memory CD8 T cell
- GZMK+ effector-memory CD8 T cell
- GNLY/PRF1 cytotoxic CD8 T cell
- TCF7+ precursor-exhausted CD8 T cell
- TOX/CXCL13 exhausted CD8 T cell
- CD69/ITGAE tissue-resident memory T cell
- FOXP3/CTLA4 regulatory T cell
- CCR8+ effector Treg
- CXCR5/BCL6 T follicular helper cell
- ISG15+ interferon-stimulated T cell
- MKI67+ proliferating T cell
- HSPA/DNAJB1 stress-response T cell
- TRDC+ gamma-delta T cell
- KLRB1/SLC4A10 MAIT cell
