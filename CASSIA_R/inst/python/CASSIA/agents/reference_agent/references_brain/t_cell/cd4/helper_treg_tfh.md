---
id: cd4_helper_treg_tfh
category: t_cell
parent: immune_t_cell_consensus
cell_types:
  - CD4+ T cell
  - Helper T cell
  - Regulatory T cell
  - Effector regulatory T cell
  - T follicular helper cell
  - Th1 cell
  - Th17 cell
  - Th2 cell
condition:
  - human immune tissue
  - inflamed tissue
  - tumor microenvironment
  - lymphoid tissue
trigger_markers:
  - CD4
  - IL7R
  - CCR7
  - SELL
  - TCF7
  - LEF1
  - ANXA1
  - FOXP3
  - IL2RA
  - CTLA4
  - TIGIT
  - IKZF2
  - CCR8
  - TNFRSF18
  - TNFRSF4
  - ENTPD1
  - LAYN
  - CXCR5
  - BCL6
  - ICOS
  - IL21
  - CD40LG
  - TOX2
  - TBX21
  - IFNG
  - CXCR3
  - CCR5
  - RORC
  - IL17A
  - IL17F
  - CCR6
  - IL23R
  - GATA3
  - IL4
  - IL5
  - IL13
  - CCR4
exclusion_markers:
  - CD8A
  - CD8B
  - NKG7
  - GNLY
  - PRF1
  - MS4A1
  - CD79A
  - LST1
  - FCER1G
sources:
  - Szabo et al. Nat Commun 2019. https://doi.org/10.1038/s41467-019-12464-3
  - Zheng et al. Cell 2017. https://doi.org/10.1016/j.cell.2017.05.035
  - Guo et al. Nat Med 2018. https://doi.org/10.1038/s41591-018-0045-3
  - Zheng et al. Science 2021. https://doi.org/10.1126/science.abe6474
  - Chu et al. Nat Med 2023. https://doi.org/10.1038/s41591-023-02371-y
  - Zemmour et al. Nat Immunol 2018. https://doi.org/10.1038/s41590-018-0051-0
---

# CD4 Helper, Treg, and Tfh Reference

## Scope

Use this document for CD4 T-cell subclustering. It separates naive/central
memory CD4 cells, regulatory T cells, effector Tregs, T follicular helper
cells, and polarized helper programs such as Th1, Th17, and Th2.

CD4 mRNA can be weak in droplet scRNA-seq, so combine `CD4` with T-cell
identity and supporting helper/Treg/Tfh markers.

## CD4 Identity Gate

**T-cell identity:** `CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1`, or `TRBC2`.

**CD4 support:** `CD4`, `IL7R`, `CCR7`, `LTB`, `ANXA1`, `MAL`, `TCF7`,
`LEF1`.

**Monocyte caution:** Monocytes can express `CD4`; if `LST1`, `FCER1G`,
`S100A8`, `S100A9`, `FCN1`, `VCAN`, `C1QA/B/C`, or `TYROBP` dominate, do not
call CD4 T cell without clear CD3/TCR support.

## Consensus CD4 State Catalog

### CCR7/TCF7 naive-like or central-memory CD4 T cell

**Core markers:** `CCR7`, `SELL`, `TCF7`, `LEF1`, `IL7R`, `LTB`, `MAL`,
`ANXA1`, `NOSIP`, `PIK3IP1`.

**Interpretation:** Naive-like or central-memory CD4 T-cell program. This is
common in blood and lymphoid tissues and often has low effector cytokine
expression.

**Suggested label:** `CCR7/TCF7 naive-like CD4 T cell`.

### FOXP3/CTLA4 regulatory T cell

**Core markers:** `FOXP3`, `IL2RA`, `CTLA4`, `TIGIT`, `IKZF2`, `TNFRSF18`,
`TNFRSF4`, `ENTPD1`, `RTKN2`, with low `IL7R` when measured.

**Interpretation:** Regulatory T-cell program. Zemmour et al. supports human
Treg heterogeneity, and tumor T-cell atlases commonly separate Tregs from
conventional CD4 cells by `FOXP3/IL2RA/CTLA4/TIGIT`.

**Suggested label:** `FOXP3/CTLA4 regulatory T cell`.

**Pitfalls:** Activated conventional CD4 T cells can transiently express
`FOXP3` or `IL2RA`. Stable Treg annotation should combine `FOXP3`, `IL2RA`,
`CTLA4/TIGIT`, and low `IL7R` when possible.

### CCR8/TNFRSF18 effector Treg

**Core markers:** `CCR8`, `TNFRSF18`, `TNFRSF4`, `TNFRSF9`, `LAYN`, `ENTPD1`,
`IL1R2`, `BATF`, `MAGEH1`, with Treg core markers `FOXP3`, `IL2RA`, `CTLA4`,
and `TIGIT`.

**Interpretation:** Activated or effector Treg program, especially common in
tumor and inflamed tissue contexts.

**Suggested label:** `CCR8+ effector Treg`.

**Pitfalls:** `CCR8` alone is insufficient. Require Treg core markers and CD4
T-cell identity.

### CXCR5/BCL6 T follicular helper cell

**Core markers:** `CXCR5`, `BCL6`, `ICOS`, `PDCD1`, `TOX2`, `IL21`, `CD40LG`,
`SH2D1A`, `MAF`, `BTLA`.

**Interpretation:** B-cell helper or follicular helper CD4 T-cell program.
This state is most plausible in lymph node, tertiary lymphoid structure,
inflamed tissue, or B-cell-rich tumor samples.

**Suggested label:** `CXCR5/BCL6 T follicular helper cell`.

**Pitfalls:** `PDCD1` overlaps with exhausted CD8 and activated T cells.
Require `CXCR5`, `BCL6`, `ICOS`, `IL21`, or `CD40LG`, and usually CD4 lineage
support.

### Th1 helper T-cell program

**Core markers:** `TBX21`, `IFNG`, `CXCR3`, `CCR5`, `STAT1`, `CCL5`, `TNF`.

**Suggested label:** `Th1-like CD4 T cell`.

**Pitfalls:** `IFNG` and `CCL5` also appear in cytotoxic CD8 and NK cells.
Require CD4 context and avoid Th1 calls when cytotoxic genes dominate.

### Th17 helper T-cell program

**Core markers:** `RORC`, `IL17A`, `IL17F`, `CCR6`, `IL23R`, `KLRB1`, `RORA`,
`IL22`.

**Suggested label:** `Th17-like CD4 T cell`.

**Pitfalls:** MAIT cells often express `KLRB1`, `CCR6`, and `RORC`. Check for
`SLC4A10`, `TRAV1-2`, `ZBTB16`, and MAIT context before calling Th17.

### Th2 helper T-cell program

**Core markers:** `GATA3`, `IL4`, `IL5`, `IL13`, `CCR4`, `PTGDR2`, `HPGDS`.

**Suggested label:** `Th2-like CD4 T cell`.

**Pitfalls:** Th2 cytokines are sparse in droplet scRNA-seq. Avoid confident
Th2 calls from `GATA3` alone.

## Practical Distinctions

- `CCR7/SELL/TCF7/LEF1/IL7R` -> naive-like or central-memory CD4.
- `FOXP3/IL2RA/CTLA4/TIGIT` -> Treg.
- `CCR8/TNFRSF18/TNFRSF4/LAYN/ENTPD1` plus Treg core -> effector Treg.
- `CXCR5/BCL6/ICOS/IL21/CD40LG` -> Tfh.
- `TBX21/IFNG/CXCR3` -> Th1-like, with cytotoxic/NK caution.
- `RORC/IL17A/IL17F/CCR6/IL23R` -> Th17-like, with MAIT caution.
- `GATA3/IL4/IL5/IL13/CCR4` -> Th2-like.

## Recommended Output Language

- CCR7/TCF7 naive-like CD4 T cell
- FOXP3/CTLA4 regulatory T cell
- CCR8+ effector Treg
- CXCR5/BCL6 T follicular helper cell
- Th1-like CD4 T cell
- Th17-like CD4 T cell
- Th2-like CD4 T cell
