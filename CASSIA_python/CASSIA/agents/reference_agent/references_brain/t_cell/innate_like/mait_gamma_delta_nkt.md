---
id: innate_like_t_cell_mait_gamma_delta_nkt
category: t_cell
parent: immune_t_cell_consensus
cell_types:
  - MAIT cell
  - Mucosal-associated invariant T cell
  - Gamma-delta T cell
  - NKT-like cell
  - Innate-like T cell
condition:
  - human immune tissue
  - mucosal tissue
  - inflamed tissue
  - tumor microenvironment
  - peripheral blood
trigger_markers:
  - KLRB1
  - SLC4A10
  - TRAV1-2
  - ZBTB16
  - RORA
  - RORC
  - IL18RAP
  - CCR6
  - CXCR6
  - DPP4
  - NCR3
  - TRDC
  - TRGC1
  - TRGC2
  - TRDV1
  - TRDV2
  - TRGV9
  - TRGV10
  - CD3D
  - CD3E
  - TRAC
  - NKG7
  - GNLY
  - GZMB
  - PRF1
  - KLRD1
  - NCAM1
exclusion_markers:
  - MS4A1
  - CD79A
  - LST1
  - FCER1G
  - C1QA
  - C1QB
  - C1QC
sources:
  - Szabo et al. Nat Commun 2019. https://doi.org/10.1038/s41467-019-12464-3
  - Garner et al. Nat Immunol 2023. https://doi.org/10.1038/s41590-023-01575-1
  - Zheng et al. Science 2021. https://doi.org/10.1126/science.abe6474
  - Chu et al. Nat Med 2023. https://doi.org/10.1038/s41591-023-02371-y
---

# Innate-like T Cell Reference: MAIT, Gamma-delta, and NKT-like

## Scope

Use this document when a T-cell cluster has innate-like, unconventional TCR,
MAIT, gamma-delta, or NK-like features. The goal is to prevent MAIT,
gamma-delta T, and NKT-like cells from being overcalled as conventional CD8
cytotoxic T cells, Th17 cells, or NK cells.

## Identity Gate

**T-cell support:** `CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1`, or `TRBC2`.

**Innate-like caution:** Innate-like T cells can express cytotoxic and NK-like
genes such as `NKG7`, `GNLY`, `PRF1`, `KLRD1`, and `NCAM1`. Keep T-cell and
TCR evidence central to the annotation.

## Consensus State Catalog

### KLRB1/SLC4A10 MAIT cell

**Core markers:** `KLRB1`, `SLC4A10`, `TRAV1-2`, `ZBTB16`, `RORA`, `RORC`,
`IL18RAP`, `CCR6`, `CXCR6`, `DPP4`, `NCR3`, with variable cytotoxic markers
such as `NKG7`, `GZMK`, or `GZMB`.

**Interpretation:** Mucosal-associated invariant T-cell program. Garner et al.
supports human MAIT-cell diversity across tissues, activation states, and
clonotypes, while retaining core MAIT genes such as `KLRB1` and `SLC4A10`.

**Suggested label:** `KLRB1/SLC4A10 MAIT cell`.

**Pitfalls:** `KLRB1` alone is not specific. Require `SLC4A10`, `TRAV1-2`,
`ZBTB16`, or a broader MAIT program. MAIT cells can look Th17-like because of
`KLRB1/CCR6/RORC`, or cytotoxic because of `NKG7/GZMB`.

### TRDC/TRGC gamma-delta T cell

**Core markers:** `TRDC`, `TRGC1`, `TRGC2`, `TRDV1`, `TRDV2`, `TRGV9`,
`TRGV10`, with pan-T markers `CD3D/E/G` and variable cytotoxic or tissue
residency genes.

**Interpretation:** Gamma-delta T-cell lineage. These cells can show cytotoxic,
innate-like, tissue-resident, or IL17-like states depending on tissue and
stimulation.

**Suggested label:** `TRDC+ gamma-delta T cell`.

**Pitfalls:** A single TCR gamma/delta gene can be noisy in droplet data.
Prefer a multi-gene TCR-gamma/delta signal. Do not call conventional CD8 T
cell when `TRDC/TRGC` dominate over `CD8A/CD8B`.

### NKT-like or NK-like T cell

**Core markers:** T-cell identity genes plus NK/cytotoxic genes such as
`NKG7`, `GNLY`, `PRF1`, `GZMB`, `GZMH`, `KLRD1`, `KLRK1`, `FCGR3A`, `NCAM1`,
and `TYROBP` variably.

**Interpretation:** NKT-like or NK-like T-cell state. Use cautiously when both
T-cell receptor evidence and strong NK/cytotoxic genes are present.

**Suggested label:** `NK-like cytotoxic T cell` or `NKT-like cell`, depending
on TCR/CD3 strength and dataset context.

**Pitfalls:** If CD3/TCR genes are absent and NK genes dominate, call NK cell
instead. If `CD8B` and conventional CD8 markers dominate, call cytotoxic CD8 T
cell instead.

### Innate-like IL17 program

**Core markers:** `RORC`, `IL17A`, `IL17F`, `CCR6`, `KLRB1`, `IL23R`, `RORA`,
with possible MAIT or gamma-delta markers.

**Interpretation:** IL17-like innate T-cell program, often in MAIT or
gamma-delta cells rather than conventional CD4 Th17 cells.

**Suggested label:** `IL17-like MAIT cell` or `IL17-like gamma-delta T cell`
depending on TCR/MAIT evidence.

**Pitfalls:** Do not call Th17 unless conventional CD4 helper markers support
that interpretation.

## Practical Distinctions

- `KLRB1/SLC4A10/TRAV1-2/ZBTB16` -> MAIT.
- `TRDC/TRGC1/TRGC2/TRDV1/TRDV2` -> gamma-delta T cell.
- `NKG7/GNLY/PRF1/GZMB` plus strong TCR/CD3 but weak CD8B -> NK-like T cell or innate-like cytotoxic T cell.
- `RORC/IL17A/IL17F/CCR6` plus MAIT or gamma-delta evidence -> innate-like IL17 program, not automatically Th17.

## Recommended Output Language

- KLRB1/SLC4A10 MAIT cell
- TRDC+ gamma-delta T cell
- NK-like cytotoxic T cell
- IL17-like MAIT cell
- IL17-like gamma-delta T cell
