---
id: cd4_t_cell
category: t_cell
parent: t_cell_overview
cell_types:
  - CD4+ T cell
  - Helper T cell
  - Th1
  - Th2
  - Th17
  - Tfh
  - Treg
  - Regulatory T cell
trigger_markers:
  - CD4
  - IL7R
  - CCR7
  - SELL
  - TCF7
  - LEF1
  - FOXP3
  - IL2RA
exclusion_markers:
  - CD8A
  - CD8B
  - NCAM1
  - NKG7
---

# CD4+ T Cell Annotation Guide

## Overview

CD4+ T cells, also known as helper T cells, are a subset of T lymphocytes that play a central role in adaptive immunity. They help coordinate immune responses by activating other immune cells through cytokine secretion.

## Key Canonical Markers

### Essential Markers
- **CD4**: Primary lineage marker (REQUIRED for CD4+ T cell identification)
- **CD3D/CD3E/CD3G**: T cell receptor complex (confirms T cell identity)

### Supporting Markers
- **IL7R (CD127)**: Common on most CD4+ T cells, lower in Tregs
- **CCR7**: Homing receptor, high in naive and central memory
- **SELL (CD62L)**: Lymph node homing, high in naive cells

## Subtype Differentiation

### Naive CD4+ T Cells
**Key markers**: CCR7+, SELL+, TCF7+, LEF1+, CD45RA+
**Characteristics**:
- Have not encountered antigen
- High expression of homing receptors
- Low cytokine production

### Th1 Cells
**Key markers**: TBX21 (T-bet)+, IFNG+, CXCR3+, CCR5+
**Transcription factor**: T-bet (TBX21)
**Cytokines**: IFN-gamma, TNF-alpha, IL-2
**Function**: Anti-intracellular pathogen immunity

### Th2 Cells
**Key markers**: GATA3+, IL4+, IL5+, IL13+, CCR4+
**Transcription factor**: GATA3
**Cytokines**: IL-4, IL-5, IL-13
**Function**: Anti-helminth immunity, allergy

### Th17 Cells
**Key markers**: RORC (RORgt)+, IL17A+, IL17F+, CCR6+, IL23R+
**Transcription factor**: RORgt (RORC)
**Cytokines**: IL-17A, IL-17F, IL-22
**Function**: Mucosal immunity, anti-fungal

### Follicular Helper T Cells (Tfh)
**Key markers**: CXCR5+, BCL6+, PDCD1 (PD-1)+, ICOS+, IL21+
**Transcription factor**: BCL6
**Location**: Germinal centers
**Function**: B cell help, antibody production

### Regulatory T Cells (Tregs)
**Key markers**: FOXP3+ (DEFINITIVE), IL2RA (CD25)high, CTLA4+, TIGIT+
**Transcription factor**: FOXP3 (gold standard)
**Important markers**:
- FOXP3: The defining marker - REQUIRED for Treg identification
- IL2RA (CD25): High expression, but not specific
- CTLA4: Constitutively expressed
- IL7R (CD127): LOW expression (unlike other CD4+ T cells)

**Treg subtypes**:
- Natural Tregs (nTreg): Thymus-derived, HELIOS+
- Induced Tregs (iTreg): Peripherally-derived

## Common Pitfalls

1. **CD4 expression on monocytes**: CD4 can be expressed at low levels on monocytes - always confirm with CD3
2. **FOXP3 transient expression**: Activated conventional CD4+ T cells can transiently express FOXP3 - look for stable high expression
3. **IL7R in Tregs**: Tregs have LOW IL7R - don't exclude Tregs based on IL7R positivity
4. **Th subset plasticity**: T helper subsets can be plastic and convert between phenotypes
5. **Mixed signatures**: Cells may show mixed Th1/Th2 or other combined signatures

## Technical Notes

- CD4 mRNA detection in scRNA-seq can be variable due to low expression
- Use protein-level data (CITE-seq) when available for CD4 confirmation
- Transcription factors (TBX21, GATA3, RORC, FOXP3) are more reliable than cytokines for subset identification
- Consider tissue context - Th subsets vary by anatomical location
