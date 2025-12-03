---
id: cd8_t_cell
category: t_cell
parent: t_cell_overview
cell_types:
  - CD8+ T cell
  - Cytotoxic T cell
  - CTL
  - Cytotoxic T lymphocyte
trigger_markers:
  - CD8A
  - CD8B
  - GZMA
  - GZMB
  - GZMK
  - PRF1
  - NKG7
exclusion_markers:
  - CD4
---

# CD8+ T Cell Annotation Guide

## Overview

CD8+ T cells, also known as cytotoxic T lymphocytes (CTLs), are effector cells that directly kill infected or cancerous cells. They recognize antigens presented on MHC class I molecules.

## Key Canonical Markers

### Essential Markers
- **CD8A/CD8B**: Primary lineage markers (REQUIRED)
- **CD3D/CD3E/CD3G**: T cell receptor complex (confirms T cell identity)

### Cytotoxic Markers
- **GZMA**: Granzyme A - serine protease
- **GZMB**: Granzyme B - key effector molecule
- **GZMK**: Granzyme K - associated with specific subsets
- **PRF1**: Perforin - pore-forming protein
- **GNLY**: Granulysin - antimicrobial
- **NKG7**: Natural killer cell granule protein 7

## Subtype Differentiation

### Naive CD8+ T Cells
**Key markers**: CCR7+, SELL+, TCF7+, LEF1+, low cytotoxic genes
**Characteristics**: Have not encountered antigen

### Effector CD8+ T Cells
**Key markers**: GZMB high, PRF1+, IFNG+, CCR7-, SELL-
**Function**: Active killing capacity

### Memory CD8+ T Cells
**Central memory (Tcm)**: CCR7+, SELL+, IL7R+
**Effector memory (Tem)**: CCR7-, SELL-, higher cytotoxic capacity
**Tissue-resident memory (Trm)**: CD69+, ITGAE (CD103)+, tissue-specific

### Exhausted CD8+ T Cells
**Key markers**: PDCD1 (PD-1)+, CTLA4+, LAG3+, TIM3 (HAVCR2)+, TOX+
**Characteristics**: Chronic antigen exposure, reduced effector function
**Context**: Tumors, chronic infections

## Common Pitfalls

1. **CD8A on other cells**: CD8A can appear on some NK cells and dendritic cells - confirm with CD3
2. **Cytotoxic markers shared with NK**: GZMB, PRF1, NKG7 are also expressed by NK cells
3. **Exhaustion vs activation**: Some exhaustion markers overlap with activation

## Technical Notes

- CD8B is more specific to T cells than CD8A
- Use CD3 + CD8 combination for confident identification
- Cytotoxic gene expression varies with activation state
