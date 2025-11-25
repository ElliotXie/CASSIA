---
id: t_cell_overview
category: t_cell
cell_types:
  - T cell
  - T lymphocyte
  - CD4+ T cell
  - CD8+ T cell
  - gamma-delta T cell
  - NKT cell
trigger_markers:
  - CD3D
  - CD3E
  - CD3G
  - CD2
  - CD7
  - TRAC
  - TRBC1
  - TRBC2
exclusion_markers: []
---

# T Cell Annotation Guide

## Overview

T cells (T lymphocytes) are a type of white blood cell that plays a central role in cell-mediated immunity. They are distinguished from other lymphocytes by the presence of a T-cell receptor (TCR) on their cell surface.

## Key Canonical Markers

### Pan-T Cell Markers (Required for T cell identification)
- **CD3D/CD3E/CD3G**: T cell receptor complex components - essential markers
- **CD2**: T cell surface antigen, also on NK cells
- **CD7**: Present on most T cells, but can be lost in some malignancies
- **TRAC/TRBC1/TRBC2**: T cell receptor constant regions

### Important Notes
- CD3 complex markers (CD3D, CD3E, CD3G) are the most reliable pan-T cell markers
- In scRNA-seq, CD3 detection can be variable - use multiple markers for confidence

## Major T Cell Subsets

### CD4+ T Cells (Helper T Cells)
**Key markers**: CD4+, CD3+
**Function**: Coordinate immune responses, help other immune cells
**Subtypes include**: Th1, Th2, Th17, Tfh, Treg

### CD8+ T Cells (Cytotoxic T Cells)
**Key markers**: CD8A+, CD8B+, CD3+
**Function**: Kill infected or cancerous cells
**Cytotoxic markers**: GZMA, GZMB, GZMK, PRF1, GNLY

### gamma-delta T Cells
**Key markers**: TRGC1+, TRDC+, CD3+, typically CD4- CD8-
**Function**: Bridge innate and adaptive immunity

### NKT Cells
**Key markers**: CD3+, KLRB1+, characteristics of both T and NK cells
**Function**: Rapid cytokine production

## Activation and Memory States

### Naive T Cells
**Markers**: CCR7+, SELL (CD62L)+, TCF7+, LEF1+
**Characteristics**: Have not encountered antigen, reside in lymphoid organs

### Memory T Cells
**Central memory**: CCR7+, SELL+, IL7R+
**Effector memory**: CCR7-, SELL-, higher cytotoxic capacity

### Activated/Effector T Cells
**Markers**: CD69+, HLA-DR+, CD38+, reduced CCR7/SELL
**Note**: Activation markers can vary by tissue and context

## Common Pitfalls

1. **CD4 on monocytes**: CD4 can be expressed at low levels on monocytes - always confirm T cell identity with CD3
2. **CD8 expression**: CD8A alone can appear on some NK cells and dendritic cells
3. **TCR detection in scRNA-seq**: TCR genes may have variable detection due to VDJ recombination
4. **Tissue-resident vs circulating**: Marker profiles differ between blood and tissue T cells

## Technical Notes

- In droplet-based scRNA-seq, CD4 and CD8 protein detection is often better than mRNA
- For definitive T cell identification, prioritize CD3 complex genes
- Consider the tissue context - T cell subsets vary greatly between blood, lymph nodes, and tissues
- Activation state affects many marker genes - consider the biological context
