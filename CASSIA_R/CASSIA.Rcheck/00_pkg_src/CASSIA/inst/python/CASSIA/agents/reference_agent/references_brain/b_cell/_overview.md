---
id: b_cell_overview
category: b_cell
cell_types:
  - B cell
  - B lymphocyte
  - Plasma cell
  - Plasmablast
  - Memory B cell
  - Naive B cell
trigger_markers:
  - CD19
  - CD79A
  - CD79B
  - MS4A1
  - PAX5
  - CD27
  - IGHM
  - IGHD
exclusion_markers:
  - CD3D
  - CD3E
---

# B Cell Annotation Guide

## Overview

B cells (B lymphocytes) are responsible for humoral immunity through antibody production. They develop in the bone marrow and differentiate into antibody-secreting plasma cells upon activation.

## Key Canonical Markers

### Pan-B Cell Markers
- **CD19**: Essential B cell marker
- **CD79A/CD79B**: B cell receptor signaling components
- **MS4A1 (CD20)**: B cell surface marker
- **PAX5**: B cell transcription factor

### Stage-Specific Markers
- **IGHM/IGHD**: Naive B cells (surface IgM/IgD)
- **CD27**: Memory B cells
- **SDC1 (CD138)**: Plasma cells
- **PRDM1 (BLIMP1)**: Plasma cell transcription factor

## Subtype Differentiation

### Naive B Cells
**Key markers**: IGHD+, IGHM+, CD27-, MS4A1+
**Characteristics**: Have not encountered antigen

### Memory B Cells
**Key markers**: CD27+, class-switched (IgG, IgA)
**Subtypes**: IgM memory, switched memory

### Plasma Cells
**Key markers**: SDC1 (CD138)+, PRDM1+, XBP1+, MZB1+
**Characteristics**: Low/absent CD19 and MS4A1, high Ig secretion

### Plasmablasts
**Key markers**: Intermediate between activated B and plasma cells
**Characteristics**: Actively proliferating, antibody-secreting

## Common Pitfalls

1. **Plasma cells lose B cell markers**: CD19 and MS4A1 are often low/absent in plasma cells
2. **CD27 on T cells**: CD27 is also expressed on some T cells
3. **Detection of secreted Ig**: Plasma cells have very high Ig mRNA that can contaminate other cells

## Technical Notes

- Use multiple B cell markers (CD19 + CD79A/B) for confident identification
- Plasma cells require different markers (SDC1, PRDM1) than other B cells
- BCR sequences (IGHV genes) can help identify B cell lineage

*This is a placeholder reference. Add detailed content as needed.*
