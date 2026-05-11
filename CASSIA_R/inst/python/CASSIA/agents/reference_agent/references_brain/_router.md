# Reference Library Router

Only select paths that are listed as available below. Prefer the most specific
subtype document over a lineage overview when the marker set supports it.

## Available Reference Files

```
references_brain/
├── b_cell/
│   └── _overview.md
├── myeloid/
│   ├── _overview.md
│   └── macrophage/
│       ├── _overview.md
│       ├── inflammatory_interferon.md
│       ├── resident_like.md
│       └── tam_pan_cancer.md
└── t_cell/
    ├── _overview.md
    ├── cd4/
    │   └── _overview.md
    └── cd8/
        └── _overview.md
```

## Routing Rules

### Macrophage / Monocyte / TAM Subclustering

Use these documents when the caller gives a macrophage/myeloid hint or the
marker set contains macrophage genes such as `CD68`, `CSF1R`, `C1QA`, `C1QB`,
`C1QC`, `APOE`, `APOC1`, `LST1`, `LYZ`, `TYROBP`, `FCER1G`, `MS4A7`, `MRC1`,
`CD163`, `FOLR2`, `SPP1`, `TREM2`, `CXCL9`, `CXCL10`, `IL1B`, or `ISG15`.

- `myeloid/macrophage/tam_pan_cancer.md`
  - Best for tumor-associated macrophage subtyping.
  - Trigger markers: `SPP1`, `AREG`, `TREM2`, `APOE`, `APOC1`, `GPNMB`,
    `LIPA`, `MMP9`, `VEGFA`, `COL1A1`, `COL1A2`, `SPARC`, `HMOX1`, `HAMP`,
    `FOLR2`, `SLC40A1`.
- `myeloid/macrophage/resident_like.md`
  - Best for resident-like macrophages, tissue macrophages, phagocytic
    complement macrophages, and monocyte-vs-macrophage separation.
  - Trigger markers: `C1QA`, `C1QB`, `C1QC`, `APOE`, `APOC1`, `RNASE1`,
    `FOLR2`, `SELENOP`, `SLC40A1`, `STAB1`, `LYVE1`, `MRC1`, `CD163`,
    `MARCO`, `FABP4`, `PPARG`, `TMEM119`, `P2RY12`, `CLEC4F`, `VSIG4`.
- `myeloid/macrophage/inflammatory_interferon.md`
  - Best for inflammatory, interferon-stimulated, inflammasome-like, and
    heat-stress macrophage states.
  - Trigger markers: `CXCL9`, `CXCL10`, `GBP1`, `GBP5`, `STAT1`, `ISG15`,
    `IFIT1`, `IFITM1`, `IFITM3`, `MX1`, `IL1B`, `TNF`, `CXCL8`, `CCL3`,
    `CCL4`, `NLRP3`, `HSPA1A`, `HSPA1B`, `HSPA6`, `DNAJB1`.
- `myeloid/macrophage/_overview.md`
  - Use with one or more specific macrophage files, or alone when the marker
    set is myeloid/macrophage but subtype evidence is weak.

### Broad Myeloid

- `myeloid/_overview.md`
  - Use for generic myeloid/monocyte/DC distinction when macrophage subtype
    markers are not dominant.

### Lymphoid Overviews

- `t_cell/_overview.md`, `t_cell/cd4/_overview.md`, `t_cell/cd8/_overview.md`
  - Use for T cell markers such as `CD3D`, `CD3E`, `TRAC`, `CD4`, `CD8A`,
    `CD8B`, `IL7R`, `CCR7`, `FOXP3`, `GZMB`, `PRF1`, `NKG7`, `TOX`,
    `PDCD1`, `LAG3`, `HAVCR2`.
- `b_cell/_overview.md`
  - Use for B/plasma cell markers such as `MS4A1`, `CD79A`, `CD79B`, `CD19`,
    `CD27`, `MZB1`, `JCHAIN`, `XBP1`, `SDC1`, immunoglobulin heavy or light
    chain genes.

## Important Selection Guidance

- If the library does not contain a good match, return an empty
  `selected_references` list.
- For macrophage subclustering, select 1-3 macrophage files. Use
  `myeloid/macrophage/_overview.md` plus the most specific subtype file when
  there is enough context budget.
- Do not select planned but unavailable paths.
