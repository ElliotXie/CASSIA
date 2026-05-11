---
id: macrophage_human_cancer_consensus
category: myeloid
cell_types:
  - Tumor-associated macrophage
  - TAM
  - Human cancer macrophage
  - Macrophage
condition:
  - human cancer
  - solid tumor
  - tumor microenvironment
trigger_markers:
  - SPP1
  - AREG
  - TREM2
  - APOE
  - APOC1
  - C1QA
  - C1QB
  - C1QC
  - C3
  - CXCL9
  - CXCL10
  - ISG15
  - IFIT1
  - IFI27
  - IL1B
  - TNF
  - FOLR2
  - SELENOP
  - SLC40A1
  - HMOX1
  - MT1G
  - COL1A1
  - COL1A2
  - SPARC
exclusion_markers:
  - EPCAM
  - KRT8
  - KRT18
  - KRT19
  - PECAM1
sources:
  - Cheng et al. Cell 2021. https://doi.org/10.1016/j.cell.2021.01.010
  - Mulder et al. Immunity 2021. https://doi.org/10.1016/j.immuni.2021.07.007
  - Coulton et al. Nat Commun 2024. https://doi.org/10.1038/s41467-024-49885-8
  - Li et al. Exp Mol Med 2023. https://doi.org/10.1038/s12276-023-01115-9
  - Li et al. Nat Commun 2024. https://doi.org/10.1038/s41467-024-50478-8
---

# Human Cancer Macrophage Consensus Reference

## Scope

This document is a consensus layer for human cancer macrophage/TAM
subclustering. It is not a universal macrophage ontology and it is not a table
of one paper's cluster names. Use the consensus program as the primary subtype
label, then mention paper-specific labels as traceability evidence when useful.

Good output pattern:

`APOE/TREM2 lipid-phagolysosomal TAM (Coulton 4_ICIMac2-like)`

Bad output pattern:

`ICIMac2` with no explanation of the broader marker program or paper source.

## Evidence Model

- **Consensus name**: preferred subtype wording for CASSIA output.
- **Core markers**: marker program that should drive the call.
- **Paper aliases**: paper-specific labels that support traceability but should
  not be presented as universal cell ontology.
- **Confidence**: strength of cross-paper support for using this as a human
  cancer macrophage program.
- **Pitfalls**: common overcalls and contaminating cell programs.

Cheng et al. Cell 2021 supports broad pan-cancer TIM/TAM states such as SPP1+,
C1QC+, ISG15+, FN1/ECM-like, LYVE1+, NLRP3+, and INHBA+ macrophage programs.
Mulder et al. Immunity 2021 supports conserved cross-tissue macrophage and
resident-like marker programs. Coulton et al. Nat Commun 2024 provides a large
pan-cancer TAM atlas across many studies and cancer types, useful for mapping
consensus programs to author-defined cluster aliases. Li et al. Exp Mol Med
2023 provides an additional tumor-specific macrophage example in uveal
melanoma. Li et al. Nat Commun 2024 provides an independent pan-cancer ICB
myeloid atlas with author-defined Macro_FOLR2/APOE, Macro_ISG15,
Macro_NLRP3, Macro_OLFML3, and Macro_IFI27 states.

## Consensus State Catalog

### SPP1/AREG inflammatory angiogenic TAM

**Confidence:** high

**Core markers:** `SPP1`, `AREG`, `EREG`, `CXCL3`, `CXCL2`, `CXCL8`, `CXCL1`,
`CCL20`, `IL1B`, `TIMP1`, `SOD2`, `IL1RN`, `MMP9`, `VEGFA`.

**Interpretation:** Tumor-conditioned inflammatory, angiogenic, epithelial
growth factor ligand, and matrix-remodeling macrophage program. This is a
better label than generic M2 or generic TAM when SPP1/AREG/EREG and neutrophil
chemokines co-occur.

**Paper aliases:** Cheng SPP1+ angiogenesis-associated TAM; Coulton
`6_SPP1AREGMac`.

**Pitfalls:** `SPP1` alone is not enough. Require macrophage context and a
co-program such as `AREG/EREG`, `CXCL1/2/3/8`, `TIMP1`, or `MMP9`. If
`FCN1/S100A8/S100A9/VCAN` dominate, consider monocyte-like inflammatory
macrophage instead.

### APOE/TREM2 lipid-phagolysosomal TAM

**Confidence:** high

**Core markers:** `TREM2`, `APOE`, `APOC1`, `GPNMB`, `LIPA`, `CTSD`, `LGMN`,
`PLA2G7`, `ACP5`, `PLD3`, `PSAP`.

**Interpretation:** Lipid-handling, phagolysosomal, often immunosuppressive
tumor macrophage program. Use this label when TREM2/APOE/APOC1 and lysosomal
lipid genes dominate.

**Paper aliases:** TREM2+ lipid-associated TAM; Coulton `4_ICIMac2`.

**Pitfalls:** Do not collapse this into resident-like FOLR2 TAM if
`TREM2/GPNMB/LIPA/CTSD/PLA2G7` dominate. Do not call it foam cell without tumor
or macrophage context.

### AREG/THBS1 angiogenic remodeling TAM

**Confidence:** medium-high

**Core markers:** `AREG`, `THBS1`, `EREG`, `NAMPT`, `ZNF331`, `IL1B`,
`GPR183`, `NR4A3`, `G0S2`, `BTG1`.

**Interpretation:** EGF-ligand, angiogenic, and tissue-remodeling TAM program.
This is related to but separable from the more chemokine-dominant SPP1/AREG
inflammatory angiogenic program. Use this label when `AREG/EREG/THBS1/NAMPT`
dominate and `SPP1` or the full CXCL1/2/3/8 inflammatory chemokine block is
weak.

**Paper aliases:** EGF-ligand angiogenic TAM; Coulton `9_AngioMac`.

**Pitfalls:** `AREG` and `EREG` overlap with SPP1/AREG inflammatory TAM. Keep
the angiogenic remodeling label when `THBS1`, `NAMPT`, `GPR183`, and `NR4A3`
support vascular or tissue-remodeling biology.

### FOLR2/SELENOP resident-like iron-handling TAM

**Confidence:** high

**Core markers:** `FOLR2`, `SELENOP`, `SLC40A1`, `F13A1`, `STAB1`, `RNASE1`,
`CD163`, `MRC1`, `LYVE1`, `DAB2`, `PLTP`, `LGMN`.

**Interpretation:** Tissue-adapted, resident-like, iron-handling macrophage
program within tumors. It is often better described as resident-like or
metabolic TAM than as simple M2.

**Paper aliases:** FOLR2+ resident-like TAM; Coulton `1_MetM2Mac`; related to
Mulder conserved macrophage/resident programs.

**Pitfalls:** If `HMOX1/CD163/SLC40A1` dominate, consider heme-iron TAM. If
`C1QA/B/C` and MHC genes dominate, consider complement antigen-presenting TAM.

### C1QC/C3 complement antigen-presenting TAM

**Confidence:** high

**Core markers:** `C1QA`, `C1QB`, `C1QC`, `C3`, `APOE`, `APOC1`, `CD74`,
`HLA-DRA`, `HLA-DPA1`, `HLA-DPB1`, `CX3CR1`, `PLD4`, `GPR34`.

**Interpretation:** Complement-rich, phagocytic, antigen-presenting macrophage
state. This program often separates from SPP1/angiogenic TAM programs in tumor
atlases.

**Paper aliases:** Cheng C1QC+ phagocytic TAM; Coulton `2_C3Mac`.

**Pitfalls:** MHC genes can suggest dendritic cells, but macrophage markers
such as `C1QA/B/C`, `APOE`, `CX3CR1`, and `PLD4` support macrophage identity.

### CXCL9/CXCL10 IFNG-response TAM

**Confidence:** high

**Core markers:** `CXCL9`, `CXCL10`, `GBP1`, `GBP5`, `STAT1`, `WARS1`,
`SLAMF7`, `LGALS2`, `VAMP5`, `CALHM6`, `MMP9`.

**Interpretation:** IFN-gamma response macrophage/TAM, often linked to T cell
interaction, antigen presentation, and immune-active tumor regions.

**Paper aliases:** IFNGMac; Coulton `8_IFNGMac`.

**Pitfalls:** Separate from type-I IFN TAM by `GBP1/GBP5/STAT1/CXCL9`
dominance rather than `ISG15/IFIT/MX1/RSAD2` dominance.

### ISG15/IFIT type-I interferon TAM

**Confidence:** high

**Core markers:** `ISG15`, `IFIT1`, `IFIT2`, `IFIT3`, `IFITM1`, `IFITM2`,
`IFITM3`, `MX1`, `OAS1`, `OAS2`, `RSAD2`, `CXCL10`, `TNFSF10`.

**Interpretation:** Type-I interferon or viral-mimicry macrophage program.
This can reflect tumor nucleic acid sensing, cGAS-STING/TLR stimulation, or
therapy-induced interferon signaling.

**Paper aliases:** Cheng ISG15+ TAM; Coulton `17_IFNMac3` or other IFNMac
aliases depending on marker details.

**Pitfalls:** IFN markers occur in many immune cells. Require macrophage or
parent-cluster context. Do not label this IFNG-response if `CXCL9/GBP/STAT1`
are weak.

### IFI27/APOE/C1Q interferon-lipid TAM

**Confidence:** medium

**Core markers:** `IFI27`, `APOE`, `APOC1`, `C1QA`, `C1QB`, `C1QC`, `GPNMB`,
`TREM2`, `A2M`, `FTL`, `NUPR1`.

**Interpretation:** Hybrid macrophage state combining lipid/phagolysosomal and
complement-rich TAM features with a focused interferon-exposed signal. Use this
label when `IFI27` is prominent but the broader `ISG15/IFIT/MX1/RSAD2` block is
not dominant, and `APOE/C1Q/GPNMB/TREM2` still anchor macrophage identity.

**Paper aliases:** Li et al. Nat Commun 2024 `Macro_IFI27`.

**Pitfalls:** Do not force this into pure APOE/TREM2 lipid TAM if `IFI27` is a
top marker. Do not force it into ISG15/IFIT type-I IFN TAM unless multiple
canonical ISGs (`ISG15`, `IFIT1/2/3`, `MX1`, `RSAD2`) dominate. `C1QA/B/C`
overlap with complement antigen-presenting TAM; keep the interferon-lipid label
when `IFI27` and APOE/APOC1/GPNMB are both present.

### IL1B/TNF inflammatory TAM

**Confidence:** high

**Core markers:** `IL1B`, `TNF`, `CXCL8`, `CXCL1`, `CXCL2`, `CXCL3`, `CCL3`,
`CCL4`, `CCL3L3`, `CCL4L2`, `DUSP2`, `IER3`, `BCL2A1`, `NFKBIA`.

**Interpretation:** Acute inflammatory, NF-kB-active macrophage/monocyte-
derived macrophage state with cytokine and chemokine output.

**Paper aliases:** inflammatory TAM; Coulton `10_InflamMac`.

**Pitfalls:** Overlaps with SPP1/AREG inflammatory angiogenic TAM. `TNF` plus
`CCL3/CCL4` supports InflamMac-like inflammatory TAM; `AREG/EREG/TIMP1` and
SPP1-like markers support SPP1/AREG inflammatory angiogenic TAM.

### HMOX1/SLC40A1 heme-iron handling TAM

**Confidence:** medium-high

**Core markers:** `HMOX1`, `SLC40A1`, `HAMP`, `CD163`, `CCL18`, `LGMN`,
`STAB1`, `SELENOP`, `CTSB`, `CTSL`.

**Interpretation:** Heme-processing, erythrophagocytic, iron-exporting
macrophage program, often plausible in hemorrhagic or necrotic tumor regions.

**Paper aliases:** HemeMac; Coulton `21_HemeMac`.

**Pitfalls:** Distinguish from FOLR2/SELENOP resident-like TAM when `HMOX1`
and `CD163` are weak. Distinguish from APOE/TREM2 lipid TAM when TREM2/GPNMB
dominate.

### MT1/MT2 metallothionein stress TAM

**Confidence:** medium-high

**Core markers:** `MT1G`, `MT1X`, `MT2A`, `MT1E`, `MT1H`, `MT1F`, `MT1M`,
`MIF`, `SPP1`, `LDHA`, `LGALS1`.

**Interpretation:** Metal-ion, oxidative-stress, or metabolic-stress
macrophage/TAM program. It should be labeled by metallothionein dominance
rather than as generic SPP1 TAM.

**Paper aliases:** MetalloMac; Coulton `11_MetalloMac`.

**Pitfalls:** Metallothioneins can be stress responses in many cells. Require
macrophage context or parent macrophage clustering.

### HSPA heat-shock stress macrophage

**Confidence:** medium

**Core markers:** `HSPA6`, `HSPA1B`, `HSPA1A`, `DNAJB1`, `HSPB1`, `BAG3`,
`HSPH1`, `HSP90AA1`, `ZFAND2A`, `HSPD1`, `HSPE1`, `IER5`.

**Interpretation:** Heat-shock or proteotoxic stress macrophage state. This
can be biological stress, therapy pressure, hypoxia, or tissue dissociation.

**Paper aliases:** StressMac; Coulton `5_StressMac`.

**Pitfalls:** Treat as a state label, not a stable lineage. If heat-shock
genes appear across many cell types, consider technical stress.

### MKI67/TOP2A proliferating macrophage

**Confidence:** medium

**Core markers:** `MKI67`, `TOP2A`, `STMN1`, `PCLAF`, `UBE2C`, `TYMS`,
`CENPF`, `HMGB2`, `TUBB`, `H2AZ1`.

**Interpretation:** Cycling/proliferating macrophage state. Use this as a
cell-cycle state layered on macrophage identity, not as a separate lineage.

**Paper aliases:** ProliMac; Coulton `14_ProliMac`.

**Pitfalls:** Require macrophage identity markers elsewhere. Cell-cycle genes
alone cannot establish macrophage identity.

### FCN1/S100A8 monocyte-like inflammatory macrophage

**Confidence:** medium

**Core markers:** `FCN1`, `S100A8`, `S100A9`, `S100A12`, `VCAN`, `LYZ`,
`LST1`, `CCR2`, `IL1B`, `EREG`, `RETN`, `PLAC8`.

**Interpretation:** Recruited monocyte-like or early monocyte-derived
macrophage state in tumor. Often inflammatory and less resident-like.

**Paper aliases:** classical monocyte-like TAM; Coulton `19_ClassMono`.

**Pitfalls:** If `FCN1/S100A8/S100A9` dominate without macrophage maturation
markers, call monocyte-like rather than mature macrophage.

### COL1A1/SPARC ECM-remodeling macrophage-like state

**Confidence:** ambiguous

**Core markers:** `COL1A1`, `COL1A2`, `COL3A1`, `SPARC`, `LUM`, `DCN`,
`POSTN`, `BGN`, `MGP`, `SFRP2`, `CALD1`.

**Interpretation:** ECM-remodeling macrophage-like or macrophage-associated
matrix program in tumor data. This is paper-supported but biologically
ambiguous because the marker set is also highly compatible with fibroblasts.

**Paper aliases:** ECMMac; METAM-like; Coulton `18_ECMMac`.

**Pitfalls:** If macrophage identity markers (`PTPRC`, `LST1`, `TYROBP`,
`C1QA/B/C`, `APOE`, `CD68`) are absent, flag fibroblast/CAF contamination or
doublet as an important alternative. Do not force macrophage identity from ECM
genes alone.

## Tissue Modifiers

Do not split the first-pass consensus by tissue. Use tissue as a modifier:

- Lung or lung metastasis: `FABP4`, `MARCO`, `PPARG`, `MCEMP1` support
  alveolar-like macrophage features.
- Liver or liver metastasis: `CLEC4F`, `VSIG4`, `MARCO`, `TIMD4`, `CD5L`
  support Kupffer-like features.
- Brain tumor: `TMEM119`, `P2RY12`, `SALL1`, `CX3CR1`, `CSF1R` support
  microglia-like or glioma-associated macrophage context.
- Bone or giant-cell lesion: `CTSK`, `ACP5`, `MMP9`, `SIGLEC15` support
  osteoclast-like macrophage.

## Output Rules

- Prefer consensus labels over M1/M2 labels.
- Add a paper alias only after the consensus label, e.g. "Coulton
  8_IFNGMac-like".
- State uncertainty explicitly when the marker program is tissue-specific,
  technical-stress-like, or contamination-prone.
- For ECM and heat-shock programs, mention the main biological alternative
  because marker-only annotation can be misleading.
