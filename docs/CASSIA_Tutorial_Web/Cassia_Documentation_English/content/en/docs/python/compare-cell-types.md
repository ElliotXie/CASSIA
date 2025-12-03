---
title: Symphony Compare (Optional)
---

The `symphonyCompare` agent acts as a virtual panel of experts to resolve ambiguous cell type annotations. It orchestrates multiple AI models (a "Symphony" of agents) to compare potential cell types, debate their findings in discussion rounds, and reach a consensus based on marker gene evidence.

### Usage

This agent is particularly useful after running the default CASSIA pipeline if you are unsure about a specific cluster's identity. For example, distinguishing between different subtypes of Plasma Cells.

```python
# The marker here are copy from CASSIA's previous results.
marker = "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

# Run the Symphony Compare analysis
results = CASSIA.symphonyCompare(
    tissue = "large intestine",
    celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"],
    marker_set = marker,
    species = "human",
    model_preset = "premium",  # Options: "premium", "budget"
    output_basename = "plasma_cell_comparison",
    enable_discussion = True
)

print(f"Consensus: {results['consensus']} (confidence: {results['confidence']:.1%})")
```

### Parameter Details

- **`tissue`**: The tissue type being analyzed (e.g., "large intestine").
- **`celltypes`**: A list of 2-4 cell types to compare.
- **`marker_set`**: A string of comma-separated marker genes.
- **`species`**: The species of the sample (default: "human").
- **`model_preset`**: Configuration of models to use.
    - `"premium"` (Default): High-performance ensemble (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4).
    - `"budget"`: Cost-effective models (DeepSeek V3.2, Grok 4 Fast, Kimi K2, Gemini 2.5 Flash).
- **`output_basename`**: Base name for output files.
- **`enable_discussion`**: Whether to enable multi-round debate between models (default: `True`).
- **`max_discussion_rounds`**: Maximum number of discussion rounds (default: 2).
- **`consensus_threshold`**: Fraction of models required for consensus (default: 0.8).

### Output Files
- `{output_basename}.csv`: Detailed comparison results, reasoning, and scores from all models and rounds.
- `{output_basename}_report.html`: An interactive HTML report visualizing the debate and consensus process.
