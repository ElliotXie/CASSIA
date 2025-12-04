---
title: Extended Analysis
---

This tutorial covers advanced features of CASSIA, including uncertainty quantification, annotation boost, and specialized agents.

## Optional Step: Uncertainty Quantification

This could be useful to study the uncertainty of the annotation, and potentially improve the accurracy.

Note: This is step could be costy, since multiple iteration will be performed.

```python
output_name="intestine_detailed"

# Run multiple iterations
iteration_results = CASSIA.runCASSIA_batch_n_times(
    n=2,
    marker=unprocessed_markers,
    output_name=output_name + "_Uncertainty",
    model="anthropic/claude-sonnet-4.5",
    provider="openrouter",
    tissue="large intestine",
    species="human",
    max_workers=6,
    batch_max_workers=3  # Conservative setting for API rate limits
)


# Calculate similarity scores
similarity_scores = CASSIA.runCASSIA_similarity_score_batch(
    marker=unprocessed_markers,
    file_pattern=output_name + "_Uncertainty_*_full.csv",
    output_name="intestine_uncertainty",
    max_workers=6,
    model="openai/gpt-5.1",
    provider="openrouter",
    main_weight=0.5,
    sub_weight=0.5
)
```

## Optional Step: Annotation Boost on Selected Cluster

The monocyte cluster is sometimes annotated as mixed population of immune cell and neuron/glia cells.

Here we use annotation boost agent to test these hypotheses in more detail.

```python
# Run validation plus for the high mitochondrial content cluster
CASSIA.runCASSIA_annotationboost(
    full_result_path = output_name + "_full.csv",
    marker = unprocessed_markers,
    output_name = "monocyte_annotationboost",
    cluster_name = "monocyte",
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

## Optional Step: Retrieve Augmented Generation (RAG)

This is particularly useful if you have a very specific and detialed annottaion to work with. It can significantly imrpove the granularity and accuracy of the annotation. It automatically extract marker information and genearte a report as additional informatyion for default CASSIA pipeline.

```bash
pip install cassia-rag
```

```python
from cassia_rag import run_complete_analysis
import os

os.environ["ANTHROPIC_API_KEY"] = "your-anthropic-key"
os.environ["OPENAI_API_KEY"] = "your-openai-key"

run_complete_analysis(
        tissue_type="Liver", # tissue you are analyzing
        target_species="Tiger", # species you are analyzing
        reference_species="Human", # either Human or mouse, if other species, then use Human instead of mouse
        model_choice='claude', # either claude or gpt, highly recommend claude
        compare=True,  # if you want to compare with reference species, for example fetal vs human, then set to True
        db_path="~/Canonical_Marker (1).csv", # path to the database
        max_workers=8
)
```

## Optional Step: Compare the Subtypes Using Multiple LLMs (Symphony Compare)

This agent can be used after you finish the default CASSIA pipeline, and are still unsure about a celltype. You can use this agent to get a more confident subtype annotation. Here we use the Plasma Cells cluster as examples. To distinguish if it is more like a general plasma cell or other celltypes.

```python
# The markers here are copied from CASSIA's previous results.
marker = "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

# Run Symphony Compare
results = CASSIA.symphonyCompare(
    tissue = "large intestine",
    celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"],
    marker_set = marker,
    species = "human",
    model_preset = "premium",
    output_basename = "plasma_cell_subtype",
    enable_discussion = True
)

print(f"Consensus: {results['consensus']} (confidence: {results['confidence']:.1%})")
```

## Optional Step: Subclustering

This agent can be used to study subclustered population, such as a T cell population or a Fibroblast cluster. We recommend to apply the default cassia first, and on a target cluster, apply Seurat pipeline to subcluster the cluster and get the findallmarke results to be used here. Here we present the results for the cd8-positive, alpha-beta t cell cluster as example. This cluster is a cd8 population mixed with other celltypes.

```python
CASSIA.runCASSIA_subclusters(marker = subcluster_results,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter")
```

It is recommended to run the CS score for the subclustering to get a more confident answer.

```python
CASSIA.runCASSIA_n_subcluster(
    n=5, 
    marker=subcluster_results,
    major_cluster_info="cd8 t cell", 
    base_output_name="subclustering_results_n",
    model="anthropic/claude-sonnet-4.5",
    temperature=0,
    provider="openrouter",
    max_workers=5,
    n_genes=50
)

# Calculate similarity scores
CASSIA.runCASSIA_similarity_score_batch(
    marker = subcluster_results,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

