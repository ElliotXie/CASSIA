---
title: RAG (Optional)
---

This is particularly useful if you have a very specific and detailed annotation to work with. It can significantly improve the granularity and accuracy of the annotation. It automatically extracts marker information and generates a report as additional information for the default CASSIA pipeline.

### Canonical Marker Database

The RAG agent uses a canonical marker database to provide curated, literature-backed marker genes for cell type identification. This database contains established markers for various cell types across different tissues and species, compiled from published research.

**Download the Canonical Marker Database:** [Google Drive Link](https://drive.google.com/drive/folders/1LHha0BL4LCh4noQXvMAyvIjYVJ14QaMp?usp=drive_link)

The database file (`Canonical_Marker.csv`) should be downloaded and the path provided to the `db_path` parameter.

### Installation

```bash
pip install cassia-rag
```

### Usage

```python
from cassia_rag import run_complete_analysis
import os

# Set up the API keys if you have not done so.
os.environ["ANTHROPIC_API_KEY"] = "your-anthropic-key"
os.environ["OPENAI_API_KEY"] = "your-openai-key"

# Run the wrapper function to trigger a multiagent pipeline.
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

### Parameters

| Parameter | Description |
|-----------|-------------|
| `tissue_type` | The tissue type you are analyzing (e.g., "Liver", "Brain", "Heart") |
| `target_species` | The species of your dataset (e.g., "Tiger", "Mouse", "Human") |
| `reference_species` | The reference species for marker comparison. Use "Human" or "Mouse". For other species, use "Human" |
| `model_choice` | LLM to use: `'claude'` or `'gpt'`. Claude is highly recommended |
| `compare` | When set to `True`, enables cross-species comparison mode. This is useful when analyzing developmental stages (e.g., fetal vs adult) or comparing cell types between your target species and the reference species. The comparison identifies species-specific markers and highlights differences in cell type definitions |
| `db_path` | Path to the Canonical Marker database CSV file |
| `max_workers` | Number of parallel workers for processing |

### Output Files

All the outputs (intermediate and final) are saved as a txt file in the "TissueType_Species" folder. In our example, it is Liver_Tiger folder. 
Final output is in `summary_clean.txt` file. And the content in this file can be used as additional information in CASSIA pipeline later.

There are also some other files in the folder, which are intermediate outputs. 
Use the tutorial input as example, the files are:
1. `liver_tiger_marker_analysis.txt` # marker analysis and interpretation from the database
2. `final_ontology.txt` # ontology related to the tissue type and target species
3. `cell_type_patterns_claude.txt` # cell type patterns analysis from the
4. `summary.txt` # raw summary file
5. `additional_considerations.txt` # additional considerations if we have different species than reference species.

