---
title: RAG (Optional)
---

This is particularly useful if you have a very specific and detialed annottaion to work with. It can significantly imrpove the granularity and accuracy of the annotation. It automatically extract marker information and genearte a report as additional informatyion for default CASSIA pipeline.

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

