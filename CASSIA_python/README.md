# CASSIA

**CASSIA** is a Python and R package designed for **automated, accurate, and interpretable single-cell RNA-seq cell type annotation** using a modular **multi-agent LLM framework**. CASSIA provides comprehensive annotation workflows that incorporate reasoning, validation, quality scoring, and reporting—alongside optional agents for refinement, uncertainty quantification, and retrieval-augmented generation (RAG).

## Highlights

- 🔬 **Reference-free and interpretable** LLM-based cell type annotation  
- 🧠 Multi-agent architecture with dedicated agents for annotation, validation, formatting, quality scoring, and reporting  
- 📈 **Quality scores (0–100)** and optional consensus scoring to quantify annotation reliability  
- 📊 Detailed **HTML reports** with reasoning and marker validation  
- 💬 Supports OpenAI, Anthropic, OpenRouter APIs and open-source models (e.g., LLaMA 3.2 90B)  
- 🧬 Compatible with markers from Seurat (`FindAllMarkers`) and Scanpy (`tl.rank_genes_groups`)  
- 🚀 Optional agents: Annotation Boost, Subclustering, RAG (retrieval-augmented generation), Uncertainty Quantification  
- 🌎 Cross-species annotation capabilities, validated across human, mouse, and non-model organisms  
- 🧪 Web UI also available: [https://www.cassia.bio](https://www.cassia.bio/)

## Installation

Install the core CASSIA framework:

```bash
pip install CASSIA
```

To enable optional RAG functionality:

```bash
pip install CASSIA_rag
```

**Note**: For R users, see the R package on [GitHub](https://github.com/ElliotXie/CASSIA-SingleCell-LLM-Annotation).

## Quick Start

```python
# Run the CASSIA pipeline in fast mode
CASSIA.runCASSIA_pipeline(
    output_file_name = "FastAnalysisResults",
    tissue = "large intestine",
    species = "human",
    marker = unprocessed_markers,
    max_workers = 6,  # Matches the number of clusters in dataset
    annotation_model = "openai/gpt-4o-2024-11-20", #openai/gpt-4o-2024-11-20
    annotation_provider = "openrouter",
    score_model = "anthropic/claude-3.5-sonnet",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model="anthropic/claude-3.5-sonnet",
    annotationboost_provider="openrouter"
)
```

For detailed workflows and agent customization, see the [Documentation](https://docs.cassia.bio/en/vignette/python/introduction).

## Contributing

We welcome contributions! Please submit pull requests or open issues via [GitHub](https://github.com/ElliotXie/CASSIA/issues).

## License

MIT License © 2024 Elliot Xie and contributors.

## Support

Open an issue on [GitHub](https://github.com/ElliotXie/CASSIA/issues) or visit [cassia.bio](https://www.cassia.bio/) for help.
