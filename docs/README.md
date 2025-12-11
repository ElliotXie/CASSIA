<div align="center">

<img src="logo2.png" width="200" style="vertical-align: middle;" />

[English](README.md) | [中文](README_CN.md)

</div>

**CASSIA** (Collaborative Agent System for Single-cell Interpretable Annotation) is a tool that enhances cell type annotation using multi-agent Large Language Models (LLMs).


🌐 [CASSIA Web UI (cassia.bio)](https://cassia.bio/) - Try CASSIA's core features online. For a comprehensive experience with all advanced features, use our R or Python package.

📚 [Complete Documentation/Vignette (docs.cassia.bio)](https://docs.cassia.bio/en)

🤖 [LLMs Annotation Benchmark (sc-llm-benchmark.com)](https://sc-llm-benchmark.com/methods/cassia)



## 📰 News

> **2025-11-29**
>🎇 **Major update with new features and improvements!**
> - **Python Documentation**: Complete Python docs and vignettes now available
> - **Annotation Boost Improvements**: Sidebar navigation, better reports, bug fixes
> - **Better Scanpy Support**: Fixed marker processing, improved R/Python sync
> - **Symphony Compare Update**: Improved comparison module
> - **Batch Output & Ranking**: Updated HTML output for runCASSIA_batch with new ranking method option
> - **Fuzzy Model Aliases**: Easier model selection without remembering exact names

<details>
<summary>📜 Previous Updates (click to expand)</summary>

> **2025-05-05**
> 📊 **CASSIA annotation benchmark is now online!**
> The latest update introduces a new benchmarking platform that evaluates how different LLMs perform on single-cell annotation tasks, including accuracy and cost.
> **LLaMA4 Maverick, Gemini 2.5 Flash, and DeepSeekV3** are the top three most balanced options—nearly free!
> 🔧 A new **auto-merge** function unifies CASSIA output across different levels, making subclustering much easier.
> 🐛 Fixed a bug in the annotation boost agent to improve downstream refinement.

> **2025-04-19**
> 🔄 **CASSIA adds a retry mechanism and optimized report storage!**
> The latest update introduces an automatic retry mechanism for failed tasks and optimizes how reports are stored for easier access and management.
> 🎨 **The CASSIA logo has been drawn and added to the project!**

> **2025-04-17**
> 🚀 **CASSIA now supports automatic single-cell annotation benchmarking!**
> The latest update introduces a new function that enables fully automated benchmarking of single-cell annotation. Results are evaluated automatically using LLMs, achieving performance on par with human experts.
> **A dedicated benchmark website is coming soon—stay tuned!**

</details>


## 🏗️ Installation

```R
# Install dependencies
install.packages("devtools")
install.packages("reticulate")

# Install CASSIA
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

***Note: If the environment is not set up correctly the first time, please restart R and run the code below***

```R
library(CASSIA)
setup_cassia_env()
```



### 🔑 Set Up API Keys

It should take about 3 minutes to get your API key.

We recommend starting with OpenRouter since it provides access to most models through a single API key.

```R
# For OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)

# For OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# For Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)

```


- **API Provider Guides:**
	- [How to get an OpenAI api key](https://platform.openai.com/api-keys)
	- [How to get an Anthropic api key](https://console.anthropic.com/settings/keys)
	- [How to get an OpenRouter api key](https://openrouter.ai/settings/keys)
    - [OpenAI API Documentation](https://beta.openai.com/docs/)
    - [Anthropic API Documentation](https://docs.anthropic.com/)
    - [OpenRouter API documentation](https://openrouter.ai/docs/quick-start)


## 🧬 Example Data

CASSIA includes example marker data in two formats:
```R
# Load example data
markers_unprocessed <- loadExampleMarkers(processed = FALSE)  # Direct Seurat output
markers_processed <- loadExampleMarkers(processed = TRUE)     # Processed format
```

## ⚙️ Pipeline Usage


```R
# The default provider is set to OpenRouter.

runCASSIA_pipeline(
    output_file_name = "cassia_test",            # Base name for output files
    tissue = "Large Intestine",                   # Tissue type (e.g., "brain")
    species = "Human",              		 # Species (e.g., "human")
    marker = markers_unprocessed,               # Marker data from findallmarker
    max_workers = 4                              # Number of parallel workers
)
```

## 🤖 Supported Models

You can choose any model for annotation and scoring. Some classic models are listed below. OpenRouter supports most of the current popular models, although some have not been extensively benchmarked in the CASSIA paper — feel free to experiment with them.



### OpenAI
- `gpt-5.1`: Balanced option (Recommended)
- `gpt-4o`: Used in the benchmark

### OpenRouter
- `google/gemini-2.5-flash`: One of the best-performing low-cost models, comparable with models like gpt-4o (Recommended)
- `deepseek/deepseek-chat-v3-0324`: One of the best-performing open-source models, which gives very detailed annotations
- `x-ai/grok-4-fast`: One of the best-performing low-cost models.

### Anthropic
- `claude-sonnet-4-5`: The latest best-performing model (Most recommended)

## 📤 Output

The pipeline generates four key files:
1. Complete Annotation Results CSV File
2. Annotation Summary HTML Report
3. Annotation Boost Agent Report for Low Quality Annotation

## 🧰 Troubleshooting

### Authentication (Error 401)
```R
# Check if API key is set correctly
key <- Sys.getenv("ANTHROPIC_API_KEY")
print(key)  # Should not be empty

# Reset API key if needed
setLLMApiKey("your_api_key", provider = "anthropic", persist = TRUE)
```

### File Errors
- Use absolute paths when necessary
- Check file permissions
- Ensure files aren't open in other programs

### Best Practices
- Keep API keys secure
- Maintain sufficient API credits


Note: This README covers only basic CASSIA functionality. For a complete tutorial including advanced features and detailed examples, please visit:
[CASSIA Complete Tutorial](https://docs.cassia.bio/en).

## 📖 Citation

📖 [Read our paper in Nature Communications](https://doi.org/10.1038/s41467-025-67084-x)

Xie, E., Cheng, L., Shireman, J. et al. CASSIA: a multi-agent large language model for automated and interpretable cell annotation. Nat Commun (2025). https://doi.org/10.1038/s41467-025-67084-x

## 📬 Contact

If you have any questions or need help, feel free to email us. We are always happy to help:
**xie227@wisc.edu**
If you find this project helpful, please share it with your friends, and give this repo a star ⭐
Many thanks!
