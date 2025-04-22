# CASSIA <img src="CASSIA_python/logo2.png" align="right" width="150" style="vertical-align: middle;" />

CASSIA (Collaborative Agent System for Single-cell Interpretable Annotation) is a tool that enhances cell type annotation using multi-agent Large Language Models (LLMs).

📖 [Read our preprint (v2, latest)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2)
 
📖 [Original preprint (v1, historical)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v1)

📝 [Example R workflow](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_tutorial_final.Rmd)

📝 [Example Python workflow](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_python_tutorial.ipynb)

📖 [Read our preprint (v2, latest)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2)
 
📖 [Original preprint (v1, historical)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v1)


## 📰 News

> **2025-04-19**  
> 🔄 **CASSIA adds a retry mechanism and optimized report storage!**  
> The latest update introduces an automatic retry mechanism for failed tasks and optimizes how reports are stored for easier access and management.  
> 🎨 **The CASSIA logo has been drawn and added to the project!**

> **2025-04-17**  
> 🚀 **CASSIA now supports automatic single-cell annotation benchmarking!**  
> The latest update introduces a new function that enables fully automated benchmarking of single-cell annotation. Results are evaluated automatically using LLMs, achieving performance on par with human experts.  
> **A dedicated benchmark website is coming soon—stay tuned!**


## 🏗️ Installation (R)

Install from GitHub
```R
# Install dependencies
install.packages("devtools")
install.packages("reticulate")

# Install CASSIA
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

### 🔑 Set Up API Keys

We recommend starting with OpenRouter since it provides access to most models through a single API key. While slightly more expensive and occasionally unstable, it offers greater convenience. For production use, direct access via OpenAI or Anthropic provides better stability.

Note that in certain countries, OpenAI and Anthropic may be banned. In these cases, users can use OpenRouter instead.

```R
# For OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# For Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)

# For OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```


- **API Provider Guides:**
	- [How to get an OpenAI api key](https://platform.openai.com/api-keys)
	- [How to get an Anthropic api key](https://console.anthropic.com/settings/keys)
	- [How to get an OpenRouter api key](https://openrouter.ai/settings/keys)
    - [OpenAI API Documentation](https://beta.openai.com/docs/)
    - [Anthropic API Documentation](https://docs.anthropic.com/)
    - [OpenRouter API documentatioon](https://openrouter.ai/docs/quick-start)


## 🧬 Example Data

CASSIA includes example marker data in two formats:
```R
# Load example data
markers_unprocessed <- loadExampleMarkers(processed = FALSE)  # Direct Seurat output
markers_processed <- loadExampleMarkers(processed = TRUE)     # Processed format
```

## ⚙️ Pipeline Usage

```R
runCASSIA_pipeline(
    output_file_name,     # Base name for output files
    tissue,               # Tissue type (e.g., "brain")
    species,              # Species (e.g., "human")
    marker,               # Marker data from findallmarker
    max_workers = 4,      # Number of parallel workers
    annotation_model = "gpt-4o",                    # Model for annotation
    annotation_provider = "openai",                 # Provider for annotation
    score_model = "anthropic/claude-3.5-sonnet",    # Model for scoring
    score_provider = "openrouter",                  # Provider for scoring
    annotationboost_model="anthropic/claude-3.5-sonnet", #model for annotation boost
    annotationboost_provider="openrouter", #provider for annotation boost
    score_threshold = 75,                          # Minimum acceptable score
    additional_info = NULL                         # Optional context information
)
```

## 🤖 Supported Models

You can choose any model for annotation and scoring. Some classic models are listed below. Most current popular models are supported by OpenRouter, though they have not been extensively benchmarked in the CASSIA paper — feel free to experiment with them.

### OpenAI (Most Common)
- `gpt-4o` (recommended): Balanced performance and cost
- `o1-mini`: Advanced reasoning capabilities (higher cost)

### Anthropic
- `claude-3-5-sonnet-20241022`: High-performance model
- `claude-3-7-sonnet-latest`: The latest model

### OpenRouter
- `anthropic/claude-3.5-sonnet`: High rate limit access to Claude
- `openai/gpt-4o-2024-11-20`: Alternative access to GPT-4o
- `meta-llama/llama-3.2-90b-vision-instruct`: Cost-effective open-source option
- `deepseek/deepseek-chat-v3-0324`: Very cost-effective and comparable to GPT-4o

## 📤 Output

The pipeline generates four key files:
1. Initial annotation results
2. Quality scores with reasoning
3. Summary report
4. Annotation boost report

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
- Verify sufficient disk space

### Best Practices
- Keep API keys secure
- Maintain sufficient API credits
- Backup data before overwriting files
- Double-check file paths and permissions

## 🛠️ To-Do List

- [ ] Imporve the MOE (Mixture of Experts) system
- [ ] Demo video
- [x] Better integration with Seurat
- [ ] Better integration with Scanpy
- [ ] Integration with clustering pipeline  
- [ ] Marker selection process optimization
- [ ] Multiomics integration
- [ ] Laucnch the benchmark website
- [x] More robust workflow with auto-retry mechanisms  
- [x] Improved output file management
- [x] Automatic evalutate annotation performance
- [x] Draw the CASSIA logo

Note: This README covers basic CASSIA functionality. For a complete tutorial including advanced features and detailed examples, please visit:
[CASSIA Complete Tutorial](https://cassia-true-final-4.vercel.app/).

## 📖 Citation

CASSIA: a multi-agent large language model for reference free, interpretable, and automated cell annotation of single-cell RNA-sequencing data  
Elliot Xie, Lingxin Cheng, Jack Shireman, Yujia Cai, Jihua Liu, Chitrasen Mohanty, Mahua Dey, Christina Kendziorski  
bioRxiv 2024.12.04.626476; doi: https://doi.org/10.1101/2024.12.04.626476
