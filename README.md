# CASSIA

CASSIA (Collaborative Agent System for Single cell Interpretable Annotation) is a tool that enhances cell type annotation using multi-agent Large Language Models (LLMs).

📖 [Read our paper](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v1) for detailed methodology and benchmarking results.

📝 [Example R workflow](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_tutorial_final.Rmd)

📚 [Complete R Documentation](https://cassia-true-final-4.vercel.app/)

📝 [Example Python workflow](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_python_tutorial.ipynb)

🌐 [Try CASSIA Web UI](https://cassiacell.com/) - A web interface for basic CASSIA functionality

## 📰 News

> **2025-04-17**  
> 🚀 **CASSIA now supports automatic single cell annotation benchmarking!**  
> The latest update introduces a new function that enables fully automated benchmarking of single cell annotation. Results are evaluated automatically using LLMs, achieving performance on par with human experts.  
> **A dedicated benchmark website is coming soon—stay tuned!**

### Installation

Option 1: Install from GitHub
```R
# Install dependencies
install.packages("devtools")
install.packages("reticulate")

# Install CASSIA
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

Option 2: Install from source
```R
install.packages("reticulate")
install.packages("remotes")
remotes::install_url("https://github.com/ElliotXie/CASSIA/raw/main/CASSIA_source_R/CASSIA_0.1.0.tar.gz")
```

### Set Up API Keys

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


## Example Data

CASSIA includes example marker data in two formats:
```R
# Load example data
markers_unprocessed <- loadExampleMarkers(processed = FALSE)  # Direct Seurat output
markers_processed <- loadExampleMarkers(processed = TRUE)     # Processed format
```

## Pipeline Usage

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

## Supported Models

### OpenAI (Most Common)
- `gpt-4o` (recommended): Balanced performance and cost
- `gpt-4o-mini`: Faster, more economical option
- `o1-mini`: Advanced reasoning capabilities (higher cost)

### Anthropic
- `claude-3-5-sonnet-20241022`: High-performance model

### OpenRouter
- `anthropic/claude-3.5-sonnet`: High rate limit access to Claude
- `openai/gpt-4o-2024-11-20`: Alternative access to GPT-4o
- `meta-llama/llama-3.2-90b-vision-instruct`: Cost-effective open-source option

## Output

The pipeline generates four key files:
1. Initial annotation results
2. Quality scores with reasoning
3. Summary report
4. Annotation boost report

## Troubleshooting

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

- [ ] Better integration with Seurat directly  
- [ ] Integration with Clustering Pipeline  
- [ ] Marker selection process optimization  
- [ ] More robust workflow with auto-retry mechanisms  
- [ ] Improved output file management
- [x] Automatic benchmarking

Note: This README covers basic CASSIA functionality. For a complete tutorial including advanced features and detailed examples, please visit:
[CASSIA Complete Tutorial](https://cassia-true-final-4.vercel.app/).
