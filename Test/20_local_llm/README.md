# Test 20: Local LLM Support

## Purpose
Tests CASSIA's ability to work with local LLM servers (Ollama, LM Studio, vLLM) without requiring an API key. Includes a full batch annotation test using real marker data.

## Tests Performed

| Test | Description |
|------|-------------|
| 1. Localhost Detection | Verifies `localhost` and `127.0.0.1` are recognized |
| 2. API Key Bypass | Confirms no API key required for localhost URLs |
| 3. Remote Requires Key | Ensures non-localhost URLs still require API key |
| 4. Ollama Connection | Tests actual connection if Ollama is running (optional) |
| 5. Batch Annotation | Runs `runCASSIA_batch` with 2 clusters via local Ollama (optional) |

## Running the Tests

```bash
# Python
python test_local_llm.py

# R
Rscript test_local_llm.R
```

## Optional: Start Ollama

To run the full test including actual LLM calls and batch annotation:

```bash
# Install Ollama: https://ollama.ai
ollama serve  # Start server
ollama pull llama2  # Download a model (or any other model)
```

## Expected Results

- Tests 1-3 should always pass
- Tests 4-5 will be skipped if Ollama is not running
- Test 5 annotates 2 clusters: `monocyte` and `plasma cell`

## Supported Local LLM Servers

| Server | Default URL |
|--------|-------------|
| Ollama | `http://localhost:11434/v1` |
| LM Studio | `http://localhost:1234/v1` |
| vLLM | `http://localhost:8000/v1` |

All localhost URLs (any port) will work without requiring an API key.
