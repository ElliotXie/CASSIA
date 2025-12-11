# Test 17: Provider Defaults

Tests the `overall_provider` parameter in `runCASSIA_pipeline()` which automatically sets provider-specific default models for each pipeline stage.

## What This Tests

The `overall_provider` parameter allows users to specify a single provider and have all pipeline stages automatically use appropriate default models:

| Stage | OpenAI | Anthropic | OpenRouter |
|-------|--------|-----------|------------|
| Annotation | gpt-5.1 | claude-sonnet-4-5 | openai/gpt-5.1 |
| Score | gpt-5.1 | claude-sonnet-4-5 | anthropic/claude-sonnet-4.5 |
| Merge | gpt-5-mini | claude-haiku-4-5 | google/gemini-2.5-flash |
| AnnotationBoost | gpt-5.1 | claude-sonnet-4-5 | anthropic/claude-sonnet-4.5 |

## Test Files

| File | Description |
|------|-------------|
| `test_provider_defaults.py` | Python development mode test |
| `test_provider_defaults_install.py` | Python pip-installed package test |
| `test_provider_defaults.R` | R development mode test |
| `test_provider_defaults_install.R` | R installed package test |

## Usage

### Python Development Mode
```bash
python test_provider_defaults.py
```

### Python Install Mode
```bash
python test_provider_defaults_install.py
```

### R Development Mode
```bash
Rscript test_provider_defaults.R
```

### R Install Mode
```bash
Rscript test_provider_defaults_install.R
```

## Test Configuration

- **Provider tested**: `openrouter` (default)
- **Clusters**: monocyte, plasma cell (2 clusters for faster testing)
- **Score threshold**: 99 (to avoid triggering annotation boost)
- **Data**: Uses marker data from `Test/16_cassia_pipeline/data/`

## Expected Output

The test validates that:
1. Pipeline runs successfully with `overall_provider` parameter
2. Default models are automatically selected for each stage
3. Output directory structure is created correctly
4. FINAL_RESULTS.csv is generated

## Results

Results are saved to:
- `results/python/development/` - Python dev mode results
- `results/python/install/` - Python pip install results
- `results/r/development/` - R dev mode results
- `results/r/install/` - R install mode results
