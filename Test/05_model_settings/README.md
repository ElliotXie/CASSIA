# Test 05: Model Settings

## Purpose
Tests model name resolution and provider shortcuts.

## What it Tests
- Simple name resolution (gemini -> google/gemini-2.5-flash)
- Quality shortcuts (best, cheap, flash)
- Practical annotation with OpenRouter + gemini-2.5-flash

## Resolution Tests
| Input | Provider | Expected |
|-------|----------|----------|
| gemini | openrouter | Contains "gemini" |
| best | openrouter | Contains "gemini" |
| cheap | openrouter | Resolves to valid model |
| flash | openrouter | Contains "flash" |

## Expected Output
- All resolution tests pass
- Practical annotation returns valid result
- Main cell type identified

## Running the Test

### Python
```bash
python test_model_settings.py
```

### R
```bash
Rscript test_model_settings.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Resolution and practical test results
