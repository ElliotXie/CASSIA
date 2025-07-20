# CASSIA Test Code

This directory contains focused test notebooks for specific CASSIA features.

## Available Tests

### `test_model_settings.ipynb`
Comprehensive test of the CASSIA Model Settings System:

- **Configuration Testing**: Verifies model settings load from package data directory
- **Model Resolution**: Tests simple name resolution (gpt4 → gpt-4o)
- **Provider Requirements**: Ensures provider must be specified for API key control
- **Quality Shortcuts**: Tests 'best', 'cheap', 'premium', 'fast' per provider
- **Practical Examples**: Shows real-world usage in CASSIA functions
- **Migration Guide**: Demonstrates how to upgrade from complex to simple names

## How to Run

1. **Open the notebook**: Navigate to the test you want to run
2. **No installation required**: Tests use the local CASSIA code directly
3. **Run all cells**: Execute the notebook from top to bottom
4. **Review results**: Each test shows pass/fail status and detailed output

## Key Features Tested

### Model Settings System
- ✅ Provider control (openai, anthropic, openrouter)
- ✅ Simple names (gpt4, claude, gemini)
- ✅ Quality shortcuts (best, cheap, premium, fast)
- ✅ Package data directory integration
- ✅ Error prevention (no accidental API switching)
- ✅ Backward compatibility

### Usage Examples
```python
# Simple names
runCASSIA_batch(
    marker=markers,
    model="gpt4",         # Instead of "gpt-4o"
    provider="openai"     # Required for API key control
)

# Quality shortcuts
runCASSIA_batch(
    marker=markers,
    model="cheap",        # Cheapest option for provider
    provider="openrouter"
)

# Mixed providers
runCASSIA_pipeline(
    annotation_model="gemini",
    annotation_provider="openrouter",
    score_model="gpt4",
    score_provider="openai"
)
```

## Benefits

- **Focused Testing**: Each notebook tests specific features
- **Clear Output**: Easy to see what works and what doesn't
- **No Setup**: Uses local code without package installation
- **Comprehensive Coverage**: Tests all major functionality
- **Practical Examples**: Shows real usage patterns

## Next Steps

After tests pass:
1. Use simple model names in your CASSIA workflows
2. Always specify provider for API key control
3. Use quality shortcuts for easy model selection
4. Enjoy simplified model management!