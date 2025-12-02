## R CMD check results

0 errors | 0 warnings | X notes

(Will be updated after running R CMD check)

## First submission

This is the first submission of CASSIA to CRAN.

## Package notes

### Python dependency

CASSIA uses Python via the reticulate package for its core LLM-based annotation functionality. The package:

- Gracefully handles cases where Python is not available
- All examples use `\dontrun{}` blocks since they require API keys
- Tests use `skip_on_cran()` to avoid Python-dependent tests on CRAN

### API key requirements

The package requires LLM API keys (OpenAI, Anthropic, or OpenRouter) for core functionality. This is documented in the package description and all examples appropriately use `\dontrun{}`.

## Test environments

* Local: Windows 11, R 4.4.1
* (Additional platforms will be listed after rhub testing if performed)
