## R CMD check results

0 errors | 0 warnings | 2 notes

## First submission

This is the first submission of CASSIA to CRAN.

## Publication

This package accompanies a paper accepted for publication in Nature Communications.

## Notes explanation

### Package size (NOTE)

The installed package size is approximately 17 MB. This is because CASSIA bundles a complete Python module in `inst/python/` that provides the LLM-based cell type annotation engine. This Python code is essential for the package's core functionality and is loaded via the reticulate package. The size is necessary to provide a self-contained annotation system that works across platforms.

### Non-standard file at top level (NOTE)

The `cran-comments.md` file is intentionally included for CRAN reviewer communication.

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
