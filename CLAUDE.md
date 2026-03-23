# CLAUDE.md

## Project Overview

CASSIA (Collaborative Agent System for Single-cell Interpretable Annotation) is a multi-agent LLM framework for automated, reference-free cell type annotation in single-cell RNA-seq data. Maintains both Python and R implementations sharing core Python logic.

## Project Structure

| Directory | Description |
|-----------|-------------|
| `CASSIA_python/CASSIA/` | Python package (PyPI) - primary source |
| `CASSIA_R/R/` | R package wrappers (GitHub) |
| `CASSIA_R/inst/python/` | Python code synced to R package |
| `Test/` | Test suite with `shared/` utilities |
| `docs/cassia_javascript/cassia-web` | Next.js web app for online use | cassia.bio
| `docs/CASSIA_Tutorial_Web/` | Documentation and vignettes |doc.cassia.bio
| `docs/README.md`, `README_CN.md` | English/Chinese READMEs |
| 'doccs/CASSIA-Benchmark' | CASSIA Benchmark |

## Multi-Agent System

**Core Agents:** Annotation (`main_function_code.py`) → Validation → Formatting → Scoring (`LLM_evalution.py`) → Reporting (`generate_reports.py`)

**Optional Agents:** Annotation Boost, Subclustering, Uncertainty Quantification, Symphony Compare, Merging

**Key Modules:** `tools_function.py` (API), `llm_utils.py` (LLM interface), `model_settings.py` (config)

## Development Workflow

1. Edit Python code in `CASSIA_python/CASSIA/`
2. Add R wrapper in `CASSIA_R/R/` if needed
3. Create tests in `Test/XX_feature_name/` (both `.py` and `.R`)
4. Update `__init__.py` exports
5. **Run `sync_python_files.py`** to sync Python code to R package
6. **Update docs** in `docs/CASSIA_Tutorial_Web/` — parameter tables, code examples in both en/zh and python/r versions
7. **Bump version** for PyPI release — update version in both `setup.py` and `__init__.py` (patch version for small changes)
