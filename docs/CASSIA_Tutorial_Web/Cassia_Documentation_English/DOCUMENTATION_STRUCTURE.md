# CASSIA Documentation Structure

## Overview

- **Languages**: English (`en/`) and Chinese (`zh/`) - always update both
- **Packages**: R and Python - keep in sync when feature exists in both

## Directory Structure

```
content/
├── en/                     # English
│   ├── docs/               # API reference (r/, python/)
│   └── vignette/           # Tutorials (r/, python/)
└── zh/                     # Chinese (same structure)
```

## Key Rules

1. Update both `en/` and `zh/` for any change
2. Use ` ```bash ` for shell commands, not ` ```python ` with `!pip`
3. Model presets: `premium` and `budget` (see `CASSIA_python/CASSIA/comparison/model_config.json`)

## File Naming Pattern

Files are named consistently across R and Python:
- `setting-up-cassia.md` - Installation
- `basic-annotation.md` - Core annotation
- `quality-scoring-and-report-generation.md` - Scoring
- `introduction-to-optional-agents.md` - Agent overview
- `annotation-boost.md`, `compare-cell-types.md`, `subclustering.md`, `uncertainty-quantification.md` - Optional agents
- `ragagent.md` - Python only
