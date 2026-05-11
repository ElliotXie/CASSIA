# Reference Subtype Benchmarks

This folder stores benchmarks for CASSIA reference-assisted subtype annotation.

## Macrophage Subtype Benchmark

`macrophage_subtype_benchmark.py` compares three modes on the same macrophage
subtype marker panel:

1. `cassia_baseline`: `runCASSIA_subclusters(..., use_reference=False)`
2. `cassia_reference`: `runCASSIA_subclusters(..., use_reference=True)`
3. `direct_kimi`: one-shot `moonshotai/kimi-k2.6` call with no CASSIA
   annotation pipeline and no reference agent

The benchmark uses expected subtype term groups with synonyms. A case receives
one point per expected term group found in `main_cell_type`, `sub_cell_type`, or
`reason`; each case has a `min_score` threshold. Literature cases also report
`paper_label_hit`, a stricter metric that checks whether the output preserved
the paper's author-defined subtype label such as `8_IFNGMac`.

Literature-derived case CSVs are generated from locally downloaded paper
supplementary tables:

```bash
python Benchmark/reference_subtype/build_literature_cases.py
```

The first real benchmark panel is:

- `literature_cases/coulton_2024_pan_cancer_tam_covered.csv`: 12 Coulton 2024
  pan-cancer TAM clusters covered by the current macrophage reference brain.
- `literature_cases/coulton_2024_pan_cancer_tam_all.csv`: broader non-doublet,
  non-NA Coulton 2024 TAM/monocyte-like cluster set.
- `literature_cases/literature_case_inventory.csv`: local paper usability notes.

Run:

```bash
OPENROUTER_API_KEY=... python Benchmark/reference_subtype/macrophage_subtype_benchmark.py
```

Run with the literature panel:

```bash
OPENROUTER_API_KEY=... python Benchmark/reference_subtype/macrophage_subtype_benchmark.py \
  --cases-csv Benchmark/reference_subtype/literature_cases/coulton_2024_pan_cancer_tam_covered.csv
```

Outputs are written to `Benchmark/reference_subtype/results/<timestamp>/`:

- `cases.csv`
- `inputs.csv`
- `cassia_baseline.csv`
- `cassia_reference.csv`
- `direct_kimi.csv`
- `scores.csv`
- `usage.csv`
- `summary.json`

For literature benchmarks, treat these automated scores as triage. The final
call should include a human-readable evaluation that compares each generated
answer with the paper subtype and marker evidence, for example
`results/20260511_152028/human_evaluation.md`.

Cost tracking uses the provider-native `usage.cost` field when returned by
OpenRouter. Token and cost numbers should be compared within one run because
model routing and provider pricing can change over time.
