# Core `runCASSIA_batch` Benchmark

This suite benchmarks the default CASSIA annotation pipeline through the public
`runCASSIA_batch()` API. It is separate from `Benchmark/reference_subtype/`,
which targets reference-assisted subclustering.

## Scope

The default case set has 50 marker panels selected from the bundled historical
marker-list workbook:

```text
CASSIA_example/Benchmark/marker100.xlsx
```

The selection prioritizes GTEx and Azimuth, then adds selected Tabula Sapiens
and HCL cases for broader tissue/subtype coverage:

| Source | Tissue | Cases | Role |
| --- | --- | ---: | --- |
| GTEx | lung | 15 | high-quality lung epithelial, stromal, endothelial, and immune cases |
| Azimuth | kidney | 22 | high-quality kidney immune, epithelial, stromal, endothelial, and renal-segment cases |
| Tabula Sapiens | large intestine | 8 | non-eye intestinal epithelial/immune stress cases |
| HCL | fetal skin | 5 | fetal skin stromal, endothelial, immune, and melanocyte cases |

The suite intentionally mixes difficulty levels:

- `easy`: canonical lineages such as B cell, mast cell, NK cell, macrophage,
  neutrophil, fibroblast, podocyte, and proximal tubule.
- `medium`: tissue-specialized epithelial, endothelial, myeloid, fibroblast,
  melanocyte, and renal tubule subtypes.
- `hard`: ambiguous or fine-grained cases such as DC/macrophage, Pericyte/SMC,
  CD8/NKT distinction, myofibroblast, kidney endothelial subtype resolution,
  connecting tubule, intestinal tuft/stem/transit-amplifying/Paneth states.

The selected 50 rows have historical baseline metadata from
`Supplementary_Data 4.xlsx` when available:

| Method | Weighted correctness | Full-correct rows |
| --- | ---: | ---: |
| CASSIA | 46.5/50 = 0.930 | 43/50 |
| GPTCellType-4o | 37/50 = 0.740 | 30/50 |
| GPTCellType-4 | 32/50 = 0.640 | 24/50 |

## Layout

```text
Benchmark/core_run_cassia_batch/
  README.md
  .openrouter_key                 # optional local key, ignored by git
  cases/
    core_50.csv                   # stable benchmark case manifest
    source_inventory.csv          # all 100 marker100 rows with selected flag
  scripts/
    build_core_cases.py           # rebuilds core_50.csv from marker100.xlsx
    run_core_batch_benchmark.py   # grouped runCASSIA_batch runner
    score_core_batch_results.py   # deterministic term-hit scorer
  results/
    <run_id>/
```

## Case Manifest

`cases/core_50.csv` is the source of truth for benchmark execution. Important
columns:

- `case_id`: stable cluster ID passed into `runCASSIA_batch`.
- `dataset`, `tissue`, `species`: grouping metadata. The runner calls
  `runCASSIA_batch` once per dataset/tissue/species group because the core API
  accepts one tissue/species context per batch.
- `expected_cell_type`, `broad_cell_type`: hidden ground truth for scoring and
  review. These are not sent to CASSIA.
- `difficulty`, `case_type`: analysis strata.
- `expected_terms`, `min_score`: deterministic triage scoring rubric.
- `source_cassia_correctness`, `source_gptcelltype_4o_correctness`,
  `source_gptcelltype_4_correctness`: historical baseline correctness copied
  from `Supplementary_Data 4.xlsx` when matched.
- `marker_list`: comma-separated markers used as the batch input.

## Rebuild Cases

The checked-in manifest is a stable snapshot. To rebuild it from
`marker100.xlsx`:

```bash
python Benchmark/core_run_cassia_batch/scripts/build_core_cases.py
```

This also writes `cases/source_inventory.csv`, which keeps all 100 source rows
and marks whether each row was selected for the 50-case core suite.

## Run

Set an API key with either an environment variable or a local ignored key file:

```bash
export OPENROUTER_API_KEY=...
# or
printf '%s\n' '...' > Benchmark/core_run_cassia_batch/.openrouter_key
```

Dry-run the selected groups without API calls:

```bash
python Benchmark/core_run_cassia_batch/scripts/run_core_batch_benchmark.py --dry-run
```

Run the full 50-case benchmark:

```bash
python Benchmark/core_run_cassia_batch/scripts/run_core_batch_benchmark.py \
  --model openai/gpt-5.6-terra \
  --provider openrouter \
  --max-workers 5
```

Useful smaller runs:

```bash
python Benchmark/core_run_cassia_batch/scripts/run_core_batch_benchmark.py --difficulty hard
python Benchmark/core_run_cassia_batch/scripts/run_core_batch_benchmark.py --limit 5
```

Outputs are written to `results/<run_id>/`:

- `cases.csv`: exact cases used.
- `inputs.csv`: two-column CASSIA batch input.
- `<group>/cassia_batch_summary.csv`: raw `runCASSIA_batch` summary per group.
- `combined_predictions.csv`: merged predictions plus case metadata.
- `scores.csv`: deterministic term-hit scoring.
- `score_summary.csv` and `summary.json`: run-level summaries.
- `usage.csv`: CASSIA usage accounting by group.

## Evaluation Notes

The deterministic scorer is triage only. It answers: did the prediction contain
the expected broad/specific biological terms? It should not be used as the final
paper-quality judgement for ambiguous cases. Hard cases should be manually
reviewed, especially when the model gives a biologically plausible broader call
without the exact subtype term.

The marker lists come from `marker100.xlsx`. The older
`Supplementary_Data 4.xlsx` table is used only for historical baseline metadata,
because it does not contain the marker lists needed by `runCASSIA_batch`.
