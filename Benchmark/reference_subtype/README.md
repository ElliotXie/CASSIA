# Reference-Assisted Subclustering Benchmarks

This suite benchmarks CASSIA subtype annotation when subclustering is assisted
by the reference agent. Everything here is for reference + subclustering
experiments; future unrelated benchmark tasks should live in sibling folders
under `Benchmark/`.

## Layout

```text
Benchmark/reference_subtype/
  README.md
  .openrouter_key                 # optional local key, ignored by git
  scripts/
    build_literature_cases.py
    macrophage_subtype_benchmark.py
    manual_evaluation.py
    t_cell_subcluster_benchmark.py
    hard_t_cell_marker_benchmark.py
    auto_split_t_cell_experiment.py
  cases/
    literature/
      coulton_2024_pan_cancer_tam_covered.csv
      coulton_2024_pan_cancer_tam_all.csv
      li_2023_uveal_melanoma_macrophage.csv
      li_2024_pan_cancer_icb_macrophage.csv
      wang_2023_prenatal_macrophage.csv
      literature_case_inventory.csv
    heldout/
      qi_2022_crc_macrophage.csv
  source_data/                    # downloaded benchmark source files, ignored
  results/
    <run_id>/
```

## Macrophage Benchmark

`scripts/macrophage_subtype_benchmark.py` compares three modes on the same
macrophage subtype marker panel:

1. `cassia_baseline`: `runCASSIA_subclusters(..., use_reference=False)`
2. `cassia_reference`: `runCASSIA_subclusters(..., use_reference=True)`
3. `direct_kimi`: one-shot `moonshotai/kimi-k2.6` call with no CASSIA pipeline
   and no reference agent

Generate literature-derived case CSVs:

```bash
python Benchmark/reference_subtype/scripts/build_literature_cases.py
```

The case builder uses up to the top 30 positive marker genes per paper subtype.
This is intentionally harder than top10 marker lookup and closer to realistic
subclustering inputs.

Current macrophage literature panels:

- `cases/literature/coulton_2024_pan_cancer_tam_covered.csv`: 12 Coulton 2024
  pan-cancer TAM clusters covered by the current human cancer macrophage
  consensus reference.
- `cases/literature/coulton_2024_pan_cancer_tam_all.csv`: broader non-doublet,
  non-NA Coulton 2024 TAM/monocyte-like cluster set.
- `cases/literature/li_2024_pan_cancer_icb_macrophage.csv`: 9 Li 2024
  pan-cancer ICB myeloid atlas macrophage clusters parsed from Supplementary
  Data 2. This is the best independent pan-cancer macrophage validation panel.
- `cases/literature/li_2023_uveal_melanoma_macrophage.csv`: 4 Li 2023 uveal
  melanoma macrophage clusters; useful for tissue-specific coverage gaps.
- `cases/literature/wang_2023_prenatal_macrophage.csv`: developmental
  macrophage panel, kept as a future-domain stress test rather than a human
  cancer consensus benchmark.
- `cases/heldout/qi_2022_crc_macrophage.csv`: 4 Qi 2022 colorectal cancer
  myeloid/macrophage cases from Supplementary Data 2. This paper is intentionally
  not included in the reference brain and should be used as a held-out check.
- `cases/literature/literature_case_inventory.csv`: local paper usability notes.

Run the macrophage benchmark:

```bash
OPENROUTER_API_KEY=... python Benchmark/reference_subtype/scripts/macrophage_subtype_benchmark.py \
  --cases-csv Benchmark/reference_subtype/cases/literature/li_2024_pan_cancer_icb_macrophage.csv
```

Outputs are written to `Benchmark/reference_subtype/results/<timestamp>/`.

For held-out evaluation, run the same benchmark on
`Benchmark/reference_subtype/cases/heldout/qi_2022_crc_macrophage.csv` and do
not add Qi 2022 to the macrophage reference brain.

## Evaluation

Automated metrics are triage only:

- `expected_terms` score: one point per expected term group found in
  `main_cell_type`, `sub_cell_type`, or `reason`.
- `paper_label_hit`: whether the output preserved the paper's author-defined
  subtype label such as `8_IFNGMac`.

The final benchmark call uses the fixed manual rubric in
`scripts/manual_evaluation.py`, with one row per case and mode:

- `program_accuracy` (0-4): biological subtype/program match.
- `subtype_specificity` (0-2): subtype-level label rather than generic
  macrophage/TAM.
- `marker_evidence` (0-2): correct use of defining and conflicting markers.
- `ambiguity_handling` (0-1): handles hybrid states and avoids overclaiming.
- `traceability` (0-1): output is traceable to a consensus program or paper
  alias.

Pass rule: total score at least 7/10 and `program_accuracy` at least 3/4.
Excellent rule: total score at least 9/10.

Summarize manual scores:

```bash
python Benchmark/reference_subtype/scripts/manual_evaluation.py \
  Benchmark/reference_subtype/results/<timestamp>/manual_scores.csv
```

A finished literature run should contain:

- `cases.csv`
- `inputs.csv`
- mode outputs such as `cassia_baseline.csv`, `cassia_reference.csv`,
  `direct_kimi.csv`
- `scores.csv` and `summary.json` for automated triage
- `usage.csv` for token and cost tracking
- `manual_scores.csv`, `manual_scores_scored.csv`,
  `manual_score_summary.csv`, `manual_score_summary.json`
- a human-readable evaluation markdown file

Cost tracking uses the provider-native `usage.cost` field when returned by
OpenRouter. Token and cost numbers should be compared within one run because
model routing and provider pricing can change over time.
