# Held-out validation protocol

This directory is a development-only evaluation layer for the integrated
clustering + annotation agent. It is not imported by the CASSIA package and is
not part of a release path.

## What the existing benchmark can validate

The existing marker benchmarks remain valuable, but they answer only the
annotation half of the question:

- `Benchmark/cli_release/cases/full_392.csv`: 392 cases from 39 datasets; use
  as the main annotation regression set.
- `Benchmark/reference_subtype/heldout/zhang_crc_2018`: 20 genuinely blinded
  T-cell subtype cases with separate marker input and truth files.
- `Benchmark/reference_subtype/cases/heldout/qi_2022_crc_macrophage.csv`: four
  additional macrophage/myeloid subtype cases. Split marker inputs from truth
  before a run.
- `Benchmark/core_run_cassia_batch/cases/core_50.csv`: fast smoke/regression
  subset.

These files contain cluster marker lists, not cell-by-gene matrices or cell
memberships. They cannot establish whether a merge or split was biologically
correct and must not be used to report clustering accuracy.

## Integrated cell-level cohort

Restore a small, preregistered subset of the source matrices referenced by
`full_392.csv`. A useful first cohort is:

1. `pbmc_hao_citeseq`: human blood, orthogonal CITE-seq labels, many related
   immune subtypes; primary topology stress test.
2. `intestine_ts_2022`: human non-immune tissue with expert/consensus labels;
   guards against optimizing only for PBMC structure.
3. One mouse atlas such as `cerebellum_mouse_2021` or
   `pancreas_mouse_2023`; checks cross-species behavior.

Start with one dataset while the harness is being stabilized, then freeze all
three before comparing policies. The checked-in repository records source
accessions and intended `data/raw/*.h5ad` paths, but does not distribute the
large matrices.

For the first cohort, avoid downloading the 2.43 GiB source H5AD. The installed
R Census client can fetch only a deterministic balanced slice (24 labels, up to
150 cells per label):

```bash
Rscript Benchmark/integrated_agent/heldout/prepare_census_cohort.R \
  --dataset-id ed5d841d-6346-47d4-ab2f-7119ad7e3a35 \
  --labels Benchmark/integrated_agent/heldout/cohorts/pbmc_hao_citeseq_labels.csv \
  --out-dir data/heldout/pbmc_hao_citeseq
```

This writes `blind/input.rds`, `evaluator/truth.csv`, and a hashed cohort
manifest. The allowlist is derived from the existing `full_392.csv` PBMC cases.
Because the slice is label-balanced, treat it as a topology stress test rather
than a prevalence-weighted estimate of performance on unselected PBMCs.

## Leakage boundary

For every dataset, construct a blind Seurat RDS containing expression and only
allowed covariates (for example donor/batch), with paper labels, ontology
labels, source cluster names, expected terms, and label-provenance columns
removed. The agent receives only that RDS plus species/tissue. Keep the truth
file out of the agent run directory and do not mention truth column names in
the prompt. Run the scorer only after the automation process exits.

For a contamination-resistant benchmark, the evaluation operator should keep
truth outside the agent's readable sandbox. Merely hiding a column name is not
an adversarial isolation boundary.

## Comparison arms

Use identical input RDS, initial expression-derived partition, seed, variable
features, PCs, and model settings for all arms:

- `fixed`: annotation without topology edits.
- `conservative`: at most two guarded local edits; current candidate.
- `adaptive`: bounded exploratory comparator, not an assumed winner.

Repeat stochastic LLM arms at least three times. Report every run rather than
selecting the best completion.

## Metrics

The post-run scorer joins held-out truth to
`provenance/final_memberships.csv` by `cell_id` and propagates cluster labels
from `outputs/annotation.tsv` to cells. It reports:

- topology: ARI, NMI, homogeneity, completeness, V-measure, and purity;
- annotation: fine and broad cell-level accuracy, macro-F1, balanced accuracy,
  and coverage (skipped clusters receive no credit);
- fragmentation: predicted/truth cluster ratio, smallest cluster, tiny-cluster
  count, and fraction of cells in tiny clusters.

Free-form predicted labels should be mapped after the run by a locked judge or
ontology resolver. Store that frozen adjudication as:

```csv
predicted_label,canonical_fine_label,canonical_broad_label
CD14+ monocyte,CD14-positive monocyte,Monocyte
```

Without `--label-map`, the scorer deliberately falls back to conservative
normalized exact matching and records that weaker scoring mode.

```bash
python Benchmark/integrated_agent/heldout/score_integrated.py \
  --truth /evaluator/pbmc_hao_truth.csv \
  --run runs/pbmc_hao/conservative/seed_1 \
  --label-map /evaluator/pbmc_hao_label_map.csv \
  --out evaluations/pbmc_hao/conservative/seed_1
```

Truth CSV schema:

```csv
cell_id,truth_fine_label,truth_broad_label
AAAC...,CD14-positive monocyte,Monocyte
```
