# CASSIA Benchmarks

This directory is organized by benchmark task. Each task folder owns its case
files, runners, scoring code, and results.

## Suites

- `cli_release/`: lightweight, auditable benchmark snapshot for the CASSIA
  1.4 CLI release. It contains fixed 50-, 100-, and 392-case manifests plus
  compact result records; raw expression matrices, full marker tables, model
  transcripts, and large run directories are intentionally excluded.
- `core_run_cassia_batch/`: core `runCASSIA_batch()` benchmark with 50 curated
  marker panels, a deterministic scorer, and a grouped public API runner.
- `reference_subtype/`: reference-assisted subclustering and subtype annotation
  benchmarks. This currently includes macrophage/TAM and T-cell subcluster
  experiments, plus the fixed manual rubric for subtype-call evaluation.
- `integrated_agent/`: pilot experiments for the transactional, no-daemon
  clustering-and-annotation controller, including fixed/conservative/adaptive
  ablations and the policy iteration that prevents split fragmentation.

Future benchmark tasks should be added as sibling folders rather than mixed
into an existing suite.
