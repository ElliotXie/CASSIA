# Integrated clustering + annotation pilot

This pilot exercised `cassia agent auto` on Seurat's 80-cell `pbmc_small`
teaching object with Codex CLI, `gpt-6-astra`, and high reasoning effort. Its
purpose was to test controller behavior, transaction safety, and topology
policy design. It is not a biological-accuracy benchmark: the fixture has no
held-out truth labels, and a CASSIA QA pass verifies evidence/rules rather than
ground-truth correctness.

## Iterations

| arm | agent commands | committed edits | clusters | labeled clusters | labeled cells | smallest cluster | outcome |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| fixed | 13 | 0 | 3 → 3 | 2/3 | 44/80 | 19 | Mixed platelet/T/NK parent correctly left unresolved |
| conservative, manual resolution | 14 | 0 | 3 → 3 | 2/3 | 44/80 | Resolutions 0.3 and 0.6 were no-ops; controller stopped early |
| adaptive, unguarded | 47 | 3 | 3 → 10 | 6/10 | 60/80 | Recovered detail but produced four clusters below 5 cells |
| conservative, guarded auto-resolution | 31 | 2 | 3 → 5 | 5/5 | 80/80 | Recovered platelet, NK, and CCR7/IL7R T groups; smallest cluster 9 cells |

All four completed runs passed the rule-based QA gate. The final conservative
candidate was selected for the default policy because it resolved the mixed
parent without the adaptive arm's fragmentation. That conclusion is limited to
this pilot and should be re-tested on larger datasets with truth labels.
The compact machine-readable record is in `pilot_results.json`.

## Design changes driven by the pilot

- Default execution is `agent → CLI → one-shot R transaction → checkpoint`;
  daemon mode is optional only.
- A failed command cannot advance the state pointer or spend topology budget.
- `subcluster --auto-resolution` tries `0.3,0.6,0.8,1.0` and commits the first
  policy-valid result.
- Conservative mode allows at most two local topology edits, rejects children
  below 5 cells, and rejects a split into more than 3 children.
- Baseline/final memberships and cluster registries are exported for exact
  provenance; comparison metrics are computed from the audit and annotation
  TSV rather than trusted from model prose.

## Reproduce

```bash
cassia agent auto object.rds \
  --out runs/conservative \
  --strategy conservative \
  --backend codex-cli \
  --model gpt-6-astra \
  --reasoning-effort high \
  --max-topology-edits 2 \
  --min-child-cells 5 \
  --max-children-per-split 3

cassia agent compare runs/fixed runs/adaptive runs/conservative \
  --out runs/comparison.md
```

Each run records its input SHA-256, model settings, prompt, immutable state
versions, audit log, recovery checkpoint, exact memberships, QA result, and
final annotated RDS/TSV/report. Large run directories and the Seurat object are
intentionally not checked into the repository.

## Next validation stage

The pilot above validates execution and policy behavior only. The development
protocol and post-run cell-level scorer for biological held-out evaluation are
in [`heldout/README.md`](heldout/README.md). Existing marker-only benchmarks
remain annotation regression gates; source matrices are required before ARI,
NMI, or cell-level integrated accuracy can be claimed.
