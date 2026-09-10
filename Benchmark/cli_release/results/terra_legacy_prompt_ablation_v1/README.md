# Terra High legacy-prompt ablation — Full392

Status: **complete, published release benchmark record**
Record version: `terra-legacy-prompt-ablation-full392-v1`

## Question

Does the verbatim 1,602-character default CASSIA step-by-step prefix add accuracy
to the active-evidence Fused Boost v2 workflow?

## Locked design

- Dataset: 392 fixed cases; published selection SHA-256
  `69ed284887cb5575c48af26cc7374e79f35c27566fc8fd9c4bac7158cb02caa3`.
  The manifest also retains the original local source-selection hash from before
  absolute paths were converted to repository-relative provenance paths.
- Control: Terra High Fused Boost v2.
- Treatment: remove only the legacy prefix. Active evidence, breadth search,
  unlimited gene requests, five-round limit, and isolated workspace are unchanged.
- Annotation: `gpt-5.6-terra`, high reasoning, four workers.
- Judge: `stable-judge-v1.1`, `gpt-5.6-sol`, medium reasoning.
- Scoring unit: one paired `case_id`; main cell type and top-1 subtype determine
  core credit. Top-2/top-3 are diagnostic only.
- Existing control and completed Boost-100 checkpoints were reused. Only the
  remaining 292 treatment annotations were newly run.

## Results

| Arm | Correct | Partial | Wrong | Mean score | Strict | Calls | Queries | USD equiv. | USD/cluster |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fused Boost v2 | 296 | 36 | 60 | 0.8277 | 209 | 1164 | 754 | $39.93 | $0.1019 |
| Without legacy CoT block | 296 | 39 | 57 | 0.8296 | 215 | 1374 | 979 | $33.41 | $0.0852 |

Paired ordinal transitions were 25 improved,
22 regressed, and 345 tied.
Exact Correct discordance was 19 treatment-only
versus 19 control-only
(`p=1.000000`). Strict discordance was
40 versus 34
(`p=0.561381`).

## Decision

The legacy step-by-step prefix has **no measurable accuracy benefit** on Full392:
both arms produced 296/392 Correct calls and the mean-score
difference was only
`+0.0019`.
The compact prompt is therefore the cleaner release default.

The compact arm's annotation credit-equivalent cost was
`16.3%` lower, despite more calls and marker-query rounds, because
it emitted substantially fewer output tokens. Judge cost is excluded, and no
incremental subscription charge is inferred from this estimate.

## Files

- `metrics.csv` — two-arm accuracy, resource, token, and cost metrics.
- `paired_summary.csv` — paired transitions and exact McNemar tests.
- `manifest.json` — protocol, published and source selection hashes, raw-file hashes, and
  rebuild command.
- Raw per-case predictions, scores, and transcripts remain under
  `results/terra_legacy_cot_ablation_*`; those large run directories are
  intentionally ignored by Git.

## Important comparison rule

Do not compare the 296 Correct count here directly with the Luna-Max model
comparison's Terra count. This ablation used a separate Sol-Medium judge run.
Only within-record paired comparisons are valid.

## Limitation

Annotation arms were produced at different times, so model nondeterminism remains
a residual confounder. The paired result nevertheless shows no directional exact-
accuracy advantage: Correct discordance was exactly 19 versus 19.
