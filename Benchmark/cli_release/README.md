# CASSIA CLI release benchmark snapshot

This directory records the benchmark evidence used to prepare the CASSIA 1.4
CLI release. It is deliberately small enough for Git while preserving the
fixed case definitions, summary results, scoring protocol, and hashes needed
to audit release decisions.

## Fixed case sets

| Manifest | Cases | Purpose |
| --- | ---: | --- |
| `cases/model_card_50.csv` | 50 | Fast cross-model diagnostic card |
| `cases/boost_100.csv` | 100 | Paired Annotation Boost stress test |
| `cases/full_392.csv` | 392 | Full reference-free annotation benchmark |

The full set spans 199 cell types, 28 tissues, human and mouse, and
dissociated, multimodal, and spatial data. Each manifest retains label
provenance, evidence tier, publication timing, circularity risk, marker input,
and scoring granularity. Local absolute data paths have been converted to
repository-relative provenance paths; the large source matrices are not
distributed here.

## Release-facing findings

### Selective Annotation Boost

On the fixed, deliberately enriched Boost-100 diagnostic set, one-shot
annotation scored 50/100 exact and 67/100 exact-or-partial. Review Boost scored
71/100 exact and 85/100 exact-or-partial, with 26 strict rescues and five strict
regressions. Because the set oversamples known hard and boost-recoverable
cases, the 21-point difference is not a prevalence-weighted estimate of
real-world uplift. It supports exposing Boost as a selective escalation rather
than applying it unconditionally. See `results/boost_100.md`.

### Compact Fused Boost prompt

A paired 392-case ablation compared the legacy Fused Boost v2 prompt with a
version that removed only its duplicated 1,602-character step-by-step prefix.
Both arms produced 296/392 exact-correct calls. Correct discordance was 19
versus 19 (exact McNemar p=1.0), while estimated annotation cost fell from
$39.93 to $33.41 (16.3%). CASSIA 1.4 therefore uses `v2-compact` as the release
default and retains `v2` for historical reproduction. See
`results/terra_legacy_prompt_ablation_v1/`.

### Cross-model model card

The 50-case card includes single-run results for Composer 2.5, Claude Opus 4.8
xhigh, and GPT-5.5 high under one fixed Composer 2.5 referee. These results are
a capability profile, not a permanent leaderboard: model nondeterminism is
material, and gaps of roughly five points or less should be treated as noise.
See `results/model_card_comparison.md`.

## Reproducing with the public CLI

The release CLI provides the building blocks used by these runs:

```bash
cassia annotate --help
cassia boost run --help
cassia judge --help
```

For a new comparison, keep the case manifest, annotation prompt version,
backend/model, reasoning effort, marker ranking, and Judge configuration fixed
across arms. `cassia judge` evaluates broad lineage and rank-1 subtype for core
credit; lower-ranked subtypes remain diagnostic.

The checked-in manifests are sufficient for one-shot tests. Reproducing active
full-marker queries also requires the original target-vs-rest marker tables,
which are excluded because they are hundreds of megabytes. The summary records
include immutable selection and raw-artifact hashes so local source data can be
verified without publishing it.

## Interpretation limits

- All reported LLM runs are subject to model nondeterminism.
- Judge-relative scores should only be compared within the same locked Judge
  protocol and model configuration.
- Subscription-credit or API-equivalent costs are estimates, not necessarily
  the user's incremental bill.
- Exact labels can disagree with transcriptomic evidence when ground truth is
  protein-, spatial-, or developmental-stage-defined.
