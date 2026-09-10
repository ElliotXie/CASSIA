# CASSIA 1.4.0 release notes

CASSIA 1.4.0 makes the command-line interface a supported package entry point
and packages the benchmark-tested agent workflows that were previously kept in
the research workspace.

## Highlights

- `cassia annotate` supports API providers and local Codex, Claude Code,
  Cursor, OpenCode, or custom-shell agent backends.
- Three explicit annotation modes: `one-shot`, `validated`, and `fused-boost`.
- `v2-compact` is the Fused Boost release default. Historical `v2` remains
  available for exact legacy-prompt reproduction.
- Stable `cassia.annotation.v1` result normalization across standard and Fused
  Boost outputs.
- `cassia judge` implements the blinded `stable-judge-v1.1` protocol with
  rank-1 subtype determining core credit.
- `cassia guide`, nested `cassia help`, `doctor`, `validate`, `examples`,
  `report`, and `resume` cover setup, preflight, and run recovery.
- `cassia boost`, `subcluster`, and deterministic `consensus` workflows are
  available from the installed executable.
- `cassia agent` provides an optional long-lived R/Seurat bridge with marker
  queries, cluster operations, labels, QA, finalization, and audit history.

## Benchmark evidence

The release snapshot under `Benchmark/cli_release/` contains fixed 50-, 100-,
and 392-case manifests and compact result records.

- Boost-100: one-shot scored 50/100 exact and review Boost scored 71/100 exact
  on an intentionally hard, boost-enriched diagnostic set. There were 26 strict
  rescues and five strict regressions, supporting selective rather than
  unconditional Boost.
- Full392 prompt ablation: legacy Fused Boost v2 and the compact treatment both
  scored 296/392 exact. Correct discordance was 19 versus 19 (exact McNemar
  p=1.0); estimated annotation cost was 16.3% lower for the compact arm.
- The 50-case cross-model card is included as a capability profile. It is a
  single-run, fixed-referee comparison and should not be read as a permanent
  leaderboard.

## Compatibility and packaging

- Python 3.9 through 3.13 remain supported.
- The wheel and source distribution include the CLI guide, agent system prompt,
  and R daemon resource.
- Existing Python APIs remain available; the CLI is additive.

## Release verification

Before upload, the release artifacts must pass:

1. the full CLI pytest suite;
2. wheel and source-distribution metadata/content inspection;
3. installation into a fresh virtual environment;
4. offline example, validate, one-shot/validated shell-agent, Fused Boost
   dry-run, Judge help, and packaged-resource smoke tests;
5. one low-cost live OpenRouter annotation smoke test;
6. one-cluster live annotation smoke tests through Codex CLI, Claude Code,
   Cursor Agent, and OpenCode.
