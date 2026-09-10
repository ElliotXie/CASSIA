# CASSIA CLI Guide for Coding Agents

Use this guide when a user asks you to annotate single-cell RNA-seq marker
clusters with CASSIA through Codex CLI, Claude CLI, Cursor Agent, OpenCode, or
another installed agent backend.

## Operating rules

1. Work only with marker files and run directories the user places in scope.
2. Run `cassia doctor` and `cassia validate INPUT.csv` before annotation.
3. Use the user's existing agent subscription through the requested local
   backend. Do not request a CASSIA API key for `codex-cli`, `claude-cli`,
   `cursor-agent`, or `opencode`.
4. Do not edit CASSIA annotation, validator, or Fused Boost prompts.
5. Do not invent marker statistics. Fused Boost obtains additional evidence
   through CASSIA's marker-query loop.
6. Do not set `--max-genes-per-round` unless the user explicitly wants a cap.
7. Treat saved JSON as the annotation source of truth. HTML and Markdown are
   deterministic views and do not add an LLM call.
8. Never publish, upload, commit, or push results unless the user separately
   asks for that action.
9. For benchmark evaluation, use `cassia judge` and keep its protocol, judge
   model, reasoning effort, context, and inputs fixed across compared arms.

## Choose a mode

| User intent | Mode | Input |
| --- | --- | --- |
| Fast baseline or one-shot benchmark | `one-shot` | One marker list per cluster, or a differential-expression table |
| Faithful CASSIA annotation plus validator/revision loop | `validated` | Same as one-shot |
| Primary annotation with active full-marker evidence queries | `fused-boost` | Raw target-vs-rest differential-expression table and one `--cluster` |

Use `--prompt-version v1` when the user requests the original CASSIA-based
one-shot/validated prompt. `fused-boost` defaults to the benchmark-backed
`v2-compact` prompt. Use `--fused-prompt-version v2` only to reproduce the
legacy prompt with its duplicated step-by-step prefix. The legacy `--workflow`
option is an alias for `--mode`.

## Preflight

```bash
cassia --version
cassia doctor
cassia backends list
cassia help annotate
cassia validate INPUT.csv
```

If validation cannot infer columns, pass `--celltype-column` and
`--gene-column`. For a differential-expression table, retain useful statistics
such as `avg_log2FC`, `pct.1`, `pct.2`, and `p_val_adj`.

## Run one-shot annotation

```bash
cassia annotate \
  --input markers.csv \
  --backend codex-cli \
  --mode one-shot \
  --prompt-version v1 \
  --tissue brain \
  --species human \
  --out runs/brain_one_shot
```

Replace `codex-cli` with `claude-cli`, `cursor-agent`, or `opencode` when
requested. A model can be selected with `--model`, for example
`--model composer-2.5` for Cursor Agent or
`--model opencode/ling-3.0-flash-fin-free` for OpenCode when available in the
user's installation.

## Run validated annotation

```bash
cassia annotate \
  --input markers.csv \
  --backend cursor-agent \
  --model composer-2.5 \
  --mode validated \
  --prompt-version v1 \
  --validator-involvement v1 \
  --tissue brain \
  --species human \
  --out runs/brain_validated
```

Validated mode starts from the same annotation prompt and then runs the CASSIA
validator/revision loop. Do not compare it with one-shot unless both runs use
the same input clusters, tissue/species context, prompt version, backend/model,
and marker ranking configuration.

## Run Fused Boost primary annotation

```bash
cassia annotate \
  --input raw_findallmarkers.csv \
  --celltype-column cluster \
  --gene-column gene \
  --cluster 3 \
  --backend cursor-agent \
  --model composer-2.5 \
  --mode fused-boost \
  --tissue brain \
  --species human \
  --out runs/brain_cluster_3_fused
```

Fused Boost does not consume or trust a previous annotation. It must complete
at least one marker-query round before finalizing. Run one target cluster per
command. Use a unique output directory for every cluster and model/mode.

## Run integrated clustering + annotation on Seurat

`cassia agent auto` lets a coding agent investigate markers and make bounded
local merge/subcluster decisions before annotation. The default architecture
does **not** use a daemon: every command is a short-lived R transaction that
loads the current checkpoint and atomically commits a new one only on success.

```bash
cassia agent auto object.rds \
  --out runs/object_conservative \
  --strategy conservative \
  --backend codex-cli \
  --model gpt-6-astra \
  --reasoning-effort high
```

The strategies are intentionally distinct:

| Strategy | Topology behavior |
| --- | --- |
| `fixed` | Never changes the supplied clustering |
| `conservative` | At most two evidence-backed local edits; rejects children below 5 cells and splits above 3 children by default |
| `adaptive` | Allows the configured bounded edit budget and fragmentation limit |

For a mixed cluster, use `cassia agent subcluster UUID --auto-resolution`.
The CLI tries a small increasing resolution ladder and commits only the first
policy-valid split. One-child, undersized-child, and over-fragmented results
roll back without spending the topology budget. Exact baseline/final cell
memberships, immutable state versions, the audit log, model settings, and QA
outputs are retained in the run directory.

Compare completed experiments deterministically, without another model call:

```bash
cassia agent compare runs/fixed runs/conservative runs/adaptive \
  --out runs/comparison.md
```

A persistent R process is available only as an explicit accelerator via
`cassia agent init object.rds --daemon`; it is not required for correctness or
resume.

The development module boundaries and persisted-state contract are documented
in `cli/agent/ARCHITECTURE.md`.

## Read results

Standard one-shot/validated runs write:

- `results.json` — structured annotations by cluster
- `summary.csv` — compact table
- `report.md` and `report.html` — deterministic reports
- `run_manifest.json` — parameters, status, schema version, and output paths

Fused Boost runs write:

- `final.json` — structured final annotation
- `boost_manifest.json` — parameters, query history summary, schema version,
  status, and output paths
- `queries/*.csv` and `queries/*.json` — returned marker statistics
- `transcript.md` — prompt/response trace
- `report.md`, `summary.html`, and `summary_tags.txt` — deterministic reports;
  `summary.html` uses the dedicated Fused Boost primary-annotation template
  rather than the legacy revision-style Annotation Boost layout

The shared core schema is `cassia.annotation.v1`:

```json
{
  "schema_version": "cassia.annotation.v1",
  "annotation_mode": "standard or fused_boost",
  "cluster_id": "3",
  "main_cell_type": "...",
  "sub_cell_types": ["..."],
  "possible_mixed_cell_types": [],
  "confidence": "...",
  "evidence": "..."
}
```

Fused Boost keeps additional `final_*`, `ranked_sub_cell_types`,
`checked_genes`, `supporting_markers`, `refuting_markers`, `alternatives`, and
`recommended_next_steps` fields.

## Resume or rebuild reports

```bash
cassia resume RUN_DIR
cassia report RUN_DIR
```

Use `resume` for an interrupted standard agent-CLI run. Use `report` to rebuild
Markdown and HTML from saved artifacts. Report rebuilding never changes the
annotation and never calls an LLM.

## Evaluate benchmark outputs

Use the public stable Judge when the user provides held-out labels:

```bash
cassia judge \
  --truth truth.csv \
  --prediction baseline:runs/baseline/summary.csv \
  --prediction candidate:runs/candidate/summary.csv \
  --truth-marker-column marker_list \
  --backend codex-cli \
  --judge-model gpt-5.6-luna \
  --judge-reasoning-effort high \
  --out runs/judge_luna_high
```

The protocol scores `main_cell_type` for lineage and only rank-1 subtype for
the primary verdict. Rank 2/3 are diagnostic only. Never remove failed/missing
predictions, edit labels, or change the Judge model/effort between arms.

## Completion checklist

- Confirm the command exit status.
- Read the run manifest and structured JSON.
- Confirm the expected report files exist.
- Report failed clusters or missing final JSON explicitly.
- Give the user the run directory and per-cluster cost only when reliable
  usage/cost metadata is actually available. Agent CLIs do not expose this
  consistently; do not invent a measured cost.
- Do not claim success from an HTML file alone; the structured result and
  manifest status must both be valid.
