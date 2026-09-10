# Integrated agent development architecture

This package is in development. Nothing in this directory publishes a package,
creates a release, pushes Git state, or uploads a Seurat object.

## Python boundaries

| Module | Owns |
| --- | --- |
| `automation_cli.py` | `agent auto` / `agent compare` argparse and thin handlers |
| `automation.py` | End-to-end orchestration and recoverable run manifest lifecycle |
| `policy.py` | Fixed, conservative, and adaptive topology contracts |
| `prompting.py` | Coding-agent prompt construction |
| `experiment.py` | Audit parsing and deterministic saved-run comparison |
| `direct.py` | Locked, versioned, no-daemon transaction store |
| `commands.py` | Interactive scientific CLI commands |
| `backends.py` | External coding-agent process adapter |

`automation.py` re-exports the original prompt/report symbols so code written
against the first development iteration remains compatible.

## R boundaries

| Module | Owns |
| --- | --- |
| `transaction_worker.R` | Default one-request/one-process entrypoint |
| `daemon.R` | Explicit optional localhost daemon entrypoint |
| `runtime.R` | Shared Seurat operations and JSON dispatch |
| `policy.R` | Topology budget and split-fragmentation validation |

The default path is:

```text
coding agent
  → cassia agent command
  → direct.py lock + request
  → transaction_worker.R
  → runtime.R operation + policy.R guard
  → new immutable RDS checkpoint
  → atomic state.json pointer update
```

The daemon does not participate unless the user explicitly passes `--daemon`.

## Persisted run contract

`.cassia/state.json` selects the current checkpoint. Successful state-changing
commands append a file in `.cassia/state_versions/` and a matching
`transaction_commit` record in `.cassia/audit.jsonl`. Failed or interrupted
transactions cannot advance the pointer. Baseline/final memberships, recovery
checkpoints, prompts, model settings, QA, and deterministic reports stay inside
the chosen development run directory.
