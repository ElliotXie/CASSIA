"""Agent state files under <workdir>/.

Layout:
    .cassia/
      state.json      # selected immutable checkpoint and revision pointer
      state_versions/ # one RDS per successful state-changing transaction
      audit.jsonl     # append-only operation and transaction log
      session.json    # only present when the optional daemon is running
      daemon.log      # optional daemon/worker diagnostics
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional, Union


WORKDIR_NAME = ".cassia"
SESSION_FILE = "session.json"
DAEMON_LOG = "daemon.log"
AUDIT_LOG = "audit.jsonl"


def workdir_for(base: Union[Path, str]) -> Path:
    base = Path(base)
    return base / WORKDIR_NAME if base.name != WORKDIR_NAME else base


def session_path(workdir: Path) -> Path:
    return workdir / SESSION_FILE


def read_session(workdir: Path) -> Optional[Dict[str, Any]]:
    path = session_path(workdir)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def ensure_workdir(base: Union[Path, str]) -> Path:
    wd = workdir_for(base)
    wd.mkdir(parents=True, exist_ok=True)
    return wd


def require_session(workdir: Path) -> Dict[str, Any]:
    s = read_session(workdir)
    if s is None:
        raise SystemExit(
            f"no active daemon at {workdir} — rerun `cassia agent init <rds> --daemon`"
        )
    return s
