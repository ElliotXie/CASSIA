"""Transactional, no-daemon execution for ``cassia agent`` commands.

Each operation runs in a short-lived R process.  State-changing operations
write a new immutable checkpoint and advance a small JSON pointer only after
the R transaction succeeds.  Interrupted commands therefore leave the last
complete state selected.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, Optional, Union

try:  # pragma: no cover - platform branch
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None  # type: ignore[assignment]
    import msvcrt

from .client import DaemonError
from .state import ensure_workdir


STATE_SCHEMA_VERSION = "cassia.agent.direct-state.v1"
STATE_MANIFEST = "state.json"
STATE_VERSIONS_DIR = "state_versions"
_PERSISTING_OPS = {
    "init",
    "restore",
    "policy",
    "markers",
    "preprocess",
    "merge",
    "subcluster",
    "label",
    "unlabel",
    "finalize",
}


def _operation_persists(op: str, op_args: Mapping[str, Any]) -> bool:
    if op != "markers":
        return op in _PERSISTING_OPS
    cluster = str(op_args.get("cluster") or "")
    versus = op_args.get("vs")
    # FindAllMarkers and one-vs-rest refresh the reusable marker cache.
    # Pairwise comparisons are read-only evidence and should not clone the RDS.
    return cluster.lower() == "all" or versus in (None, "", "rest")


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _append_audit(workdir: Path, value: Mapping[str, Any]) -> None:
    with (workdir / "audit.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(value, ensure_ascii=False, default=str) + "\n")


@contextmanager
def _state_lock(workdir: Path) -> Iterator[None]:
    """Serialize operations; OS releases this lock if the process crashes."""
    lock_path = workdir / "state.lock"
    with lock_path.open("a+", encoding="utf-8") as handle:
        if fcntl is not None:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        else:  # pragma: no cover - Windows
            handle.seek(0)
            if handle.read(1) == "":
                handle.write("0")
                handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            else:  # pragma: no cover - Windows
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)


def _manifest_path(workdir: Path) -> Path:
    return workdir / STATE_MANIFEST


def _read_manifest(workdir: Path) -> Optional[Dict[str, Any]]:
    path = _manifest_path(workdir)
    if not path.exists():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"invalid direct-state manifest: {path}") from exc
    if value.get("schema_version") != STATE_SCHEMA_VERSION:
        raise RuntimeError(f"unsupported direct-state manifest: {path}")
    return value


def has_state(workdir_base: Union[Path, str]) -> bool:
    return _manifest_path(ensure_workdir(workdir_base)).exists()


def state_status(workdir_base: Union[Path, str]) -> Dict[str, Any]:
    workdir = ensure_workdir(workdir_base)
    manifest = _read_manifest(workdir)
    if manifest is None:
        raise SystemExit(
            f"no persisted state at {workdir} — run `cassia agent init <rds>` first"
        )
    return manifest


def _worker_script() -> Path:
    return Path(__file__).resolve().with_name("transaction_worker.R")


def _run_worker(
    workdir: Path,
    *,
    op: str,
    op_args: Mapping[str, Any],
    input_state: Optional[Path],
    output_state: Optional[Path],
    timeout: float,
) -> Dict[str, Any]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise FileNotFoundError("could not find Rscript on PATH")
    request = {
        "id": str(uuid.uuid4()),
        "op": op,
        "args": dict(op_args),
        "input_state": str(input_state.resolve()) if input_state else None,
        "output_state": str(output_state.resolve()) if output_state else None,
    }
    request_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".json",
            prefix="cassia-agent-",
            dir=str(workdir),
            delete=False,
        ) as handle:
            json.dump(request, handle, ensure_ascii=False)
            handle.write("\n")
            request_path = Path(handle.name)
        completed = subprocess.run(
            [
                rscript,
                "--vanilla",
                str(_worker_script()),
                str(workdir),
                str(request_path),
            ],
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    finally:
        if request_path is not None:
            request_path.unlink(missing_ok=True)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip() or "no output"
        raise RuntimeError(f"CASSIA R transaction exited {completed.returncode}: {detail}")
    response = None
    for line in reversed(completed.stdout.splitlines()):
        try:
            candidate = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(candidate, dict) and "ok" in candidate:
            response = candidate
            break
    if response is None:
        detail = completed.stderr.strip()
        raise RuntimeError(
            f"CASSIA R transaction returned invalid JSON: {completed.stdout!r}\n{detail}"
        )
    if not response.get("ok"):
        raise DaemonError(str(response.get("error") or "unknown R transaction error"))
    return response


def call(
    workdir_base: Union[Path, str],
    op: str,
    *,
    timeout: float = 600.0,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Execute one CLI operation against versioned state without a daemon."""
    workdir = ensure_workdir(workdir_base)
    with _state_lock(workdir):
        manifest = _read_manifest(workdir)
        if op == "init":
            if manifest is not None:
                raise DaemonError(
                    f"persisted state already exists at {workdir}; use a new workdir"
                )
            rds_path = Path(str(kwargs.get("rds_path", ""))).expanduser().resolve()
            if not rds_path.exists():
                raise DaemonError(f"rds path not found: {rds_path}")
            manifest = {
                "schema_version": STATE_SCHEMA_VERSION,
                "created_at": _utc_now(),
                "updated_at": _utc_now(),
                "input": {
                    "path": str(rds_path),
                    "sha256": _sha256_file(rds_path),
                    "size_bytes": rds_path.stat().st_size,
                },
                "revision": 0,
                "current_state": None,
                "history": [],
            }
            input_state = None
        else:
            if manifest is None:
                raise SystemExit(
                    f"no persisted state at {workdir} — run `cassia agent init <rds>` first"
                )
            current = manifest.get("current_state")
            if not current:
                raise RuntimeError(f"direct-state manifest has no current checkpoint: {workdir}")
            input_state = Path(str(current))
            if not input_state.exists():
                raise RuntimeError(f"current direct-state checkpoint is missing: {input_state}")

        should_persist = _operation_persists(op, kwargs)
        revision = int(manifest.get("revision", 0)) + (1 if should_persist else 0)
        output_state = None
        if should_persist:
            versions = workdir / STATE_VERSIONS_DIR
            versions.mkdir(parents=True, exist_ok=True)
            suffix = uuid.uuid4().hex[:8]
            output_state = versions / f"{revision:06d}_{op}_{suffix}.rds"

        started = time.monotonic()
        try:
            response = _run_worker(
                workdir,
                op=op,
                op_args=kwargs,
                input_state=input_state,
                output_state=output_state,
                timeout=timeout,
            )
        except Exception:
            # The R process may have installed an output checkpoint before a
            # transport/parsing failure. It is not selected by state.json and
            # is safe to remove as an uncommitted transaction artifact.
            if output_state is not None:
                output_state.unlink(missing_ok=True)
            raise
        if output_state is not None:
            manifest["revision"] = revision
            manifest["current_state"] = str(output_state.resolve())
        manifest["updated_at"] = _utc_now()
        manifest.setdefault("history", []).append(
            {
                "revision": revision if should_persist else manifest.get("revision", 0),
                "op": op,
                "persisted": should_persist,
                "elapsed_s": round(time.monotonic() - started, 4),
                "completed_at": _utc_now(),
            }
        )
        _atomic_json(_manifest_path(workdir), manifest)
        if should_persist:
            _append_audit(
                workdir,
                {
                    "event": "transaction_commit",
                    "id": response.get("id"),
                    "op": op,
                    "revision": revision,
                    "state_path": str(output_state),
                    "ts": _utc_now(),
                },
            )
        return response


def copy_current_checkpoint(
    workdir_base: Union[Path, str], destination: Union[Path, str]
) -> Dict[str, Any]:
    """Copy the selected immutable state checkpoint to an explicit path."""
    workdir = ensure_workdir(workdir_base)
    with _state_lock(workdir):
        manifest = _read_manifest(workdir)
        if manifest is None or not manifest.get("current_state"):
            raise SystemExit(
                f"no persisted state at {workdir} — run `cassia agent init <rds>` first"
            )
        source = Path(str(manifest["current_state"]))
        target = Path(destination).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(f".{target.name}.{os.getpid()}.tmp")
        shutil.copyfile(source, temporary)
        temporary.replace(target)
        return {
            "schema_version": "cassia.agent.checkpoint.v1",
            "path": str(target),
            "sha256": _sha256_file(target),
            "revision": manifest.get("revision"),
        }


def restore_checkpoint(
    workdir_base: Union[Path, str], checkpoint: Union[Path, str], *, timeout: float = 1800.0
) -> Dict[str, Any]:
    """Validate an exported checkpoint and select a new immutable copy of it."""
    workdir = ensure_workdir(workdir_base)
    source = Path(checkpoint).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"checkpoint not found: {source}")
    with _state_lock(workdir):
        manifest = _read_manifest(workdir)
        if manifest is None:
            manifest = {
                "schema_version": STATE_SCHEMA_VERSION,
                "created_at": _utc_now(),
                "updated_at": _utc_now(),
                "input": {
                    "checkpoint_path": str(source),
                    "checkpoint_sha256": _sha256_file(source),
                },
                "revision": 0,
                "current_state": None,
                "history": [],
            }
        revision = int(manifest.get("revision", 0)) + 1
        versions = workdir / STATE_VERSIONS_DIR
        versions.mkdir(parents=True, exist_ok=True)
        output_state = versions / f"{revision:06d}_restore_{uuid.uuid4().hex[:8]}.rds"
        try:
            response = _run_worker(
                workdir,
                op="restore",
                op_args={"path": str(source)},
                input_state=None,
                output_state=output_state,
                timeout=timeout,
            )
        except Exception:
            output_state.unlink(missing_ok=True)
            raise
        manifest["revision"] = revision
        manifest["current_state"] = str(output_state.resolve())
        manifest["updated_at"] = _utc_now()
        manifest.setdefault("history", []).append(
            {
                "revision": revision,
                "op": "restore",
                "persisted": True,
                "source": str(source),
                "completed_at": _utc_now(),
            }
        )
        _atomic_json(_manifest_path(workdir), manifest)
        _append_audit(
            workdir,
            {
                "event": "transaction_commit",
                "id": response.get("id"),
                "op": "restore",
                "revision": revision,
                "state_path": str(output_state),
                "ts": _utc_now(),
            },
        )
        return response


__all__ = [
    "STATE_SCHEMA_VERSION",
    "call",
    "copy_current_checkpoint",
    "has_state",
    "restore_checkpoint",
    "state_status",
]
