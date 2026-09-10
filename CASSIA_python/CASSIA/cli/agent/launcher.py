"""Launch (and shut down) the detached R daemon."""

from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Union

from .state import DAEMON_LOG, ensure_workdir, read_session, session_path


def _daemon_script() -> Path:
    here = Path(__file__).resolve().parent
    return here / "daemon.R"


def _pick_free_port() -> int:
    """Bind to port 0 to let the OS choose, then close. Race-prone by design
    (the port could be claimed by another process in the gap) but the window
    is microseconds; the daemon will fail fast if so and the user retries."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _is_daemon_alive(session: Dict) -> bool:
    """Cheap liveness probe: try a 1-line ping over TCP."""
    port = session.get("port")
    if not port:
        return False
    try:
        with socket.create_connection(("127.0.0.1", int(port)), timeout=1.0) as s:
            s.sendall(b'{"id":"probe","op":"ping","args":{}}\n')
            data = s.recv(4096)
            return b'"ok":true' in data
    except (OSError, socket.timeout):
        return False


def spawn(workdir_base: Union[Path, str], rscript: str = "Rscript") -> Dict:
    """Spawn a detached R daemon. Returns the session dict once the daemon
    has written session.json and answered a ping.

    Raises:
        FileNotFoundError: if Rscript or daemon.R is missing.
        RuntimeError: if the daemon did not come up within 30s, or if a live
            session already exists.
    """
    wd = ensure_workdir(workdir_base)
    existing = read_session(wd)
    if existing is not None and _is_daemon_alive(existing):
        raise RuntimeError(
            f"a daemon is already alive in {wd} (pid {existing.get('pid')}, "
            f"port {existing.get('port')}). Use `cassia agent close` first."
        )
    if existing is not None:
        # stale session.json — daemon dead. Clear it before respawning.
        try:
            session_path(wd).unlink()
        except OSError:
            pass

    script = _daemon_script()
    if not script.exists():
        raise FileNotFoundError(script)
    if shutil.which(rscript) is None:
        raise FileNotFoundError(f"could not find {rscript} on PATH")

    log_path = wd / DAEMON_LOG
    log_fh = log_path.open("ab")

    kwargs: Dict = {
        "stdin": subprocess.DEVNULL,
        "stdout": log_fh,
        "stderr": log_fh,
        "close_fds": True,
    }
    if sys.platform.startswith("win"):
        # DETACHED_PROCESS so the daemon survives parent exit; new process
        # group so Ctrl-C in the parent doesn't kill it.
        DETACHED = getattr(subprocess, "DETACHED_PROCESS", 0x00000008)
        CREATE_NEW_PROCESS_GROUP = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x00000200)
        kwargs["creationflags"] = DETACHED | CREATE_NEW_PROCESS_GROUP
    else:
        kwargs["start_new_session"] = True

    port = _pick_free_port()
    ready_path = wd / "daemon.ready"
    if ready_path.exists():
        try:
            ready_path.unlink()
        except OSError:
            pass

    cmd = [rscript, "--vanilla", str(script), str(wd), str(port)]
    proc = subprocess.Popen(cmd, **kwargs)

    deadline = time.time() + 30.0
    while time.time() < deadline:
        if proc.poll() is not None:
            return _spawn_failed(proc.returncode, log_path)
        if ready_path.exists():
            break
        time.sleep(0.1)
    else:
        if proc.poll() is None:
            try:
                proc.terminate()
            except OSError:
                pass
        raise RuntimeError(
            f"daemon failed to come up within 30s. Check {log_path}."
        )

    session = {
        "pid": proc.pid,
        "port": port,
        "workdir": str(wd),
        "started_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    session_path(wd).write_text(
        json.dumps(session, indent=2) + "\n", encoding="utf-8"
    )
    if not _is_daemon_alive(session):
        try:
            proc.terminate()
        except OSError:
            pass
        raise RuntimeError(
            f"daemon bound port {port} but failed ping. Check {log_path}."
        )
    return session


def _spawn_failed(returncode: int, log_path: Path) -> Dict:
    tail = ""
    try:
        with log_path.open("rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - 4000))
            tail = f.read().decode("utf-8", errors="replace")
    except OSError:
        pass
    raise RuntimeError(
        f"Rscript daemon exited with code {returncode} before becoming "
        f"ready. Log tail:\n{tail}"
    )


def shutdown(workdir: Path, timeout: float = 10.0) -> bool:
    """Tell the daemon to exit and wait for session.json to disappear."""
    from .client import call_raw

    session = read_session(workdir)
    if session is None:
        return False
    call_raw(workdir, "shutdown")
    deadline = time.time() + timeout
    ready_path = workdir / "daemon.ready"
    while time.time() < deadline:
        if not ready_path.exists():
            try:
                session_path(workdir).unlink()
            except OSError:
                pass
            return True
        time.sleep(0.1)
    # Best-effort: force-kill if pid is known.
    pid = int(session.get("pid", 0))
    if pid > 0:
        try:
            if sys.platform.startswith("win"):
                subprocess.run(["taskkill", "/F", "/PID", str(pid)], check=False,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                os.kill(pid, 9)
        except OSError:
            pass
    try:
        session_path(workdir).unlink()
    except OSError:
        pass
    return True
