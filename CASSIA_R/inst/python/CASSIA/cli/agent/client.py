"""TCP client for the CASSIA agent daemon.

Each `cassia agent <cmd>` invocation creates a fresh client, sends one JSON
request line, reads one JSON response line, disconnects.
"""

from __future__ import annotations

import json
import socket
import uuid
from pathlib import Path
from typing import Any, Dict

from .state import require_session


class DaemonError(RuntimeError):
    """Raised when the daemon returns ok=false."""


def call(workdir: Path, op: str, *, timeout: float = 1200.0, **args: Any) -> Dict[str, Any]:
    session = require_session(workdir)
    port = int(session["port"])
    req = {"id": uuid.uuid4().hex[:8], "op": op, "args": args}
    line = (json.dumps(req, ensure_ascii=False) + "\n").encode("utf-8")
    with socket.create_connection(("127.0.0.1", port), timeout=10.0) as s:
        s.settimeout(timeout)
        s.sendall(line)
        # Read until we see a newline (one response line).
        buf = bytearray()
        while True:
            chunk = s.recv(65536)
            if not chunk:
                break
            buf.extend(chunk)
            if b"\n" in buf:
                break
    raw = bytes(buf).decode("utf-8", errors="replace").strip()
    if not raw:
        raise DaemonError(f"empty response from daemon (op={op})")
    try:
        resp = json.loads(raw.splitlines()[0])
    except json.JSONDecodeError as e:
        raise DaemonError(f"bad json from daemon: {e}: {raw[:200]!r}")
    if not resp.get("ok"):
        raise DaemonError(resp.get("error", "<no error message>"))
    return resp


def call_raw(workdir: Path, op: str, **args: Any) -> Dict[str, Any]:
    """Like call() but returns the raw envelope (ok / error) without raising."""
    try:
        return call(workdir, op, **args)
    except DaemonError as e:
        return {"ok": False, "error": str(e)}
