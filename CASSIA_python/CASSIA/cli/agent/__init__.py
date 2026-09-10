"""Transactional CASSIA CLI for integrated clustering and annotation.

The default path executes each command in a short-lived R process against an
immutable checkpoint. A detached R daemon is retained only as an explicit
``cassia agent init --daemon`` accelerator for unusually large objects.
"""

from __future__ import annotations

from . import client, launcher, state  # noqa: F401
