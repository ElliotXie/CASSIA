"""Backend discovery and shell-agent execution for the CASSIA CLI."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence


API_BACKENDS = {
    "openrouter": {
        "description": "Use the existing CASSIA OpenRouter API provider.",
        "env": "OPENROUTER_API_KEY",
    },
    "openai": {
        "description": "Use the existing CASSIA OpenAI API provider.",
        "env": "OPENAI_API_KEY",
    },
    "anthropic": {
        "description": "Use the existing CASSIA Anthropic API provider.",
        "env": "ANTHROPIC_API_KEY",
    },
}


@dataclass(frozen=True)
class AgentBackend:
    """A local agent CLI backend that can be called non-interactively."""

    name: str
    executable: str
    argv_template: Sequence[str]
    description: str
    docs_url: str = ""

    @property
    def available(self) -> bool:
        return shutil.which(self.executable) is not None

    @property
    def executable_path(self) -> Optional[str]:
        return shutil.which(self.executable)


AGENT_BACKENDS: Dict[str, AgentBackend] = {
    "claude-cli": AgentBackend(
        name="claude-cli",
        executable="claude",
        argv_template=("claude", "-p", "{prompt}", "--output-format", "text"),
        description="Run Claude Code in non-interactive print mode.",
        docs_url="https://docs.claude.com/en/docs/claude-code/cli-reference",
    ),
    "codex-cli": AgentBackend(
        name="codex-cli",
        executable="codex",
        argv_template=("codex", "exec", "--skip-git-repo-check", "--output-last-message", "{agent_output_file}", "{prompt}"),
        description="Run OpenAI Codex CLI in non-interactive exec mode.",
        docs_url="https://help.openai.com/en/articles/11096431-openai-codex-cli-getting-started",
    ),
    "cursor-agent": AgentBackend(
        name="cursor-agent",
        executable="cursor-agent",
        argv_template=("cursor-agent", "-p", "{prompt}", "--output-format", "text"),
        description="Run Cursor Agent in non-interactive print mode.",
        docs_url="https://docs.cursor.com/en/cli/overview",
    ),
}


def list_backends() -> Dict[str, Dict[str, object]]:
    """Return all known CASSIA CLI backends with availability metadata."""
    backends: Dict[str, Dict[str, object]] = {}

    for name, info in API_BACKENDS.items():
        env_name = info["env"]
        backends[name] = {
            "name": name,
            "kind": "api",
            "available": bool(os.environ.get(env_name)),
            "requires": env_name,
            "description": info["description"],
        }

    for name, info in AGENT_BACKENDS.items():
        backends[name] = {
            "name": name,
            "kind": "agent-cli",
            "available": info.available,
            "executable": info.executable,
            "path": info.executable_path,
            "description": info.description,
            "docs_url": info.docs_url,
        }

    backends["shell"] = {
        "name": "shell",
        "kind": "agent-cli",
        "available": True,
        "executable": None,
        "path": None,
        "description": "Run a user-provided command template against each CASSIA prompt.",
    }
    return backends


def is_api_backend(name: str) -> bool:
    """Return whether a backend name should be routed to CASSIA's API provider path."""
    return name in API_BACKENDS or name.startswith("http://") or name.startswith("https://")


def is_agent_backend(name: str) -> bool:
    """Return whether a backend name is a local agent CLI backend."""
    return name in AGENT_BACKENDS or name == "shell"


def render_agent_argv(
    backend: AgentBackend,
    prompt: str,
    prompt_file: Path,
    context: Mapping[str, object],
) -> Sequence[str]:
    """Render a built-in backend argv template without invoking a shell."""
    values = {
        "prompt": prompt,
        "prompt_file": str(prompt_file),
        "input": str(context.get("input", "")),
        "out": str(context.get("out", "")),
        "cluster": str(context.get("cluster", "")),
        "agent_output_file": str(context.get("agent_output_file", "")),
    }
    return [part.format(**values) for part in backend.argv_template]


def render_shell_command(
    command_template: str,
    prompt: str,
    prompt_file: Path,
    context: Mapping[str, object],
) -> str:
    """Render a user command template with shell-quoted placeholders."""
    raw_values = {
        "prompt": prompt,
        "prompt_file": str(prompt_file),
        "input": str(context.get("input", "")),
        "out": str(context.get("out", "")),
        "cluster": str(context.get("cluster", "")),
        "agent_output_file": str(context.get("agent_output_file", "")),
    }
    quoted_values = {key: shlex.quote(value) for key, value in raw_values.items()}
    quoted_values.update({f"{key}_raw": value for key, value in raw_values.items()})
    return command_template.format(**quoted_values)


class AgentCLIBackend:
    """Execute CASSIA prompts through Claude Code, Codex CLI, Cursor Agent, or shell."""

    def __init__(
        self,
        name: str,
        command_template: Optional[str] = None,
        timeout_seconds: int = 900,
    ) -> None:
        if not is_agent_backend(name):
            raise ValueError(f"Unknown agent CLI backend: {name}")
        if name == "shell" and not command_template:
            raise ValueError("--command-template is required when --backend shell is used")
        self.name = name
        self.command_template = command_template
        self.timeout_seconds = timeout_seconds

    def run(
        self,
        prompt: str,
        prompt_file: Path,
        cwd: Path,
        context: Mapping[str, object],
    ) -> str:
        """Run one prompt and return stdout."""
        if self.command_template:
            command = render_shell_command(self.command_template, prompt, prompt_file, context)
            completed = subprocess.run(
                command,
                cwd=str(cwd),
                shell=True,
                text=True,
                capture_output=True,
                timeout=self.timeout_seconds,
            )
            display_command = command
        else:
            backend = AGENT_BACKENDS[self.name]
            if not backend.available:
                raise RuntimeError(
                    f"Backend '{self.name}' needs executable '{backend.executable}', "
                    "but it was not found on PATH."
                )
            argv = render_agent_argv(backend, prompt, prompt_file, context)
            completed = subprocess.run(
                argv,
                cwd=str(cwd),
                text=True,
                capture_output=True,
                timeout=self.timeout_seconds,
            )
            display_command = " ".join(shlex.quote(part) for part in argv)

        if completed.returncode != 0:
            stderr = completed.stderr.strip()
            stdout = completed.stdout.strip()
            detail = stderr or stdout or "no output"
            raise RuntimeError(
                f"Agent command failed with exit code {completed.returncode}: "
                f"{display_command}\n{detail}"
            )

        output_file = context.get("agent_output_file")
        if output_file:
            output_path = Path(str(output_file))
            if output_path.exists() and output_path.read_text(encoding="utf-8").strip():
                return output_path.read_text(encoding="utf-8").strip()

        return completed.stdout.strip()
