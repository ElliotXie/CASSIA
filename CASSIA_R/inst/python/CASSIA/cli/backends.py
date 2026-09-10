"""Backend discovery and shell-agent execution for the CASSIA CLI."""

from __future__ import annotations

import json
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
    # When True, the prompt is piped to the agent on stdin instead of being
    # substituted into argv via the {prompt} placeholder. Required for agents
    # (notably cursor-agent) that treat argv-passed prompts as a partial system
    # message and only execute when the user message arrives on stdin.
    prompt_via_stdin: bool = False

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
        # JSON preserves Claude's result text and provider-reported token/cost
        # metadata. AgentCLIBackend unwraps the result before returning it.
        argv_template=("claude", "-p", "{prompt}", "--output-format", "json"),
        description="Run Claude Code in non-interactive print mode.",
        docs_url="https://docs.claude.com/en/docs/claude-code/cli-reference",
    ),
    "codex-cli": AgentBackend(
        name="codex-cli",
        executable="codex",
        argv_template=("codex", "exec", "--skip-git-repo-check", "--json", "--output-last-message", "{agent_output_file}", "{prompt}"),
        description="Run OpenAI Codex CLI in non-interactive exec mode.",
        docs_url="https://help.openai.com/en/articles/11096431-openai-codex-cli-getting-started",
    ),
    "cursor-agent": AgentBackend(
        name="cursor-agent",
        executable="cursor-agent",
        # Prompt is piped via stdin (see prompt_via_stdin). cursor-agent ignores
        # the task in long argv prompts and falls back to a conversational reply;
        # stdin makes it process the prompt as the user message.
        argv_template=(
            "cursor-agent",
            "-p",
            "--output-format", "text",
            "--trust",
            "--force",
        ),
        description="Run Cursor Agent in non-interactive print mode (prompt via stdin).",
        docs_url="https://docs.cursor.com/en/cli/overview",
        prompt_via_stdin=True,
    ),
    "opencode": AgentBackend(
        name="opencode",
        executable="opencode",
        argv_template=("opencode", "run", "--pure", "{prompt}"),
        description="Run OpenCode in non-interactive mode without external plugins.",
        docs_url="https://opencode.ai/docs/cli/",
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
    argv = [part.format(**values) for part in backend.argv_template]
    # On Windows, subprocess.run with shell=False does not consult PATHEXT, so a
    # bare command name like ``cursor-agent`` won't resolve to ``cursor-agent.cmd``.
    # Resolve argv[0] via shutil.which so .cmd/.bat/.exe shims are picked up.
    if argv:
        resolved = shutil.which(argv[0])
        if resolved:
            argv[0] = resolved
    return argv


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

    # Backends whose CLI accepts an explicit ``--model <id>`` flag.
    _BACKENDS_SUPPORTING_MODEL_FLAG = frozenset(
        {"cursor-agent", "codex-cli", "claude-cli", "opencode"}
    )
    _BACKENDS_SUPPORTING_EFFORT_FLAG = frozenset({"codex-cli", "claude-cli"})
    _CODEX_REASONING_EFFORTS = frozenset({"minimal", "low", "medium", "high", "xhigh", "max"})
    _CLAUDE_REASONING_EFFORTS = frozenset({"low", "medium", "high", "xhigh", "max"})
    _CODEX_SANDBOX_MODES = frozenset({"read-only", "workspace-write", "danger-full-access"})

    def __init__(
        self,
        name: str,
        command_template: Optional[str] = None,
        timeout_seconds: int = 900,
        model: Optional[str] = None,
        agent_mode: Optional[str] = None,
        reasoning_effort: Optional[str] = None,
        sandbox_mode: Optional[str] = None,
    ) -> None:
        if not is_agent_backend(name):
            raise ValueError(f"Unknown agent CLI backend: {name}")
        if name == "shell" and not command_template:
            raise ValueError("--command-template is required when --backend shell is used")
        if reasoning_effort and name not in self._BACKENDS_SUPPORTING_EFFORT_FLAG:
            raise ValueError("reasoning_effort is supported only by codex-cli and claude-cli")
        if reasoning_effort and name == "codex-cli" and reasoning_effort not in self._CODEX_REASONING_EFFORTS:
            allowed = ", ".join(sorted(self._CODEX_REASONING_EFFORTS))
            raise ValueError(f"Unsupported Codex reasoning effort: {reasoning_effort}. Choose: {allowed}")
        if reasoning_effort and name == "claude-cli" and reasoning_effort not in self._CLAUDE_REASONING_EFFORTS:
            allowed = ", ".join(sorted(self._CLAUDE_REASONING_EFFORTS))
            raise ValueError(f"Unsupported Claude reasoning effort: {reasoning_effort}. Choose: {allowed}")
        if sandbox_mode and name != "codex-cli":
            raise ValueError("sandbox_mode is supported only by codex-cli")
        if sandbox_mode and sandbox_mode not in self._CODEX_SANDBOX_MODES:
            allowed = ", ".join(sorted(self._CODEX_SANDBOX_MODES))
            raise ValueError(f"Unsupported Codex sandbox mode: {sandbox_mode}. Choose: {allowed}")
        self.name = name
        self.command_template = command_template
        self.timeout_seconds = timeout_seconds
        self.model = model
        self.agent_mode = agent_mode
        self.reasoning_effort = reasoning_effort
        self.sandbox_mode = sandbox_mode
        self.last_run_metadata: Dict[str, object] = {}

    def run(
        self,
        prompt: str,
        prompt_file: Path,
        cwd: Path,
        context: Mapping[str, object],
    ) -> str:
        """Run one prompt and return stdout."""
        self.last_run_metadata = {}
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
            argv = list(render_agent_argv(backend, prompt, prompt_file, context))
            if self.name == "cursor-agent" and self.agent_mode:
                argv = [part for part in argv if part != "--force"]
                argv += ["--mode", self.agent_mode]
            if self.model and self.name in self._BACKENDS_SUPPORTING_MODEL_FLAG:
                argv += ["--model", self.model]
            if self.name == "codex-cli" and self.reasoning_effort:
                # Codex reads model_reasoning_effort from the normal config layer.
                # Passing a per-process override keeps benchmark arms reproducible
                # without changing the user's ~/.codex/config.toml.
                argv += ["--config", f'model_reasoning_effort="{self.reasoning_effort}"']
            if self.name == "codex-cli" and self.sandbox_mode:
                argv += ["--sandbox", self.sandbox_mode]
            if self.name == "claude-cli" and self.reasoning_effort:
                argv += ["--effort", self.reasoning_effort]
            stdin_input = prompt if backend.prompt_via_stdin else None
            completed = subprocess.run(
                argv,
                cwd=str(cwd),
                text=True,
                input=stdin_input,
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

        if self.name == "codex-cli" and not self.command_template:
            for line in completed.stdout.splitlines():
                try:
                    event = json.loads(line)
                except (json.JSONDecodeError, TypeError):
                    continue
                if event.get("type") == "thread.started":
                    self.last_run_metadata["thread_id"] = event.get("thread_id")
                elif event.get("type") == "turn.completed":
                    usage = event.get("usage") or {}
                    self.last_run_metadata["usage"] = {
                        "input_tokens": int(usage.get("input_tokens", 0) or 0),
                        "cached_input_tokens": int(usage.get("cached_input_tokens", 0) or 0),
                        "output_tokens": int(usage.get("output_tokens", 0) or 0),
                        "reasoning_output_tokens": int(usage.get("reasoning_output_tokens", 0) or 0),
                    }

        if self.name == "claude-cli" and not self.command_template:
            try:
                payload = json.loads(completed.stdout)
            except (json.JSONDecodeError, TypeError) as exc:
                raise RuntimeError("Claude CLI returned invalid JSON output") from exc
            usage = payload.get("usage") or {}
            self.last_run_metadata = {
                "thread_id": payload.get("session_id"),
                "duration_ms": int(payload.get("duration_ms", 0) or 0),
                "cost_usd": float(payload.get("total_cost_usd", 0) or 0),
                "usage": {
                    "input_tokens": int(usage.get("input_tokens", 0) or 0),
                    "cached_input_tokens": int(usage.get("cache_read_input_tokens", 0) or 0),
                    "cache_creation_input_tokens": int(usage.get("cache_creation_input_tokens", 0) or 0),
                    "output_tokens": int(usage.get("output_tokens", 0) or 0),
                    "reasoning_output_tokens": 0,
                },
            }
            return str(payload.get("result", "")).strip()

        output_file = context.get("agent_output_file")
        if output_file:
            output_path = Path(str(output_file))
            if output_path.exists() and output_path.read_text(encoding="utf-8").strip():
                return output_path.read_text(encoding="utf-8").strip()

        return completed.stdout.strip()
