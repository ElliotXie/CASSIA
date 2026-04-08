"""
Progress tracking utilities for CASSIA batch operations.

This module provides thread-safe progress tracking with visual progress bars
for batch processing operations.

Rendering modes (auto-detected per environment):

- "rich":   Native Python terminal with UTF-8 stdout. Uses Braille spinner
            (⠋⠙⠹...), block progress bar (█░), and ANSI cursor control
            for in-place updates. The default for direct Python users.
- "ascii":  Non-UTF-8 or non-TTY stdout — most importantly, R via reticulate,
            which wraps Python stdout in an ASCII-encoded stream that crashes
            with UnicodeEncodeError on the Braille characters. Uses an ASCII
            spinner (|/-\\), ASCII bar (#-), no ANSI codes, and prints each
            update on a new line instead of overwriting in place.
- "notebook": Jupyter / Colab. Uses IPython.display.clear_output for in-place
            updates and emits Braille (notebooks render UTF-8 fine).

The previous implementation only branched on notebook-vs-terminal and assumed
terminal stdout was always UTF-8. Under R/reticulate that assumption is false
and a UnicodeEncodeError on \\u2819 (Braille spinner) crashed `tracker.finish()`,
which in turn aborted `runCASSIA_batch` *after* every cluster had succeeded —
making the whole batch appear to fail.
"""

import sys
import time
import threading
import atexit


def _is_notebook():
    """Detect if running in a Jupyter notebook or Google Colab environment."""
    try:
        from IPython import get_ipython
        shell = get_ipython()
        if shell is None:
            return False
        shell_name = shell.__class__.__name__
        if shell_name == 'ZMQInteractiveShell':
            return True  # Jupyter notebook or qtconsole
        elif shell_name == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        elif 'google.colab' in str(shell):
            return True  # Google Colab
        else:
            # Check for Colab specifically
            try:
                import google.colab
                return True
            except ImportError:
                pass
            return False
    except (NameError, ImportError):
        return False


def _stdout_supports_unicode():
    """
    Return True iff sys.stdout can encode the Braille spinner + block-bar
    characters used in 'rich' mode.

    Reticulate (R) wraps Python stdout in a stream whose encoding is the
    system locale, which is typically ASCII or Latin-1 on macOS/Windows
    R installations, and that stream raises UnicodeEncodeError on \\u2819.
    Detecting this *before* writing avoids the crash in `_render()`.
    """
    enc = getattr(sys.stdout, 'encoding', None)
    if not enc:
        return False
    try:
        '⠋█░'.encode(enc)
        return True
    except (UnicodeEncodeError, LookupError):
        return False


def _stdout_is_tty():
    """
    Return True iff sys.stdout is an interactive terminal.

    ANSI cursor control codes (\\033[K, \\033[?25l, etc.) only work on TTYs.
    Reticulate-wrapped stdout is not a TTY, so we must avoid those codes
    even if encoding happened to be UTF-8.
    """
    try:
        return bool(sys.stdout.isatty())
    except Exception:
        return False


class BatchProgressTracker:
    """Thread-safe progress tracker for batch processing with visual progress bar."""

    # Spinner animation frames — Unicode (Braille) for rich mode
    SPINNER_FRAMES = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
    # ASCII fallback spinner for non-UTF-8 stdout (e.g. R via reticulate)
    SPINNER_FRAMES_ASCII = ['|', '/', '-', '\\']

    def __init__(self, total, bar_width=40, refresh_rate=0.1, title="CASSIA Batch Analysis"):
        self.total = total
        self.completed = 0
        self.in_progress = set()
        self.lock = threading.Lock()
        self.bar_width = bar_width
        self.title = title
        self._lines_printed = 0
        self._spinner_idx = 0
        self._running = True
        self._refresh_rate = refresh_rate
        self._is_notebook = _is_notebook()

        # Decide rendering mode. Order matters: notebook first (its stdout
        # may not be a TTY but does support UTF-8), then check whether the
        # terminal stdout can both encode Unicode and accept ANSI codes.
        if self._is_notebook:
            self._mode = 'notebook'
        elif _stdout_supports_unicode() and _stdout_is_tty():
            self._mode = 'rich'
        else:
            # Catch-all: R/reticulate, redirected stdout, ASCII locales,
            # CI logs, etc. Safe everywhere; never raises UnicodeEncodeError.
            self._mode = 'ascii'

        # For notebook environments, use slower refresh to reduce flicker
        if self._is_notebook:
            self._refresh_rate = max(refresh_rate, 0.5)  # At least 0.5s in notebooks
            self._last_render_time = 0
        elif self._mode == 'ascii':
            # In ascii mode every render appends a new line (no in-place
            # update), so refresh slowly to avoid flooding the log.
            self._refresh_rate = max(refresh_rate, 1.0)

        # Start background thread for continuous spinner animation
        self._animation_thread = threading.Thread(target=self._animate, daemon=True)
        self._animation_thread.start()

        # Hide cursor during animation to prevent flashing.
        # Only safe to emit ANSI in rich mode (real TTY + UTF-8).
        if self._mode == 'rich':
            sys.stdout.write('\033[?25l')
            sys.stdout.flush()

        # Register cleanup to ensure cursor is restored on unexpected exit
        atexit.register(self._restore_cursor)

    def _animate(self):
        """Background thread that continuously updates the spinner."""
        while self._running:
            time.sleep(self._refresh_rate)
            with self.lock:
                if self._running and len(self.in_progress) > 0:
                    self._render()

    def _restore_cursor(self):
        """Restore cursor visibility. Called on exit or finish."""
        # Only emit the ANSI cursor-show sequence in rich mode; in ascii
        # mode the sequence would be printed verbatim, and in notebook
        # mode IPython handles its own cursor.
        if self._mode == 'rich':
            try:
                sys.stdout.write('\033[?25h')
                sys.stdout.flush()
            except Exception:
                pass

    def start_task(self, name):
        """Mark a task as started/in-progress."""
        with self.lock:
            self.in_progress.add(name)
            self._render()

    def complete_task(self, name):
        """Mark a task as completed."""
        with self.lock:
            self.in_progress.discard(name)
            self.completed += 1
            self._render()

    def _render(self):
        """Render the progress display."""
        # Calculate progress percentage
        pct = self.completed / self.total if self.total > 0 else 0
        filled = int(self.bar_width * pct)

        # Calculate counts
        processing = len(self.in_progress)
        pending = self.total - self.completed - processing

        # Pick characters and spinner frames per mode. ASCII mode uses
        # only 7-bit characters so it's safe under any stdout encoding,
        # most importantly the ASCII-wrapped stream that R/reticulate
        # exposes.
        if self._mode == 'ascii':
            bar = '#' * filled + '-' * (self.bar_width - filled)
            done_marker = 'OK'
            idle_marker = '.'
            spinner_frames = self.SPINNER_FRAMES_ASCII
        else:
            bar = '█' * filled + '░' * (self.bar_width - filled)
            done_marker = '✓'
            idle_marker = '○'
            spinner_frames = self.SPINNER_FRAMES

        # Get spinner character (only animate when processing)
        if processing > 0:
            spinner = spinner_frames[self._spinner_idx % len(spinner_frames)]
            self._spinner_idx += 1
        else:
            spinner = done_marker if self.completed == self.total else idle_marker

        # Truncate active task names if too many
        active_names = list(self.in_progress)[:3]
        active_str = ', '.join(str(name) for name in active_names)
        if len(self.in_progress) > 3:
            active_str += f', ... (+{len(self.in_progress)-3} more)'

        # Build display lines
        lines = [
            f"{self.title} {spinner}",
            f"[{bar}] {pct*100:.0f}%",
            f"Completed: {self.completed} | Processing: {processing} | Pending: {pending}",
            f"Active: {active_str if active_str else 'None'}"
        ]

        if self._mode == 'notebook':
            # In notebook environments, use clear_output to update in place.
            # Throttle updates to reduce flicker.
            current_time = time.time()
            if current_time - self._last_render_time < 0.3 and self.completed < self.total:
                return  # Skip this render to reduce flicker
            self._last_render_time = current_time

            try:
                from IPython.display import clear_output
                clear_output(wait=True)
            except ImportError:
                pass

            for line in lines:
                print(line)

        elif self._mode == 'rich':
            # Real terminal: use ANSI escape codes for in-place updates.
            if self._lines_printed > 0:
                sys.stdout.write(f'\033[{self._lines_printed}A')

            for line in lines:
                sys.stdout.write(f'\033[K{line}\n')

            sys.stdout.flush()

        else:
            # ASCII fallback: no ANSI, no Unicode. Print one compact summary
            # line per render so logs (and reticulate's wrapped stdout)
            # remain readable. Render is throttled by _refresh_rate to
            # avoid flooding.
            try:
                summary = (
                    f"[{self.title}] {spinner} {pct*100:3.0f}% "
                    f"({self.completed}/{self.total}) "
                    f"processing={processing} pending={pending}"
                )
                if active_str:
                    summary += f" active={active_str}"
                print(summary, flush=True)
            except Exception:
                # Last-resort safety net: never let a render error abort
                # the surrounding batch operation. (This is what bit us:
                # a UnicodeEncodeError in _render aborted runCASSIA_batch
                # *after* every cluster had already succeeded.)
                pass

        self._lines_printed = len(lines)

    def finish(self):
        """Finalize the progress display."""
        # Stop the animation thread
        self._running = False
        self._animation_thread.join(timeout=0.5)

        with self.lock:
            # Force final render regardless of throttling
            if self._is_notebook:
                self._last_render_time = 0  # Reset throttle for final render
            # Final render to show 100% with checkmark
            self._render()
            # Add blank line after completion (only in terminal)
            if not self._is_notebook:
                sys.stdout.write('\033[K\n')
                sys.stdout.flush()
            else:
                print()  # Simple newline in notebooks

        # Restore cursor visibility
        self._restore_cursor()

        # Unregister atexit handler since we've cleaned up normally
        try:
            atexit.unregister(self._restore_cursor)
        except Exception:
            pass  # Ignore if already unregistered
