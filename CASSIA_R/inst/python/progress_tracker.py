"""
Progress tracking utilities for CASSIA batch operations.

This module provides thread-safe progress tracking with visual progress bars
for batch processing operations.
"""

import sys
import time
import threading
import atexit


class BatchProgressTracker:
    """Thread-safe progress tracker for batch processing with visual progress bar."""

    # Spinner animation frames
    SPINNER_FRAMES = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']

    def __init__(self, total, bar_width=40, refresh_rate=0.1):
        self.total = total
        self.completed = 0
        self.in_progress = set()
        self.lock = threading.Lock()
        self.bar_width = bar_width
        self._lines_printed = 0
        self._spinner_idx = 0
        self._running = True
        self._refresh_rate = refresh_rate

        # Start background thread for continuous spinner animation
        self._animation_thread = threading.Thread(target=self._animate, daemon=True)
        self._animation_thread.start()

        # Hide cursor during animation to prevent flashing
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
        sys.stdout.write('\033[?25h')
        sys.stdout.flush()

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
        """Render the progress display, updating in place."""
        # Calculate progress percentage
        pct = self.completed / self.total if self.total > 0 else 0
        filled = int(self.bar_width * pct)
        bar = '█' * filled + '░' * (self.bar_width - filled)

        # Calculate counts
        processing = len(self.in_progress)
        pending = self.total - self.completed - processing

        # Get spinner character (only animate when processing)
        if processing > 0:
            spinner = self.SPINNER_FRAMES[self._spinner_idx % len(self.SPINNER_FRAMES)]
            self._spinner_idx += 1
        else:
            spinner = '✓' if self.completed == self.total else '○'

        # Truncate active task names if too many
        active_names = list(self.in_progress)[:3]
        active_str = ', '.join(str(name) for name in active_names)
        if len(self.in_progress) > 3:
            active_str += f', ... (+{len(self.in_progress)-3} more)'

        # Build display lines
        lines = [
            f"CASSIA Batch Analysis {spinner}",
            f"[{bar}] {pct*100:.0f}%",
            f"Completed: {self.completed} | Processing: {processing} | Pending: {pending}",
            f"Active: {active_str if active_str else 'None'}"
        ]

        # Move cursor up to overwrite previous output
        if self._lines_printed > 0:
            sys.stdout.write(f'\033[{self._lines_printed}A')

        # Print each line, clearing to end of line
        for line in lines:
            sys.stdout.write(f'\033[K{line}\n')

        sys.stdout.flush()
        self._lines_printed = len(lines)

    def finish(self):
        """Finalize the progress display."""
        # Stop the animation thread
        self._running = False
        self._animation_thread.join(timeout=0.5)

        with self.lock:
            # Final render to show 100% with checkmark
            self._render()
            # Add blank line after completion
            sys.stdout.write('\033[K\n')
            sys.stdout.flush()

        # Restore cursor visibility
        self._restore_cursor()

        # Unregister atexit handler since we've cleaned up normally
        try:
            atexit.unregister(self._restore_cursor)
        except Exception:
            pass  # Ignore if already unregistered
