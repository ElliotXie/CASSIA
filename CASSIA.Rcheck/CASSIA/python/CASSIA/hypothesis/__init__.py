# CASSIA Hypothesis Module
# Hypothesis generation

from .hypothesis_generation import generate_hypothesis, process_marker_file, run_multi_analysis
from .summarize_hypothesis_runs import summarize_runs

# Aliases for backward compatibility
runCASSIA_hypothesis = generate_hypothesis
run_hypothesis_analysis = run_multi_analysis
summarize_hypothesis_runs = summarize_runs

__all__ = [
    'generate_hypothesis',
    'process_marker_file',
    'run_multi_analysis',
    'summarize_runs',
    # Aliases
    'runCASSIA_hypothesis',
    'run_hypothesis_analysis',
    'summarize_hypothesis_runs',
]
