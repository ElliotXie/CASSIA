# CASSIA Reports Module
# Report generation

from .generate_batch_report import generate_batch_html_report_from_data
from .generate_reports import (
    generate_analysis_html_report,
    process_single_report,
    generate_score_index_page,
    runCASSIA_generate_score_report
)
from .generate_hypothesis_report import create_html_report
from .generate_report_uncertainty import generate_uq_html_report

# Alias for backward compatibility
generate_hypothesis_report = create_html_report

__all__ = [
    'generate_batch_html_report_from_data',
    'generate_analysis_html_report',
    'process_single_report',
    'generate_score_index_page',
    'runCASSIA_generate_score_report',
    'create_html_report',
    'generate_hypothesis_report',
    'generate_uq_html_report',
]
