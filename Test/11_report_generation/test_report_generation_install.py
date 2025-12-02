"""
CASSIA Test 11: Report Generation (PIP INSTALL MODE)
=====================================================
Tests the generate_batch_html_report function using pip-installed CASSIA.

Usage:
    python test_report_generation_install.py
"""

import sys
import time
import csv
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_cluster_markers, get_all_cluster_names
from test_utils import (
    load_config,
    print_test_header,
    print_test_result,
    print_config_summary,
    verify_cassia_pip_install
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA


def find_existing_batch_results():
    """Find existing batch results from Test 02."""
    test_02_dir = Path(__file__).parent.parent / "02_batch_annotation" / "results"
    if not test_02_dir.exists():
        return None
    result_folders = sorted(test_02_dir.glob("*/"), reverse=True)
    for folder in result_folders:
        full_csv = folder / "batch_results_full.csv"
        if full_csv.exists():
            return str(full_csv)
    return None


def create_sample_batch_data():
    """Create sample batch data for testing."""
    cluster_names = get_all_cluster_names()
    sample_data = []
    for cluster_name in cluster_names[:3]:
        markers = get_cluster_markers(cluster_name)
        marker_list = ", ".join(markers[:15])
        sample_data.append({
            'Cluster ID': cluster_name,
            'Predicted General Cell Type': cluster_name.title(),
            'Predicted Detailed Cell Type': f"{cluster_name} subtype 1, {cluster_name} subtype 2",
            'Possible Mixed Cell Types': '',
            'Marker Number': str(len(markers[:15])),
            'Marker List': marker_list,
            'Iterations': '1',
            'Model': 'google/gemini-2.5-flash',
            'Provider': 'openrouter',
            'Tissue': 'large intestine',
            'Species': 'human',
            'Additional Info': '',
            'Conversation History': f'Analysis of {cluster_name} markers.'
        })
    return sample_data


def run_report_generation_test(results_dir):
    """Test HTML Report Generation using pip-installed CASSIA."""
    print_test_header("11 - Report Generation (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Import report generation functions
    from CASSIA.reports.generate_batch_report import generate_batch_html_report, generate_batch_html_report_from_data

    # Results directory passed in from main()
    print(f"Results will be saved to: {results_dir['base']}")

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    report_results = {}

    try:
        # Test 1: Find existing batch results
        print(f"\n--- Test 1: Find existing batch results ---")
        existing_csv = find_existing_batch_results()

        if existing_csv:
            print(f"  Found existing results: {Path(existing_csv).name}")

            # Generate report from existing CSV
            print(f"\n--- Test 2: Generate HTML report from CSV ---")
            output_html = str(results_dir['outputs'] / "test_report_from_csv.html")

            result_path = generate_batch_html_report(
                full_csv_path=existing_csv,
                output_path=output_html,
                report_title="CASSIA Test Report (from CSV)"
            )

            if Path(result_path).exists():
                file_size = Path(result_path).stat().st_size
                print(f"  Generated: {Path(result_path).name}")
                print(f"  File size: {file_size / 1024:.1f} KB")
                report_results['csv_report'] = {
                    'source': 'existing_csv',
                    'output_file': result_path,
                    'file_size_bytes': file_size,
                }
        else:
            print(f"  No existing batch results found")

        # Test 3: Generate report from sample data
        print(f"\n--- Test 3: Generate HTML report from data ---")
        sample_data = create_sample_batch_data()
        print(f"  Created sample data: {len(sample_data)} clusters")

        output_html_data = str(results_dir['outputs'] / "test_report_from_data.html")

        result_path_data = generate_batch_html_report_from_data(
            rows=sample_data,
            output_path=output_html_data,
            report_title="CASSIA Test Report (from Data)"
        )

        if Path(result_path_data).exists():
            file_size_data = Path(result_path_data).stat().st_size
            print(f"  Generated: {Path(result_path_data).name}")
            print(f"  File size: {file_size_data / 1024:.1f} KB")
            report_results['data_report'] = {
                'source': 'sample_data',
                'output_file': result_path_data,
                'file_size_bytes': file_size_data,
            }

        # Validate HTML content
        print(f"\n--- Test 4: Validate HTML content ---")
        report_to_check = result_path_data if Path(result_path_data).exists() else (
            report_results.get('csv_report', {}).get('output_file')
        )

        if report_to_check and Path(report_to_check).exists():
            with open(report_to_check, 'r', encoding='utf-8') as f:
                html_content = f.read()

            checks = {
                'DOCTYPE': '<!DOCTYPE html>' in html_content,
                'javascript': '<script>' in html_content,
                'css_styles': '<style>' in html_content
            }

            all_passed = all(checks.values())
            report_results['validation'] = {
                'file_checked': report_to_check,
                'checks': checks,
                'all_passed': all_passed
            }

            if all_passed:
                status = "passed"
            else:
                status = "failed"
                errors.append("Some HTML validation checks failed")
        else:
            status = "failed"
            errors.append("No HTML file to validate")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="report_generation_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[],
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)

    save_test_results(results_dir['base'], {
        "report_results": report_results,
        "mode": "pip_install"
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("11_report_generation")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_report_generation_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
