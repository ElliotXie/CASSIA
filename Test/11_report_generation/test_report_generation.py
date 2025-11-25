"""
CASSIA Test 11: Report Generation
=================================
Tests the generate_batch_html_report function for generating interactive HTML
reports from batch analysis results WITHOUT re-running the analysis.

Usage:
    python test_report_generation.py

Functions tested:
- generate_batch_html_report(): Generate HTML report from CSV file
- generate_batch_html_report_from_data(): Generate HTML report from data
"""

import sys
import time
import os
import csv
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_cluster_markers, get_all_cluster_names
from test_utils import (
    setup_cassia_imports,
    load_config,
    print_test_header,
    print_test_result,
    print_config_summary
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata
)

# Setup CASSIA imports
setup_cassia_imports()


def find_existing_batch_results():
    """Find existing batch results from Test 02."""
    test_02_dir = Path(__file__).parent.parent / "02_batch_annotation" / "results"

    if not test_02_dir.exists():
        return None

    # Find most recent results folder
    result_folders = sorted(test_02_dir.glob("*/"), reverse=True)

    for folder in result_folders:
        # Look for the full CSV file
        full_csv = folder / "batch_results_full.csv"
        if full_csv.exists():
            return str(full_csv)

    return None


def create_sample_batch_data():
    """Create sample batch data for testing if no existing results."""
    cluster_names = get_all_cluster_names()
    sample_data = []

    for cluster_name in cluster_names[:3]:  # Use first 3 clusters
        markers = get_cluster_markers(cluster_name)
        marker_list = ", ".join(markers[:15])

        sample_data.append({
            'True Cell Type': cluster_name,
            'Predicted Main Cell Type': cluster_name.title(),
            'Predicted Sub Cell Types': f"{cluster_name} subtype 1, {cluster_name} subtype 2",
            'Possible Mixed Cell Types': '',
            'Marker Number': str(len(markers[:15])),
            'Marker List': marker_list,
            'Iterations': '1',
            'Model': 'google/gemini-2.5-flash',
            'Provider': 'openrouter',
            'Tissue': 'large intestine',
            'Species': 'human',
            'Additional Info': '',
            'Conversation History': f'Final Annotation Agent: Analysis of {cluster_name} markers shows characteristic expression pattern. | Coupling Validator: VALIDATION PASSED - Cell type identification is consistent. | Formatting Agent: {{"main_cell_type": "{cluster_name.title()}", "sub_cell_types": ["{cluster_name} subtype 1", "{cluster_name} subtype 2"], "possible_mixed_cell_types": []}}'
        })

    return sample_data


def run_report_generation_test():
    """Test HTML Report Generation functionality."""
    print_test_header("11 - Report Generation")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Import CASSIA function
    from generate_batch_report import generate_batch_html_report, generate_batch_html_report_from_data

    # Create results directory
    results_dir = create_results_dir("11_report_generation")
    print(f"Results will be saved to: {results_dir}")

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    report_results = {}

    try:
        # Test 1: Try to find existing batch results
        print(f"\n--- Test 1: Find existing batch results ---")
        existing_csv = find_existing_batch_results()

        if existing_csv:
            print(f"  Found existing results: {Path(existing_csv).name}")

            # Count rows in CSV
            with open(existing_csv, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            print(f"  CSV contains: {len(rows)} clusters")

            # Generate report from existing CSV
            print(f"\n--- Test 2: Generate HTML report from CSV ---")
            output_html = str(results_dir / "test_report_from_csv.html")

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
                    'input_file': str(existing_csv),
                    'output_file': result_path,
                    'file_size_bytes': file_size,
                    'clusters': len(rows)
                }
            else:
                errors.append("CSV report not created")

        else:
            print(f"  No existing batch results found")
            print(f"  Will create sample data for testing")

        # Test 3: Generate report from sample data
        print(f"\n--- Test 3: Generate HTML report from data ---")
        sample_data = create_sample_batch_data()
        print(f"  Created sample data: {len(sample_data)} clusters")

        output_html_data = str(results_dir / "test_report_from_data.html")

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
                'clusters': len(sample_data)
            }
        else:
            errors.append("Data report not created")

        # Test 4: Validate HTML content
        print(f"\n--- Test 4: Validate HTML content ---")
        report_to_check = result_path_data if Path(result_path_data).exists() else (
            report_results.get('csv_report', {}).get('output_file')
        )

        if report_to_check and Path(report_to_check).exists():
            with open(report_to_check, 'r', encoding='utf-8') as f:
                html_content = f.read()

            # Check for key HTML elements
            checks = {
                'DOCTYPE': '<!DOCTYPE html>' in html_content,
                'report_header': 'report-header' in html_content,
                'cluster_cards': 'cluster-card' in html_content,
                'search_functionality': 'search-input' in html_content,
                'filter_dropdowns': 'filter-select' in html_content,
                'modal_popups': 'modal-overlay' in html_content,
                'javascript': '<script>' in html_content,
                'css_styles': '<style>' in html_content
            }

            print(f"  Validation results:")
            all_passed = True
            for check_name, passed in checks.items():
                status_str = "PASS" if passed else "FAIL"
                print(f"    {check_name}: {status_str}")
                if not passed:
                    all_passed = False

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

        # Summary
        print(f"\n--- Report Generation Summary ---")
        if 'csv_report' in report_results:
            print(f"  CSV Report: {Path(report_results['csv_report']['output_file']).name}")
        if 'data_report' in report_results:
            print(f"  Data Report: {Path(report_results['data_report']['output_file']).name}")
        print(f"  Validation: {'PASSED' if report_results.get('validation', {}).get('all_passed') else 'FAILED'}")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="report_generation",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[],
        errors=errors
    )
    save_test_metadata(results_dir, metadata)

    save_test_results(results_dir, {
        "report_results": report_results
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_report_generation_test()
    sys.exit(0 if success else 1)
