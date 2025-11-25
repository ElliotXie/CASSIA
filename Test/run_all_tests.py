"""
CASSIA Test Runner - Run All Tests
==================================
Run all CASSIA tests sequentially and generate a summary report.

Usage:
    python run_all_tests.py
    python run_all_tests.py --skip 03,04
    python run_all_tests.py --only 01,02
"""

import sys
import os
import subprocess
import time
import json
import argparse
from pathlib import Path
from datetime import datetime


def get_test_folders() -> list:
    """Get all test folders sorted by number."""
    test_root = Path(__file__).parent

    test_folders = sorted([
        d for d in test_root.iterdir()
        if d.is_dir() and d.name[0].isdigit()
    ])

    return test_folders


def run_python_test(test_folder: Path) -> dict:
    """Run the Python test and return results."""
    test_scripts = list(test_folder.glob("test_*.py"))
    if not test_scripts:
        return {
            "name": test_folder.name,
            "status": "skipped",
            "reason": "No test script found",
            "duration": 0
        }

    test_script = test_scripts[0]
    start_time = time.time()

    try:
        result = subprocess.run(
            [sys.executable, str(test_script)],
            cwd=str(test_folder),
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout per test
        )

        duration = time.time() - start_time

        return {
            "name": test_folder.name,
            "status": "passed" if result.returncode == 0 else "failed",
            "returncode": result.returncode,
            "duration": round(duration, 2),
            "stdout": result.stdout[-2000:] if result.stdout else "",  # Last 2000 chars
            "stderr": result.stderr[-1000:] if result.stderr else ""
        }

    except subprocess.TimeoutExpired:
        return {
            "name": test_folder.name,
            "status": "timeout",
            "duration": 600,
            "reason": "Test exceeded 10 minute timeout"
        }
    except Exception as e:
        return {
            "name": test_folder.name,
            "status": "error",
            "duration": time.time() - start_time,
            "reason": str(e)
        }


def print_summary(results: list, total_duration: float):
    """Print a summary of test results."""
    print("\n" + "=" * 60)
    print("CASSIA TEST SUITE - SUMMARY")
    print("=" * 60)

    passed = sum(1 for r in results if r["status"] == "passed")
    failed = sum(1 for r in results if r["status"] == "failed")
    errors = sum(1 for r in results if r["status"] in ["error", "timeout"])
    skipped = sum(1 for r in results if r["status"] == "skipped")

    print(f"\nResults: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped")
    print(f"Total Duration: {total_duration:.1f}s")

    print("\n" + "-" * 60)
    print(f"{'Test':<40} {'Status':<10} {'Duration':<10}")
    print("-" * 60)

    for result in results:
        status_symbol = {
            "passed": "[OK]",
            "failed": "[X]",
            "error": "[ERR]",
            "timeout": "[TO]",
            "skipped": "[--]"
        }.get(result["status"], "[?]")

        duration_str = f"{result['duration']:.1f}s" if result["duration"] else "N/A"
        print(f"{result['name']:<40} {status_symbol:<10} {duration_str:<10}")

    print("-" * 60)

    if passed == len(results):
        print("\nAll tests PASSED!")
    elif failed + errors > 0:
        print(f"\n{failed + errors} test(s) FAILED or had ERRORS")

    return passed == len([r for r in results if r["status"] != "skipped"])


def save_report(results: list, total_duration: float):
    """Save a JSON report of the test results."""
    test_root = Path(__file__).parent
    report_path = test_root / "test_report.json"

    report = {
        "timestamp": datetime.now().isoformat(),
        "total_duration": round(total_duration, 2),
        "summary": {
            "passed": sum(1 for r in results if r["status"] == "passed"),
            "failed": sum(1 for r in results if r["status"] == "failed"),
            "errors": sum(1 for r in results if r["status"] in ["error", "timeout"]),
            "skipped": sum(1 for r in results if r["status"] == "skipped"),
        },
        "results": results
    }

    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"\nReport saved to: {report_path}")


def parse_test_list(test_string: str) -> list:
    """Parse a comma-separated list of test numbers."""
    if not test_string:
        return []
    return [t.strip() for t in test_string.split(',')]


def main():
    parser = argparse.ArgumentParser(description="Run CASSIA test suite")
    parser.add_argument("--skip", type=str, help="Comma-separated list of tests to skip (e.g., 03,04)")
    parser.add_argument("--only", type=str, help="Comma-separated list of tests to run (e.g., 01,02)")
    parser.add_argument("--no-report", action="store_true", help="Don't save JSON report")
    args = parser.parse_args()

    print("=" * 60)
    print("CASSIA TEST SUITE")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Get test folders
    test_folders = get_test_folders()

    # Apply filters
    skip_tests = parse_test_list(args.skip)
    only_tests = parse_test_list(args.only)

    if only_tests:
        test_folders = [
            f for f in test_folders
            if any(f.name.startswith(t) for t in only_tests)
        ]
    elif skip_tests:
        test_folders = [
            f for f in test_folders
            if not any(f.name.startswith(t) for t in skip_tests)
        ]

    print(f"\nRunning {len(test_folders)} tests:")
    for folder in test_folders:
        print(f"  - {folder.name}")

    # Run tests
    results = []
    total_start = time.time()

    for i, test_folder in enumerate(test_folders, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(test_folders)}] Running: {test_folder.name}")
        print("=" * 60)

        result = run_python_test(test_folder)
        results.append(result)

        # Print immediate status
        if result["status"] == "passed":
            print(f"\n[OK] {test_folder.name} PASSED ({result['duration']:.1f}s)")
        else:
            print(f"\n[X] {test_folder.name} {result['status'].upper()}")
            if "reason" in result:
                print(f"    Reason: {result['reason']}")

    total_duration = time.time() - total_start

    # Print summary
    all_passed = print_summary(results, total_duration)

    # Save report
    if not args.no_report:
        save_report(results, total_duration)

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
