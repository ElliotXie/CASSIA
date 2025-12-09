"""
CASSIA Test Runner - Run All Tests
==================================
Run all CASSIA tests sequentially and generate a summary report.

Note: Test 16 (manual package test) is skipped by default.

Usage:
    python D:/CASSIA/Test/run_all_tests.py                 # Run all tests using local development source
    python D:/CASSIA/Test/run_all_tests.py --install       # pip install CASSIA from PyPI, then run tests
    python D:/CASSIA/Test/run_all_tests.py --skip 03,04    # Skip specific tests (in addition to default skips)
    python D:/CASSIA/Test/run_all_tests.py --only 01,02    # Run only specific tests (ignores default skips)
    python D:/CASSIA/Test/run_all_tests.py --only 16       # Run test 16 explicitly
    python D:/CASSIA/Test/run_all_tests.py --no-report     # Don't save JSON/HTML reports
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


def run_python_test(test_folder: Path, mode: str = "development") -> dict:
    """Run the Python test and return results.

    Args:
        test_folder: Path to the test folder
        mode: 'installed' or 'development' - determines which test file to run
    """
    # For install mode, prefer test_*_install.py files
    if mode == "installed":
        install_scripts = list(test_folder.glob("test_*_install.py"))
        if install_scripts:
            test_script = install_scripts[0]
        else:
            # Fall back to regular test file
            test_scripts = [f for f in test_folder.glob("test_*.py") if "_install" not in f.name]
            if not test_scripts:
                return {
                    "name": test_folder.name,
                    "status": "skipped",
                    "reason": "No test script found",
                    "duration": 0
                }
            test_script = test_scripts[0]
    else:
        # Development mode: use regular test files (not _install.py)
        test_scripts = [f for f in test_folder.glob("test_*.py") if "_install" not in f.name]
        if not test_scripts:
            return {
                "name": test_folder.name,
                "status": "skipped",
                "reason": "No test script found",
                "duration": 0
            }
        test_script = test_scripts[0]

    print(f"  Running: {test_script.name}")
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


def save_report(results: list, total_duration: float, mode: str = "development", version_info: dict = None):
    """Save a JSON report of the test results."""
    test_root = Path(__file__).parent

    # Create organized folder structure: reports/python/{mode}/{timestamp}/
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    reports_dir = test_root / "reports" / "python" / mode / timestamp
    reports_dir.mkdir(parents=True, exist_ok=True)

    report_path = reports_dir / "report.json"
    html_path = reports_dir / "report.html"

    report = {
        "timestamp": datetime.now().isoformat(),
        "mode": mode,
        "language": "python",
        "version": version_info.get('version', 'unknown') if version_info else 'unknown',
        "location": version_info.get('location', 'unknown') if version_info else 'unknown',
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

    print(f"\nJSON report saved to: {report_path}")

    # Generate HTML report
    save_html_report(report, html_path)
    print(f"HTML report saved to: {html_path}")

    # Try to open HTML report in browser
    import webbrowser
    try:
        webbrowser.open(f"file://{html_path.resolve()}")
    except:
        pass  # Silently ignore if browser can't be opened


def save_html_report(report: dict, html_path: Path):
    """Generate a nice HTML report of test results."""
    passed = report['summary']['passed']
    failed = report['summary']['failed']
    errors = report['summary']['errors']
    skipped = report['summary']['skipped']
    total = passed + failed + errors + skipped

    # Determine overall status color
    if failed + errors == 0:
        status_color = "#28a745"  # green
        status_text = "ALL PASSED"
    else:
        status_color = "#dc3545"  # red
        status_text = f"{failed + errors} FAILED"

    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CASSIA Test Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: #f5f5f5;
            padding: 20px;
            line-height: 1.6;
        }}
        .container {{ max-width: 1000px; margin: 0 auto; }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px 10px 0 0;
            text-align: center;
        }}
        .header h1 {{ font-size: 2em; margin-bottom: 10px; }}
        .header .meta {{ opacity: 0.9; font-size: 0.9em; }}
        .summary {{
            background: white;
            padding: 25px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
            gap: 15px;
            border-bottom: 1px solid #eee;
        }}
        .stat {{
            text-align: center;
            padding: 15px;
            border-radius: 8px;
        }}
        .stat-value {{ font-size: 2em; font-weight: bold; }}
        .stat-label {{ font-size: 0.85em; color: #666; text-transform: uppercase; }}
        .stat.passed {{ background: #d4edda; color: #155724; }}
        .stat.failed {{ background: #f8d7da; color: #721c24; }}
        .stat.errors {{ background: #fff3cd; color: #856404; }}
        .stat.skipped {{ background: #e2e3e5; color: #383d41; }}
        .stat.total {{ background: #cce5ff; color: #004085; }}
        .status-banner {{
            background: {status_color};
            color: white;
            text-align: center;
            padding: 15px;
            font-size: 1.3em;
            font-weight: bold;
        }}
        .results {{
            background: white;
            border-radius: 0 0 10px 10px;
            overflow: hidden;
        }}
        .result-item {{
            padding: 15px 25px;
            border-bottom: 1px solid #eee;
            display: flex;
            align-items: center;
            gap: 15px;
        }}
        .result-item:last-child {{ border-bottom: none; }}
        .result-item:hover {{ background: #f8f9fa; }}
        .status-badge {{
            padding: 5px 12px;
            border-radius: 20px;
            font-size: 0.8em;
            font-weight: bold;
            min-width: 70px;
            text-align: center;
        }}
        .status-passed {{ background: #28a745; color: white; }}
        .status-failed {{ background: #dc3545; color: white; }}
        .status-error {{ background: #ffc107; color: #212529; }}
        .status-timeout {{ background: #fd7e14; color: white; }}
        .status-skipped {{ background: #6c757d; color: white; }}
        .test-name {{ flex: 1; font-weight: 500; }}
        .test-duration {{ color: #666; font-size: 0.9em; }}
        .error-details {{
            background: #f8f9fa;
            padding: 15px 25px 15px 100px;
            border-bottom: 1px solid #eee;
            font-family: monospace;
            font-size: 0.85em;
            white-space: pre-wrap;
            color: #721c24;
            max-height: 200px;
            overflow-y: auto;
        }}
        .info-section {{
            background: white;
            margin-top: 20px;
            padding: 20px 25px;
            border-radius: 10px;
        }}
        .info-section h3 {{ margin-bottom: 10px; color: #333; }}
        .info-grid {{
            display: grid;
            grid-template-columns: 150px 1fr;
            gap: 8px;
            font-size: 0.9em;
        }}
        .info-label {{ color: #666; }}
        .info-value {{ color: #333; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>CASSIA Test Report</h1>
            <div class="meta">
                {report['timestamp'][:19].replace('T', ' ')} | Mode: {report['mode']} | Version: {report['version']}
            </div>
        </div>
        <div class="summary">
            <div class="stat total">
                <div class="stat-value">{total}</div>
                <div class="stat-label">Total</div>
            </div>
            <div class="stat passed">
                <div class="stat-value">{passed}</div>
                <div class="stat-label">Passed</div>
            </div>
            <div class="stat failed">
                <div class="stat-value">{failed}</div>
                <div class="stat-label">Failed</div>
            </div>
            <div class="stat errors">
                <div class="stat-value">{errors}</div>
                <div class="stat-label">Errors</div>
            </div>
            <div class="stat skipped">
                <div class="stat-value">{skipped}</div>
                <div class="stat-label">Skipped</div>
            </div>
        </div>
        <div class="status-banner">{status_text}</div>
        <div class="results">
'''

    for result in report['results']:
        status = result['status']
        status_class = f"status-{status}"
        duration = f"{result['duration']:.1f}s" if result.get('duration') else "N/A"

        html_content += f'''            <div class="result-item">
                <span class="status-badge {status_class}">{status.upper()}</span>
                <span class="test-name">{result['name']}</span>
                <span class="test-duration">{duration}</span>
            </div>
'''
        # Add error details if failed/error/timeout
        if status in ['failed', 'error', 'timeout']:
            error_text = ""
            if result.get('stderr'):
                error_text = result['stderr']
            elif result.get('reason'):
                error_text = result['reason']

            if error_text:
                # Escape HTML in error message
                error_text = error_text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
                html_content += f'''            <div class="error-details">{error_text}</div>
'''

    html_content += f'''        </div>
        <div class="info-section">
            <h3>Test Environment</h3>
            <div class="info-grid">
                <span class="info-label">Mode:</span>
                <span class="info-value">{report['mode']}</span>
                <span class="info-label">Version:</span>
                <span class="info-value">{report['version']}</span>
                <span class="info-label">Location:</span>
                <span class="info-value">{report['location']}</span>
                <span class="info-label">Total Duration:</span>
                <span class="info-value">{report['total_duration']:.1f}s</span>
            </div>
        </div>
    </div>
</body>
</html>
'''

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)


def parse_test_list(test_string: str) -> list:
    """Parse a comma-separated list of test numbers."""
    if not test_string:
        return []
    return [t.strip() for t in test_string.split(',')]


def get_cassia_version_info(mode: str) -> dict:
    """
    Get CASSIA version information based on test mode.

    Args:
        mode: 'installed' or 'development'

    Returns:
        dict with version, location, and any warnings/errors
    """
    result = {
        'version': 'unknown',
        'location': 'unknown',
        'warning': None,
        'error': None
    }

    if mode == 'installed':
        # Check if CASSIA is already imported (can cause issues)
        if 'CASSIA' in sys.modules:
            result['warning'] = "CASSIA already loaded in current session. Restart Python for clean install test."

        try:
            # Try to get version from pip
            pip_result = subprocess.run(
                [sys.executable, "-m", "pip", "show", "CASSIA"],
                capture_output=True,
                text=True
            )
            if pip_result.returncode == 0:
                for line in pip_result.stdout.split('\n'):
                    if line.startswith('Version:'):
                        result['version'] = line.split(':', 1)[1].strip()
                    elif line.startswith('Location:'):
                        result['location'] = line.split(':', 1)[1].strip()
            else:
                result['error'] = "CASSIA not installed via pip. Run 'pip install CASSIA' first."
        except Exception as e:
            result['error'] = f"Failed to check pip installation: {e}"
    else:
        # Development mode - get version from local source
        try:
            cassia_path = Path(__file__).parent.parent / "CASSIA_python" / "CASSIA"
            init_file = cassia_path / "__init__.py"
            if init_file.exists():
                with open(init_file, 'r') as f:
                    content = f.read()
                    for line in content.split('\n'):
                        if '__version__' in line:
                            # Extract version from __version__ = "x.x.x"
                            import re
                            match = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', line)
                            if match:
                                result['version'] = match.group(1)
                                break
                result['location'] = str(cassia_path)
            else:
                result['error'] = f"Development source not found at {cassia_path}"
        except Exception as e:
            result['error'] = f"Failed to read development version: {e}"

    return result


def reinstall_cassia_from_pypi() -> dict:
    """
    Reinstall CASSIA from PyPI using pip.

    Returns:
        dict with success status and message
    """
    print(f"\nReinstalling CASSIA from PyPI...")
    print("-" * 60)

    try:
        # Uninstall first to ensure clean install
        print("Uninstalling existing CASSIA...")
        uninstall_result = subprocess.run(
            [sys.executable, "-m", "pip", "uninstall", "-y", "CASSIA"],
            capture_output=True,
            text=True
        )
        if uninstall_result.stdout:
            print(uninstall_result.stdout.strip())

        # Install from PyPI
        print("Installing CASSIA from PyPI...")
        install_result = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--upgrade", "CASSIA"],
            capture_output=True,
            text=True
        )

        if install_result.returncode == 0:
            # Show what was installed
            print(install_result.stdout.strip())
            print("[OK] CASSIA installed from PyPI successfully")
            return {
                'success': True,
                'message': "CASSIA installed from PyPI"
            }
        else:
            print(f"[X] Installation failed:")
            print(f"stdout: {install_result.stdout}")
            print(f"stderr: {install_result.stderr}")
            return {
                'success': False,
                'message': f"pip install failed: {install_result.stderr}"
            }

    except Exception as e:
        return {
            'success': False,
            'message': f"Reinstall failed: {e}"
        }


def main():
    parser = argparse.ArgumentParser(description="Run CASSIA test suite")
    parser.add_argument("--skip", type=str, help="Comma-separated list of tests to skip (e.g., 03,04)")
    parser.add_argument("--only", type=str, help="Comma-separated list of tests to run (e.g., 01,02)")
    parser.add_argument("--no-report", action="store_true", help="Don't save JSON report")
    parser.add_argument("--install", action="store_true", help="pip install CASSIA from PyPI and test the installed package")
    args = parser.parse_args()

    # Set test mode environment variable
    test_mode = "installed" if args.install else "development"
    os.environ["CASSIA_TEST_MODE"] = test_mode

    # When --install is used, always reinstall from PyPI first
    if args.install:
        reinstall_result = reinstall_cassia_from_pypi()
        if not reinstall_result['success']:
            print(f"\n[ERROR] {reinstall_result['message']}")
            sys.exit(1)

    # Get CASSIA version info
    version_info = get_cassia_version_info(test_mode)

    print("=" * 60)
    print("CASSIA TEST SUITE (Python)")
    print("=" * 60)
    print(f"Mode: {test_mode}")
    print(f"Version: {version_info['version']}")
    print(f"Location: {version_info['location']}")
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Show warnings or errors
    if version_info['warning']:
        print(f"\n[WARNING] {version_info['warning']}")
    if version_info['error']:
        print(f"\n[ERROR] {version_info['error']}")
        if test_mode == "installed":
            print("Aborting tests. Please install CASSIA first: pip install CASSIA")
            sys.exit(1)

    # Get test folders
    test_folders = get_test_folders()

    # Apply filters
    skip_tests = parse_test_list(args.skip)
    only_tests = parse_test_list(args.only)

    # Default: skip test 16 (manual package test) unless explicitly included
    default_skip = ["16"]
    if not only_tests:  # Only apply default skip if not using --only
        for skip in default_skip:
            if skip not in skip_tests:
                skip_tests.append(skip)

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

        result = run_python_test(test_folder, test_mode)
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
        save_report(results, total_duration, test_mode, version_info)

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
