#!/usr/bin/env python3
"""
CASSIA Test Framework - Master Test Runner
==========================================

Run all or selected CASSIA tests with comprehensive reporting.

Usage:
    python run_all_tests.py                    # Run all implemented tests
    python run_all_tests.py 01 02 03           # Run specific tests
    python run_all_tests.py --quick            # Run quick tests only
    python run_all_tests.py --verbose          # Verbose output

Test Suites:
    01: runCASSIA_batch - Basic batch annotation
    02: runCASSIA_pipeline - Full end-to-end pipeline
    03: annotation_boost - Iterative deep-dive analysis
    04: merge_annotations - Multi-level annotation merging
    05: uncertainty_quantification - Stability testing
    06: subclustering - Hierarchical annotation
    07-11: Additional tests (placeholders)
"""

import os
import sys
import subprocess
import time
from pathlib import Path
from datetime import datetime
import json
import argparse

# ANSI color codes
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
BOLD = '\033[1m'
RESET = '\033[0m'


class TestRunner:
    """Manages execution of CASSIA test suite."""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.test_dir = Path(__file__).parent
        self.results = {}
        self.start_time = None
        self.end_time = None

        # Define available tests
        self.tests = {
            '01': {
                'name': 'runCASSIA_batch',
                'dir': '01_runCASSIA_batch',
                'script': 'test_batch.py',
                'priority': 'high',
                'runtime': '3-5 min',
                'implemented': True
            },
            '02': {
                'name': 'runCASSIA_pipeline',
                'dir': '02_runCASSIA_pipeline',
                'script': 'test_pipeline.py',
                'priority': 'high',
                'runtime': '10-15 min',
                'implemented': True
            },
            '03': {
                'name': 'annotation_boost',
                'dir': '03_annotation_boost',
                'script': 'test_annotation_boost.py',
                'priority': 'high',
                'runtime': '5-10 min',
                'implemented': True
            },
            '04': {
                'name': 'merge_annotations',
                'dir': '04_merge_annotations',
                'script': 'test_merge_annotations.py',
                'priority': 'medium',
                'runtime': '3-5 min',
                'implemented': True
            },
            '05': {
                'name': 'uncertainty_quantification',
                'dir': '05_uncertainty_quantification',
                'script': 'test_uq_batch.py',
                'priority': 'medium',
                'runtime': '15-20 min',
                'implemented': True
            },
            '06': {
                'name': 'subclustering',
                'dir': '06_subclustering',
                'script': 'test_subclustering.py',
                'priority': 'medium',
                'runtime': '3-5 min',
                'implemented': True
            },
            '07': {
                'name': 'celltype_comparison',
                'dir': '07_celltype_comparison',
                'priority': 'low',
                'runtime': 'TBD',
                'implemented': False
            },
            '08': {
                'name': 'symphony_compare',
                'dir': '08_symphony_compare',
                'priority': 'low',
                'runtime': 'TBD',
                'implemented': False
            },
            '09': {
                'name': 'llm_utils',
                'dir': '09_llm_utils',
                'priority': 'low',
                'runtime': 'TBD',
                'implemented': False
            },
            '10': {
                'name': 'model_settings',
                'dir': '10_model_settings',
                'priority': 'low',
                'runtime': 'TBD',
                'implemented': False
            },
            '11': {
                'name': 'report_generation',
                'dir': '11_report_generation',
                'priority': 'low',
                'runtime': 'TBD',
                'implemented': False
            }
        }

    def print_header(self):
        """Print test runner header."""
        print("=" * 80)
        print(f"{BOLD}CASSIA Test Framework - Master Test Runner{RESET}")
        print("=" * 80)
        print(f"Test Directory: {self.test_dir}")
        print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 80)

    def list_tests(self, filter_implemented=False):
        """List all available tests."""
        print("\nAvailable Tests:")
        print("-" * 80)

        for test_id, test_info in self.tests.items():
            if filter_implemented and not test_info['implemented']:
                continue

            status = f"{GREEN}✓ Implemented{RESET}" if test_info['implemented'] else f"{YELLOW}○ Placeholder{RESET}"
            print(f"  {BOLD}{test_id}{RESET}: {test_info['name']:<30} {status}")
            print(f"      Priority: {test_info['priority']:<8} Runtime: {test_info['runtime']}")

        print("-" * 80)

    def run_test(self, test_id):
        """Run a single test."""
        if test_id not in self.tests:
            print(f"{RED}✗ Unknown test: {test_id}{RESET}")
            return False

        test_info = self.tests[test_id]

        if not test_info['implemented']:
            print(f"{YELLOW}○ Skipping {test_id} ({test_info['name']}): Not implemented{RESET}")
            self.results[test_id] = {'status': 'skipped', 'reason': 'not_implemented'}
            return None

        print(f"\n{BLUE}▶ Running Test {test_id}: {test_info['name']}{RESET}")
        print(f"  Expected runtime: {test_info['runtime']}")
        print("-" * 80)

        test_path = self.test_dir / test_info['dir'] / test_info['script']

        if not test_path.exists():
            print(f"{RED}✗ Test script not found: {test_path}{RESET}")
            self.results[test_id] = {'status': 'error', 'reason': 'script_not_found'}
            return False

        # Run test
        start_time = time.time()
        try:
            result = subprocess.run(
                [sys.executable, str(test_path)],
                cwd=test_path.parent,
                capture_output=not self.verbose,
                text=True,
                timeout=1800  # 30 minute timeout
            )

            duration = time.time() - start_time
            success = result.returncode == 0

            self.results[test_id] = {
                'status': 'success' if success else 'failed',
                'duration': duration,
                'returncode': result.returncode
            }

            if success:
                print(f"{GREEN}✓ Test {test_id} passed ({duration:.1f}s){RESET}")
            else:
                print(f"{RED}✗ Test {test_id} failed (code {result.returncode}){RESET}")
                if not self.verbose and result.stderr:
                    print(f"  Error: {result.stderr[:200]}")

            return success

        except subprocess.TimeoutExpired:
            print(f"{RED}✗ Test {test_id} timed out{RESET}")
            self.results[test_id] = {'status': 'timeout'}
            return False

        except Exception as e:
            print(f"{RED}✗ Test {test_id} error: {str(e)}{RESET}")
            self.results[test_id] = {'status': 'error', 'error': str(e)}
            return False

    def run_tests(self, test_ids=None):
        """Run multiple tests."""
        self.start_time = time.time()

        if test_ids is None:
            # Run all implemented tests
            test_ids = [tid for tid, info in self.tests.items() if info['implemented']]

        print(f"\nRunning {len(test_ids)} test(s)...")
        print("=" * 80)

        for test_id in test_ids:
            self.run_test(test_id)

        self.end_time = time.time()

    def print_summary(self):
        """Print test summary."""
        if not self.results:
            print("\nNo tests were run.")
            return

        print("\n" + "=" * 80)
        print(f"{BOLD}TEST SUMMARY{RESET}")
        print("=" * 80)

        passed = sum(1 for r in self.results.values() if r.get('status') == 'success')
        failed = sum(1 for r in self.results.values() if r.get('status') == 'failed')
        skipped = sum(1 for r in self.results.values() if r.get('status') == 'skipped')
        errors = sum(1 for r in self.results.values() if r.get('status') in ['error', 'timeout'])

        total = len(self.results)
        total_time = self.end_time - self.start_time

        print(f"Total Tests: {total}")
        print(f"  {GREEN}✓ Passed: {passed}{RESET}")
        print(f"  {RED}✗ Failed: {failed}{RESET}")
        print(f"  {YELLOW}○ Skipped: {skipped}{RESET}")
        print(f"  {RED}⚠ Errors: {errors}{RESET}")
        print(f"\nTotal Time: {total_time:.1f}s ({total_time/60:.1f}min)")

        # Individual test results
        print("\nDetailed Results:")
        print("-" * 80)
        for test_id, result in self.results.items():
            test_name = self.tests[test_id]['name']
            status = result.get('status', 'unknown')
            duration = result.get('duration', 0)

            if status == 'success':
                print(f"  {GREEN}✓{RESET} {test_id}: {test_name:<30} ({duration:.1f}s)")
            elif status == 'failed':
                print(f"  {RED}✗{RESET} {test_id}: {test_name:<30} (failed)")
            elif status == 'skipped':
                print(f"  {YELLOW}○{RESET} {test_id}: {test_name:<30} (skipped)")
            else:
                print(f"  {RED}⚠{RESET} {test_id}: {test_name:<30} ({status})")

        print("=" * 80)

        # Save results to JSON
        results_file = self.test_dir / "results" / f"test_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        results_file.parent.mkdir(exist_ok=True)

        with open(results_file, 'w') as f:
            json.dump({
                'timestamp': datetime.now().isoformat(),
                'total_time': total_time,
                'summary': {
                    'passed': passed,
                    'failed': failed,
                    'skipped': skipped,
                    'errors': errors
                },
                'results': self.results
            }, f, indent=2)

        print(f"\n📁 Results saved to: {results_file}")

        return passed == total


def main():
    parser = argparse.ArgumentParser(
        description='CASSIA Test Framework - Master Test Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_all_tests.py                 # Run all implemented tests
  python run_all_tests.py 01 02 03        # Run specific tests
  python run_all_tests.py --list          # List available tests
  python run_all_tests.py --quick         # Run quick tests only (01, 04, 06)
  python run_all_tests.py --verbose       # Show detailed output
        """
    )

    parser.add_argument('tests', nargs='*', help='Test IDs to run (e.g., 01 02 03)')
    parser.add_argument('--list', action='store_true', help='List all available tests')
    parser.add_argument('--quick', action='store_true', help='Run quick tests only (<5 min)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('--all', action='store_true', help='Run all implemented tests')

    args = parser.parse_args()

    runner = TestRunner(verbose=args.verbose)
    runner.print_header()

    if args.list:
        runner.list_tests()
        return 0

    # Determine which tests to run
    if args.quick:
        test_ids = ['01', '04', '06']  # Quick tests
    elif args.all or not args.tests:
        test_ids = [tid for tid, info in runner.tests.items() if info['implemented']]
    else:
        test_ids = args.tests

    # Validate test IDs
    invalid_ids = [tid for tid in test_ids if tid not in runner.tests]
    if invalid_ids:
        print(f"{RED}Error: Invalid test IDs: {', '.join(invalid_ids)}{RESET}")
        runner.list_tests(filter_implemented=True)
        return 1

    # Run tests
    runner.list_tests(filter_implemented=True)
    runner.run_tests(test_ids)
    success = runner.print_summary()

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
