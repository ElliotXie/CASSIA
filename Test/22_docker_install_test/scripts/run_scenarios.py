#!/usr/bin/env python3
"""
CASSIA Docker Installation Test Orchestrator
=============================================

This script orchestrates Docker-based installation tests for the CASSIA R package.
It builds Docker images for various scenarios and runs the installation test in each.

Usage:
    python run_scenarios.py                     # Run all Linux scenarios
    python run_scenarios.py --scenario conda    # Run specific scenario
    python run_scenarios.py --rebuild           # Force rebuild images
    python run_scenarios.py --api-key KEY       # Include annotation test
    python run_scenarios.py --windows           # Include Windows scenario

Scenarios:
    1. minimal        - R only, no Python (tests auto-setup)
    2. system_python  - R + system Python 3.10
    3. virtualenv     - R + Python 3.10 + virtualenv
    4. conda          - Miniconda + R via conda-forge
    5. old_python     - R + Python 3.8 (should fail gracefully)
    6. windows        - Windows Server Core (optional, requires Windows containers)
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Configuration
SCRIPT_DIR = Path(__file__).parent
TEST_DIR = SCRIPT_DIR.parent
DOCKERFILES_DIR = TEST_DIR / "Dockerfiles"
RESULTS_DIR = TEST_DIR / "results"

# Scenario definitions
LINUX_SCENARIOS = {
    "minimal": {
        "dockerfile": DOCKERFILES_DIR / "linux" / "Dockerfile.minimal",
        "description": "R only, no Python (tests auto-setup capability)",
        "expected_success": False,  # May fail - tests error messaging
    },
    "system_python": {
        "dockerfile": DOCKERFILES_DIR / "linux" / "Dockerfile.system_python",
        "description": "R + system Python 3.10",
        "expected_success": True,
    },
    "virtualenv": {
        "dockerfile": DOCKERFILES_DIR / "linux" / "Dockerfile.virtualenv",
        "description": "R + Python 3.10 + virtualenv",
        "expected_success": True,
    },
    "conda": {
        "dockerfile": DOCKERFILES_DIR / "linux" / "Dockerfile.conda",
        "description": "Miniconda with R from conda-forge",
        "expected_success": True,
    },
    "old_python": {
        "dockerfile": DOCKERFILES_DIR / "linux" / "Dockerfile.old_python",
        "description": "R + Python 3.8 (incompatible - should fail)",
        "expected_success": False,
    },
}

WINDOWS_SCENARIOS = {
    "windows": {
        "dockerfile": DOCKERFILES_DIR / "windows" / "Dockerfile.windows",
        "description": "Windows Server Core with R + Python",
        "expected_success": True,
    },
}


def log(msg: str, level: str = "INFO"):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] [{level}] {msg}")


def run_command(cmd: List[str], timeout: int = 600, capture: bool = True) -> Tuple[int, str, str]:
    """Run a command and return exit code, stdout, stderr."""
    log(f"Running: {' '.join(cmd[:5])}...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=capture,
            text=True,
            timeout=timeout,
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -1, "", f"Command timed out after {timeout}s"
    except Exception as e:
        return -1, "", str(e)


def check_docker() -> bool:
    """Check if Docker is available and running."""
    code, stdout, stderr = run_command(["docker", "info"], timeout=30)
    if code != 0:
        log("Docker is not running or not installed", "ERROR")
        log(f"Error: {stderr}", "ERROR")
        return False
    log("Docker is available")
    return True


def build_image(scenario: str, dockerfile: Path, rebuild: bool = False) -> Tuple[bool, str]:
    """Build Docker image for a scenario."""
    tag = f"cassia-test-{scenario}"

    # Check if image already exists
    if not rebuild:
        code, stdout, _ = run_command(["docker", "images", "-q", tag], timeout=30)
        if code == 0 and stdout.strip():
            log(f"Image {tag} already exists (use --rebuild to force)")
            return True, tag

    log(f"Building image: {tag}")

    # Build context is the test directory (contains scripts/)
    build_context = TEST_DIR

    cmd = [
        "docker", "build",
        "-t", tag,
        "-f", str(dockerfile),
        str(build_context)
    ]

    code, stdout, stderr = run_command(cmd, timeout=1200)  # 20 min timeout for builds

    if code != 0:
        log(f"Failed to build {tag}", "ERROR")
        log(f"Error: {stderr[:500]}", "ERROR")
        return False, stderr

    log(f"Successfully built {tag}")
    return True, tag


def run_container(
    tag: str,
    scenario: str,
    api_key: Optional[str] = None,
    provider: str = "openrouter",
    model: str = "google/gemini-2.0-flash-001",
    timeout: int = 600
) -> Dict:
    """Run a container and capture results."""
    log(f"Running container: {tag}")

    # Build docker run command
    cmd = [
        "docker", "run",
        "--rm",
        "-e", f"CASSIA_SCENARIO={scenario}",
    ]

    # Add API key if provided
    if api_key:
        cmd.extend([
            "-e", f"CASSIA_API_KEY={api_key}",
            "-e", f"CASSIA_PROVIDER={provider}",
            "-e", f"CASSIA_MODEL={model}",
        ])

    cmd.append(tag)

    start_time = time.time()

    try:
        code, stdout, stderr = run_command(cmd, timeout=timeout)
    except Exception as e:
        log(f"Error running container: {e}", "ERROR")
        code, stdout, stderr = -1, "", str(e)

    duration = time.time() - start_time

    # Ensure stdout/stderr are strings
    stdout = stdout or ""
    stderr = stderr or ""

    # Parse JSON results from output
    results = {
        "scenario": scenario,
        "exit_code": code,
        "duration": duration,
        "raw_stdout": stdout,
        "raw_stderr": stderr,
        "parsed_results": None,
    }

    # Try to extract JSON results
    try:
        json_match = re.search(r"---JSON_RESULTS_START---\s*(.*?)\s*---JSON_RESULTS_END---", stdout, re.DOTALL)
        if json_match:
            results["parsed_results"] = json.loads(json_match.group(1))
    except (json.JSONDecodeError, TypeError) as e:
        log(f"Failed to parse JSON results: {e}", "WARN")

    return results


def generate_report(all_results: Dict[str, Dict], output_dir: Path) -> Path:
    """Generate summary report."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_dir = output_dir / timestamp
    report_dir.mkdir(parents=True, exist_ok=True)

    # Summary
    summary = {
        "timestamp": timestamp,
        "total_scenarios": len(all_results),
        "passed": 0,
        "failed": 0,
        "scenarios": {},
    }

    for scenario, result in all_results.items():
        parsed = result.get("parsed_results", {})
        success = parsed.get("overall_success", False) if parsed else result["exit_code"] == 0

        summary["scenarios"][scenario] = {
            "success": success,
            "exit_code": result["exit_code"],
            "duration": result["duration"],
            "stages": parsed.get("stages", {}) if parsed else {},
        }

        if success:
            summary["passed"] += 1
        else:
            summary["failed"] += 1

    # Save summary JSON
    summary_file = report_dir / "summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    # Save full results
    full_results_file = report_dir / "full_results.json"
    # Remove raw output from saved results (too verbose)
    clean_results = {}
    for scenario, result in all_results.items():
        clean_results[scenario] = {
            k: v for k, v in result.items()
            if k not in ("raw_stdout", "raw_stderr")
        }
    with open(full_results_file, "w") as f:
        json.dump(clean_results, f, indent=2)

    # Save individual logs
    logs_dir = report_dir / "logs"
    logs_dir.mkdir(exist_ok=True)
    for scenario, result in all_results.items():
        log_file = logs_dir / f"{scenario}.log"
        with open(log_file, "w") as f:
            f.write(f"=== {scenario} ===\n")
            f.write(f"Exit code: {result['exit_code']}\n")
            f.write(f"Duration: {result['duration']:.2f}s\n\n")
            f.write("=== STDOUT ===\n")
            f.write(result.get("raw_stdout", ""))
            f.write("\n\n=== STDERR ===\n")
            f.write(result.get("raw_stderr", ""))

    # Print summary to console
    print("\n" + "=" * 60)
    print("CASSIA Docker Installation Test Summary")
    print("=" * 60)
    print(f"Total scenarios: {summary['total_scenarios']}")
    print(f"Passed: {summary['passed']}")
    print(f"Failed: {summary['failed']}")
    print("-" * 60)

    for scenario, info in summary["scenarios"].items():
        status = "PASS" if info["success"] else "FAIL"
        print(f"  {scenario:20} [{status}] ({info['duration']:.1f}s)")

    print("-" * 60)
    print(f"Results saved to: {report_dir}")
    print("=" * 60)

    return report_dir


def main():
    parser = argparse.ArgumentParser(
        description="Run CASSIA Docker installation tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--scenario", "-s",
        choices=list(LINUX_SCENARIOS.keys()) + list(WINDOWS_SCENARIOS.keys()),
        help="Run specific scenario only"
    )
    parser.add_argument(
        "--rebuild", "-r",
        action="store_true",
        help="Force rebuild Docker images"
    )
    parser.add_argument(
        "--api-key", "-k",
        help="API key for annotation test (or set CASSIA_API_KEY env var)"
    )
    parser.add_argument(
        "--provider", "-p",
        default="openrouter",
        help="LLM provider (default: openrouter)"
    )
    parser.add_argument(
        "--model", "-m",
        default="google/gemini-2.0-flash-001",
        help="Model to use (default: google/gemini-2.0-flash-001)"
    )
    parser.add_argument(
        "--windows", "-w",
        action="store_true",
        help="Include Windows scenario (requires Windows container mode)"
    )
    parser.add_argument(
        "--timeout", "-t",
        type=int,
        default=600,
        help="Container timeout in seconds (default: 600)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=RESULTS_DIR,
        help="Output directory for results"
    )

    args = parser.parse_args()

    # Get API key from args or environment
    api_key = args.api_key or os.environ.get("CASSIA_API_KEY", "")

    # Check Docker
    if not check_docker():
        sys.exit(1)

    # Determine scenarios to run
    scenarios = {}
    if args.scenario:
        if args.scenario in LINUX_SCENARIOS:
            scenarios[args.scenario] = LINUX_SCENARIOS[args.scenario]
        elif args.scenario in WINDOWS_SCENARIOS:
            scenarios[args.scenario] = WINDOWS_SCENARIOS[args.scenario]
    else:
        scenarios = LINUX_SCENARIOS.copy()
        if args.windows:
            scenarios.update(WINDOWS_SCENARIOS)

    log(f"Running {len(scenarios)} scenario(s)")
    if api_key:
        log("API key provided - will run annotation test")
    else:
        log("No API key - skipping annotation test")

    # Build and run each scenario
    all_results = {}

    for scenario, config in scenarios.items():
        print("\n" + "=" * 60)
        log(f"Scenario: {scenario}")
        log(f"Description: {config['description']}")
        print("=" * 60)

        # Build image
        success, tag_or_error = build_image(
            scenario,
            config["dockerfile"],
            rebuild=args.rebuild
        )

        if not success:
            all_results[scenario] = {
                "scenario": scenario,
                "exit_code": -1,
                "duration": 0,
                "raw_stdout": "",
                "raw_stderr": f"Build failed: {tag_or_error}",
                "parsed_results": None,
            }
            continue

        # Run container
        result = run_container(
            tag_or_error,
            scenario,
            api_key=api_key,
            provider=args.provider,
            model=args.model,
            timeout=args.timeout,
        )

        # Ensure result is not None
        if result is None:
            result = {
                "scenario": scenario,
                "exit_code": -1,
                "duration": 0,
                "raw_stdout": "",
                "raw_stderr": "run_container returned None",
                "parsed_results": None,
            }

        all_results[scenario] = result

        # Log result
        parsed = result.get("parsed_results") or {}
        if parsed.get("overall_success"):
            log(f"Scenario {scenario}: PASSED", "INFO")
        else:
            expected = config["expected_success"]
            if not expected:
                log(f"Scenario {scenario}: FAILED (expected)", "INFO")
            else:
                log(f"Scenario {scenario}: FAILED (unexpected)", "WARN")

    # Generate report
    report_dir = generate_report(all_results, args.output)

    # Return appropriate exit code
    failed_unexpected = 0
    for s, r in all_results.items():
        parsed = r.get("parsed_results") or {}
        success = parsed.get("overall_success", r.get("exit_code", -1) == 0)
        if not success and scenarios[s]["expected_success"]:
            failed_unexpected += 1

    sys.exit(0 if failed_unexpected == 0 else 1)


if __name__ == "__main__":
    main()
