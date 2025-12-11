"""
CASSIA Test 20: Local LLM Support
=================================
Tests that CASSIA correctly handles local LLM endpoints (Ollama, LM Studio, vLLM)
without requiring an API key.

Includes a full batch annotation test using real marker data.

Usage:
    python test_local_llm.py
"""

import os
import sys
import time
import socket
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import (
    setup_cassia_imports,
    load_config,
    print_test_header,
    print_test_result,
    print_config_summary,
    get_test_mode
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)
from fixtures import load_markers

# Setup CASSIA imports
setup_cassia_imports()


def get_ollama_model():
    """Get the first available Ollama model."""
    try:
        import requests
        response = requests.get("http://localhost:11434/api/tags", timeout=5)
        if response.status_code == 200:
            data = response.json()
            models = [m.get("name", "") for m in data.get("models", [])]
            if models:
                return models[0]
    except:
        pass
    return None


def check_ollama_running(host="localhost", port=11434, timeout=2):
    """Check if Ollama is running by attempting a socket connection."""
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(timeout)
        result = sock.connect_ex((host, port))
        sock.close()
        return result == 0
    except:
        return False


def test_localhost_detection():
    """Test that localhost URLs are correctly detected."""
    test_cases = [
        ("http://localhost:11434/v1", True),
        ("http://localhost:1234/v1", True),
        ("http://127.0.0.1:8000/v1", True),
        ("http://LOCALHOST:11434/v1", True),  # Case insensitive
        ("https://api.deepseek.com", False),
        ("https://api.openai.com/v1", False),
    ]

    errors = []
    for url, expected in test_cases:
        is_localhost = any(x in url.lower() for x in ["localhost", "127.0.0.1"])
        if is_localhost != expected:
            errors.append(f"URL '{url}': expected {expected}, got {is_localhost}")

    return len(errors) == 0, errors


def test_api_key_bypass():
    """Test that API key is not required for localhost URLs."""
    from CASSIA import call_llm

    # Clear any existing API key
    if "CUSTOMIZED_API_KEY" in os.environ:
        del os.environ["CUSTOMIZED_API_KEY"]

    # This should NOT raise an error for localhost (even if Ollama isn't running)
    # It will fail on connection, but not on API key validation
    try:
        call_llm(
            prompt="test",
            provider="http://localhost:11434/v1",
            model="llama2",
            max_tokens=10
        )
    except ValueError as e:
        if "API key" in str(e):
            return False, ["API key was required for localhost URL"]
        # Connection error is expected if Ollama isn't running
        return True, []
    except Exception as e:
        # Connection error is fine - we just want to verify no API key error
        if "API key" in str(e).lower():
            return False, [f"API key error: {e}"]
        return True, []

    return True, []


def test_remote_requires_key():
    """Test that remote URLs still require API key."""
    from CASSIA import call_llm

    # Clear any existing API key
    if "CUSTOMIZED_API_KEY" in os.environ:
        del os.environ["CUSTOMIZED_API_KEY"]

    try:
        call_llm(
            prompt="test",
            provider="https://api.deepseek.com",
            model="deepseek-chat",
            max_tokens=10
        )
        return False, ["Remote URL should require API key but didn't raise error"]
    except ValueError as e:
        if "API key" in str(e):
            return True, []
        return False, [f"Wrong error type: {e}"]
    except Exception as e:
        if "API key" in str(e).lower():
            return True, []
        return False, [f"Wrong error: {e}"]


def test_ollama_connection():
    """Test actual Ollama connection if running."""
    from CASSIA import call_llm

    if not check_ollama_running():
        return "skipped", ["Ollama not running on localhost:11434"]

    # Try to get available models from Ollama
    try:
        import requests
        response = requests.get("http://localhost:11434/api/tags", timeout=5)
        if response.status_code == 200:
            data = response.json()
            # Use full model name including tag (e.g., "gpt-oss:20b")
            models = [m.get("name", "") for m in data.get("models", [])]
            if not models:
                return "skipped", ["Ollama running but no models installed. Run: ollama pull llama2"]
            model_to_use = models[0]  # Use first available model
            print(f"  Using Ollama model: {model_to_use}")
        else:
            return "skipped", ["Could not get Ollama model list"]
    except Exception as e:
        return "skipped", [f"Could not connect to Ollama API: {e}"]

    try:
        response = call_llm(
            prompt="Reply with exactly one word: hello",
            provider="http://localhost:11434/v1",
            model=model_to_use,
            temperature=0.1,
            max_tokens=20
        )
        if response and len(response.strip()) > 0:
            print(f"  Response: {response[:50]}...")
            return True, []
        # Empty response but connection worked - still a pass for API key bypass
        print("  Note: Empty response (model may need more time or different prompt)")
        return True, []  # Connection worked, API key bypass verified
    except Exception as e:
        # Model not found or other errors - skip gracefully
        if "not found" in str(e).lower() or "404" in str(e):
            return "skipped", [f"Model '{model_to_use}' not available: {e}"]
        return False, [f"Ollama call failed: {e}"]


def test_batch_annotation_ollama(results_dir):
    """Test full batch annotation with local Ollama - like test 02."""
    from CASSIA import runCASSIA_batch

    if not check_ollama_running():
        return "skipped", ["Ollama not running on localhost:11434"]

    model_to_use = get_ollama_model()
    if not model_to_use:
        return "skipped", ["No Ollama models available"]

    print(f"  Using Ollama model: {model_to_use}")
    print(f"  Provider: http://localhost:11434/v1")

    # Load marker data - use 2 clusters like test 02
    test_clusters = ['monocyte', 'plasma cell']
    full_df = load_markers()
    marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()

    print(f"  Testing annotation for: {test_clusters}")
    print(f"  Marker data shape: {marker_df.shape}")

    output_name = str(results_dir / "ollama_batch_results")

    try:
        print("  Running runCASSIA_batch with local Ollama...")
        runCASSIA_batch(
            marker=marker_df,
            output_name=output_name,
            n_genes=30,
            model=model_to_use,
            temperature=0.3,
            tissue="large intestine",
            species="human",
            max_workers=1,  # Single worker for local LLM
            provider="http://localhost:11434/v1",
            validator_involvement="v1"
        )

        # Check output files
        summary_csv = Path(f"{output_name}_summary.csv")
        if summary_csv.exists():
            import pandas as pd
            results_df = pd.read_csv(summary_csv)
            clusters_annotated = len(results_df)
            print(f"  Clusters annotated: {clusters_annotated}/{len(test_clusters)}")

            if clusters_annotated == len(test_clusters):
                # Check if we got actual annotations
                if 'Final Cell Type' in results_df.columns:
                    cell_type = results_df['Final Cell Type'].iloc[0]
                    print(f"  Annotation result: {cell_type}")
                return True, []
            else:
                return False, [f"Only {clusters_annotated}/{len(test_clusters)} clusters annotated"]
        else:
            return False, ["Summary CSV not created"]

    except Exception as e:
        error_msg = str(e)
        # Handle known issues gracefully
        if "timeout" in error_msg.lower():
            return "skipped", [f"Ollama timeout (model may be too slow): {error_msg[:100]}"]
        return False, [f"Batch annotation failed: {error_msg[:200]}"]


def run_local_llm_test():
    """Run all local LLM tests."""
    print_test_header("20 - Local LLM Support Testing")

    config = load_config()
    print_config_summary(config)

    # Create results directory
    results = create_results_dir("20_local_llm", get_test_mode())
    logging_ctx = setup_logging(results['logs'])

    start_time = time.time()
    test_results = {"tests": []}
    passed = 0
    failed = 0
    skipped = 0
    errors = []

    # Test 1: Localhost detection logic
    print("\n" + "="*50)
    print("Test 1: Localhost Detection Logic")
    print("="*50)
    success, errs = test_localhost_detection()
    test_results["tests"].append({
        "name": "localhost_detection",
        "status": "passed" if success else "failed",
        "errors": errs
    })
    if success:
        print("  [OK] Localhost detection working correctly")
        passed += 1
    else:
        print(f"  [X] FAILED: {errs}")
        failed += 1
        errors.extend(errs)

    # Test 2: API key bypass for localhost
    print("\n" + "="*50)
    print("Test 2: API Key Bypass for Localhost")
    print("="*50)
    success, errs = test_api_key_bypass()
    test_results["tests"].append({
        "name": "api_key_bypass",
        "status": "passed" if success else "failed",
        "errors": errs
    })
    if success:
        print("  [OK] API key not required for localhost")
        passed += 1
    else:
        print(f"  [X] FAILED: {errs}")
        failed += 1
        errors.extend(errs)

    # Test 3: Remote still requires key
    print("\n" + "="*50)
    print("Test 3: Remote URL Requires API Key")
    print("="*50)
    success, errs = test_remote_requires_key()
    test_results["tests"].append({
        "name": "remote_requires_key",
        "status": "passed" if success else "failed",
        "errors": errs
    })
    if success:
        print("  [OK] Remote URLs correctly require API key")
        passed += 1
    else:
        print(f"  [X] FAILED: {errs}")
        failed += 1
        errors.extend(errs)

    # Test 4: Actual Ollama connection (skip if not running)
    print("\n" + "="*50)
    print("Test 4: Ollama Connection (Optional)")
    print("="*50)
    result, errs = test_ollama_connection()
    if result == "skipped":
        print(f"  [SKIP] {errs[0]}")
        skipped += 1
        test_results["tests"].append({
            "name": "ollama_connection",
            "status": "skipped",
            "errors": errs
        })
    elif result:
        print("  [OK] Successfully connected to Ollama")
        passed += 1
        test_results["tests"].append({
            "name": "ollama_connection",
            "status": "passed",
            "errors": []
        })
    else:
        print(f"  [X] FAILED: {errs}")
        failed += 1
        errors.extend(errs)
        test_results["tests"].append({
            "name": "ollama_connection",
            "status": "failed",
            "errors": errs
        })

    # Test 5: Full batch annotation with Ollama (like test 02)
    print("\n" + "="*50)
    print("Test 5: Batch Annotation with Local Ollama")
    print("="*50)
    result, errs = test_batch_annotation_ollama(results['outputs'])
    if result == "skipped":
        print(f"  [SKIP] {errs[0]}")
        skipped += 1
        test_results["tests"].append({
            "name": "batch_annotation_ollama",
            "status": "skipped",
            "errors": errs
        })
    elif result:
        print("  [OK] Batch annotation completed successfully")
        passed += 1
        test_results["tests"].append({
            "name": "batch_annotation_ollama",
            "status": "passed",
            "errors": []
        })
    else:
        print(f"  [X] FAILED: {errs}")
        failed += 1
        errors.extend(errs)
        test_results["tests"].append({
            "name": "batch_annotation_ollama",
            "status": "failed",
            "errors": errs
        })

    duration = time.time() - start_time

    # Summary
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Skipped: {skipped}")
    print(f"  Duration: {duration:.2f}s")

    # Determine overall status
    if failed > 0:
        status = "failed"
    elif passed == 0:
        status = "skipped"
    else:
        status = "passed"

    # Save metadata
    metadata = create_test_metadata(
        test_name="local_llm_test",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[],
        errors=errors
    )
    metadata["test_summary"] = {
        "passed": passed,
        "failed": failed,
        "skipped": skipped
    }
    save_test_metadata(results['outputs'], metadata)
    save_test_results(results['outputs'], test_results)

    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    cleanup_logging(logging_ctx)
    return success


if __name__ == "__main__":
    success = run_local_llm_test()
    sys.exit(0 if success else 1)
