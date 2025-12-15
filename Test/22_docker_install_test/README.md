# Test 22: Docker Installation Test

## Purpose
Simulates fresh user installations of the CASSIA R package in various Docker environments to catch installation issues before users encounter them.

## What it Tests
- GitHub installation: `devtools::install_github("ElliotXie/CASSIA")`
- Library loading: `library(CASSIA)` triggers Python setup
- Python environment auto-configuration
- Single-cluster annotation (with API key)

## Prerequisites

### Required
- **Docker Desktop** installed and running
- ~5 GB disk space for Docker images
- Internet connection (downloads from GitHub, CRAN, PyPI)

### Optional
- **API Key** (OpenRouter/Anthropic/OpenAI) for annotation test

## Test Scenarios

| Scenario | Environment | Expected Result |
|----------|-------------|-----------------|
| **minimal** | R only, no Python | Tests error messaging |
| **system_python** | R + Python 3.10 | PASS |
| **virtualenv** | R + Python 3.10 + virtualenv | PASS |
| **conda** | Miniconda + R via conda-forge | PASS (best compatibility) |
| **old_python** | R + Python 3.8 | FAIL with "Python 3.9+ required" |
| **windows** | Windows Server Core (optional) | PASS |

## Running the Tests

### Quick Start (Python)
```bash
cd Test/22_docker_install_test
python test_docker_install.py
```

### Quick Start (R)
```bash
cd Test/22_docker_install_test
Rscript test_docker_install.R
```

### Run Specific Scenario
```bash
python scripts/run_scenarios.py --scenario conda
```

### Include Annotation Test
```bash
# Set API key
export CASSIA_API_KEY="your-openrouter-api-key"

# Run with annotation
python scripts/run_scenarios.py --api-key $CASSIA_API_KEY
```

### Run All Scenarios (Bash)
```bash
./scripts/run_scenarios.sh
```

## Command Line Options

### Python Orchestrator (`run_scenarios.py`)
```
--scenario, -s    Run specific scenario (minimal, system_python, etc.)
--rebuild, -r     Force rebuild Docker images
--api-key, -k     API key for annotation test
--provider, -p    LLM provider (default: openrouter)
--model, -m       Model to use (default: google/gemini-2.0-flash-001)
--windows, -w     Include Windows scenario
--timeout, -t     Container timeout in seconds (default: 600)
--output, -o      Output directory for results
```

### Environment Variables
```bash
CASSIA_API_KEY      # API key for annotation test
CASSIA_PROVIDER     # LLM provider (default: openrouter)
CASSIA_MODEL        # Model to use
```

## Test Stages

Each scenario runs through 5 stages:

1. **Stage 1: Install devtools** - Install R dependencies from CRAN
2. **Stage 2: Install CASSIA** - Install from GitHub
3. **Stage 3: Load library** - `library(CASSIA)` triggers .onLoad()
4. **Stage 4: Verify Python** - Check Python environment setup
5. **Stage 5: Annotation** - Run single-cluster test (optional, needs API key)

## Results

Results are saved to `results/<timestamp>/`:
```
results/20241214_120000/
├── summary.json          # Pass/fail summary
├── full_results.json     # Detailed results
└── logs/
    ├── minimal.log
    ├── system_python.log
    ├── virtualenv.log
    ├── conda.log
    └── old_python.log
```

## Timing

| Phase | First Run | Subsequent |
|-------|-----------|------------|
| Image Build | 5-15 min/image | Cached |
| Container Run | 3-10 min | 3-10 min |
| **Total (5 scenarios)** | 30-60 min | 15-30 min |

## Troubleshooting

### Docker not running
```
ERROR: Docker is not available or not running
```
**Solution:** Start Docker Desktop

### Windows containers
The Windows scenario requires Docker Desktop to be in "Windows container mode":
1. Right-click Docker Desktop tray icon
2. Select "Switch to Windows containers..."

### Build failures
Check the build log in `results/<timestamp>/logs/<scenario>_build.log`

Common issues:
- Network timeout → Retry or check internet connection
- Disk space → Free up space or prune Docker images

### Clean up Docker images
```bash
# Remove all CASSIA test images
docker rmi $(docker images -q "cassia-test-*")

# Prune all unused images
docker image prune -a
```

## File Structure
```
22_docker_install_test/
├── README.md                           # This file
├── Dockerfiles/
│   ├── linux/
│   │   ├── Dockerfile.minimal          # R only
│   │   ├── Dockerfile.system_python    # R + Python 3.10
│   │   ├── Dockerfile.virtualenv       # R + Python + virtualenv
│   │   ├── Dockerfile.conda            # Miniconda + R
│   │   └── Dockerfile.old_python       # R + Python 3.8
│   └── windows/
│       └── Dockerfile.windows          # Windows Server Core
├── scripts/
│   ├── install_test.R                  # R test script (runs in container)
│   ├── run_scenarios.py                # Python orchestrator
│   └── run_scenarios.sh                # Bash runner
├── test_docker_install.py              # Python test wrapper
├── test_docker_install.R               # R test wrapper
└── results/                            # Test results (gitignored)
```

## Cross-Platform Notes

| Platform | Docker Support | Notes |
|----------|---------------|-------|
| **Linux** | Full | Native Docker, best performance |
| **Windows** | Full | Docker Desktop, Linux + Windows containers |
| **macOS** | Linux only | Docker Desktop, cannot run macOS containers |

For true macOS testing, use GitHub Actions with `macos-latest` runner.
