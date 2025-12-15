#!/bin/bash
# =============================================================================
# CASSIA Docker Installation Test Runner (Bash)
# =============================================================================
# Simple bash script to run Docker installation tests.
# For more options, use the Python orchestrator: run_scenarios.py
#
# Usage:
#   ./run_scenarios.sh                    # Run all Linux scenarios
#   ./run_scenarios.sh conda              # Run specific scenario
#   ./run_scenarios.sh --with-api-key     # Include annotation test
#   CASSIA_API_KEY=xxx ./run_scenarios.sh # With API key from env
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(dirname "$SCRIPT_DIR")"
DOCKERFILES_DIR="$TEST_DIR/Dockerfiles/linux"
RESULTS_DIR="$TEST_DIR/results"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Scenarios
SCENARIOS=("minimal" "system_python" "virtualenv" "conda" "old_python")

# Parse arguments
SPECIFIC_SCENARIO=""
WITH_API_KEY=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --with-api-key)
            WITH_API_KEY=true
            shift
            ;;
        minimal|system_python|virtualenv|conda|old_python|windows)
            SPECIFIC_SCENARIO="$1"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [scenario] [--with-api-key]"
            echo ""
            echo "Scenarios: minimal, system_python, virtualenv, conda, old_python"
            echo ""
            echo "Options:"
            echo "  --with-api-key   Include annotation test (requires CASSIA_API_KEY env var)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Log function
log() {
    echo -e "[$(date +%H:%M:%S)] $1"
}

# Check Docker
if ! docker info &> /dev/null; then
    log "${RED}ERROR: Docker is not running${NC}"
    exit 1
fi
log "${GREEN}Docker is available${NC}"

# Create results directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_RESULTS_DIR="$RESULTS_DIR/$TIMESTAMP"
mkdir -p "$RUN_RESULTS_DIR/logs"

# Determine scenarios to run
if [ -n "$SPECIFIC_SCENARIO" ]; then
    SCENARIOS=("$SPECIFIC_SCENARIO")
fi

log "Running ${#SCENARIOS[@]} scenario(s): ${SCENARIOS[*]}"

# Track results
PASSED=0
FAILED=0
declare -A RESULTS

# Run each scenario
for SCENARIO in "${SCENARIOS[@]}"; do
    log ""
    log "=============================================="
    log "Scenario: ${YELLOW}$SCENARIO${NC}"
    log "=============================================="

    TAG="cassia-test-$SCENARIO"
    DOCKERFILE="$DOCKERFILES_DIR/Dockerfile.$SCENARIO"

    # Check dockerfile exists
    if [ ! -f "$DOCKERFILE" ]; then
        log "${RED}Dockerfile not found: $DOCKERFILE${NC}"
        RESULTS[$SCENARIO]="SKIP"
        continue
    fi

    # Build image
    log "Building image: $TAG"
    if ! docker build -t "$TAG" -f "$DOCKERFILE" "$TEST_DIR" > "$RUN_RESULTS_DIR/logs/${SCENARIO}_build.log" 2>&1; then
        log "${RED}Build failed for $SCENARIO${NC}"
        RESULTS[$SCENARIO]="BUILD_FAIL"
        ((FAILED++))
        continue
    fi
    log "${GREEN}Build successful${NC}"

    # Run container
    log "Running container..."
    DOCKER_CMD="docker run --rm -e CASSIA_SCENARIO=$SCENARIO"

    if [ "$WITH_API_KEY" = true ] && [ -n "$CASSIA_API_KEY" ]; then
        DOCKER_CMD="$DOCKER_CMD -e CASSIA_API_KEY=$CASSIA_API_KEY"
        DOCKER_CMD="$DOCKER_CMD -e CASSIA_PROVIDER=${CASSIA_PROVIDER:-openrouter}"
        DOCKER_CMD="$DOCKER_CMD -e CASSIA_MODEL=${CASSIA_MODEL:-google/gemini-2.0-flash-001}"
    fi

    DOCKER_CMD="$DOCKER_CMD $TAG"

    START_TIME=$(date +%s)
    if $DOCKER_CMD > "$RUN_RESULTS_DIR/logs/${SCENARIO}.log" 2>&1; then
        EXIT_CODE=0
    else
        EXIT_CODE=$?
    fi
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))

    # Check result
    if [ $EXIT_CODE -eq 0 ]; then
        log "${GREEN}PASSED${NC} (${DURATION}s)"
        RESULTS[$SCENARIO]="PASS"
        ((PASSED++))
    else
        # Check if failure was expected
        if [ "$SCENARIO" = "minimal" ] || [ "$SCENARIO" = "old_python" ]; then
            log "${YELLOW}FAILED (expected)${NC} (${DURATION}s)"
            RESULTS[$SCENARIO]="FAIL_EXPECTED"
        else
            log "${RED}FAILED (unexpected)${NC} (${DURATION}s)"
            RESULTS[$SCENARIO]="FAIL"
            ((FAILED++))
        fi
    fi
done

# Summary
log ""
log "=============================================="
log "SUMMARY"
log "=============================================="
log "Passed: ${GREEN}$PASSED${NC}"
log "Failed: ${RED}$FAILED${NC}"
log ""

for SCENARIO in "${!RESULTS[@]}"; do
    RESULT="${RESULTS[$SCENARIO]}"
    case $RESULT in
        PASS)
            log "  $SCENARIO: ${GREEN}PASS${NC}"
            ;;
        FAIL)
            log "  $SCENARIO: ${RED}FAIL${NC}"
            ;;
        FAIL_EXPECTED)
            log "  $SCENARIO: ${YELLOW}FAIL (expected)${NC}"
            ;;
        *)
            log "  $SCENARIO: ${RED}$RESULT${NC}"
            ;;
    esac
done

log ""
log "Results saved to: $RUN_RESULTS_DIR"

# Exit with appropriate code
if [ $FAILED -gt 0 ]; then
    exit 1
fi
exit 0
