#!/bin/bash
# Benchmark runner.
#
# Usage from repo root:
#   ./run_all_benchmarks.sh          # no verbose log files
#   ./run_all_benchmarks.sh --logs   # write verbose logs too
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

if [[ -f .venv/bin/activate ]]; then
  source .venv/bin/activate
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  PYTHON_BIN="python"
fi

WRITE_LOGS=0
for arg in "$@"; do
  if [[ "$arg" == "--logs" ]]; then
    WRITE_LOGS=1
    break
  fi
done

CMD=("$PYTHON_BIN" -u src/benchmarks/run_benchmarks.py --experiments all "$@")
if command -v caffeinate >/dev/null 2>&1; then
  CMD=(caffeinate -i -s "${CMD[@]}")
fi

if [[ "$WRITE_LOGS" == "1" ]]; then
  mkdir -p results/benchmarks/logs
  LOG="results/benchmarks/logs/run_all.log"
  echo "[$(date)] Starting full benchmark run" | tee "$LOG"
  echo "Repo: $REPO_ROOT" | tee -a "$LOG"
  echo "Command: ${CMD[*]}" | tee -a "$LOG"
  echo "Log: $LOG" | tee -a "$LOG"
  echo | tee -a "$LOG"

  "${CMD[@]}" 2>&1 | tee -a "$LOG"
  status=${PIPESTATUS[0]}

  echo | tee -a "$LOG"
  echo "[$(date)] Finished with exit code $status" | tee -a "$LOG"
  exit "$status"
else
  echo "Running without verbose log files. Use --logs to write results/benchmarks/logs/run_all.log."
  exec "${CMD[@]}"
fi
