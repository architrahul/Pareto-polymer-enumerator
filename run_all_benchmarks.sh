#!/bin/bash
# Consolidated benchmark runner. Usage from repo root:
#   caffeinate -i -s ./run_all_benchmarks.sh </dev/null
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"
if [[ -f .venv/bin/activate ]]; then
  source .venv/bin/activate
fi
python -u src/benchmarks/run_benchmarks.py --experiments all
