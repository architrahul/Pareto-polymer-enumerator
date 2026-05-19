#!/bin/bash
# Master benchmark runner — paper Figures 2/3/4 + leakage analysis (5/6).
# Usage (from repo root):
#   caffeinate -i -s ./run_all_benchmarks.sh </dev/null
# All steps log to results/logs/benchmark_<timestamp>/.

set -u

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"
source .venv/bin/activate
cd src/benchmarks

TS="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$REPO_ROOT/results/logs/benchmark_$TS"
mkdir -p "$LOG_DIR"

MASTER="$LOG_DIR/00_master.log"
echo "Benchmark sweep started: $(date)" | tee "$MASTER"
echo "Logs: $LOG_DIR"               | tee -a "$MASTER"

run_step() {
    local name="$1"; shift
    local logf="$LOG_DIR/${name}.log"
    {
        echo
        echo "================================================================"
        echo "[$(date)] STEP: $name"
        echo "  cmd: $*"
        echo "================================================================"
    } | tee -a "$MASTER"

    "$@" </dev/null 2>&1 | tee "$logf"
    local rc=${PIPESTATUS[0]}

    echo "[$(date)] FINISHED $name (exit=$rc)" | tee -a "$MASTER"
    return $rc
}

# ----------------------------------------------------------------------------
# 1   Figure 2 — probe-estimated runtime vs k for damien, cascade n=7, n=8.
# ----------------------------------------------------------------------------
run_step "01_figure2_probe_curves" python experiments.py --phases 1 || true

# ----------------------------------------------------------------------------
# 2   Figure 3 + Figure 4 covering bars — cascade m=5..9, binary d=3,4,
#     dna m=4..7, t in {3..8}. Probe + full enumeration on best k.
# ----------------------------------------------------------------------------
run_step "02_figure34_covering" python experiments.py --phases 2 || true

# ----------------------------------------------------------------------------
# 3   Figure 3 + Figure 4 Full HB bars — single Normaliz run per (system,
#     size). NORMALIZ_TIMEOUT_SECONDS=3h enforces the "hrs+" truncation bar.
# ----------------------------------------------------------------------------
run_step "03_figure34_full_hb" python experiments.py --phases 3 || true

# ----------------------------------------------------------------------------
# 4   Render Figures 2/3/4 from the JSONs.
# ----------------------------------------------------------------------------
run_step "04_make_plots_234" python make_plots.py || true

# ----------------------------------------------------------------------------
# 5   Figures 5/6 — leakage equilibrium-error pipeline (cascade n=7 incomplete):
#     compute the three Hilbert bases → COFFEE → relative-error plots.
# ----------------------------------------------------------------------------
run_step "05_leakage_compute" python leakage_compute_all.py --n 7 || true
run_step "06_leakage_coffee"  python leakage_coffee.py      --n 7 || true
run_step "07_leakage_plots"   python plot_leakage_figures.py --n 7 || true

# ----------------------------------------------------------------------------
# 8   Sorted polymer dumps (CSV + top-50 text) for manual inspection.
# ----------------------------------------------------------------------------
run_step "08_dump_polymers"   python dump_sorted_polymers.py --n 7 --top 50 || true

echo                                                | tee -a "$MASTER"
echo "================================================================" | tee -a "$MASTER"
echo "Benchmark sweep complete: $(date)"                | tee -a "$MASTER"
echo "================================================================" | tee -a "$MASTER"
