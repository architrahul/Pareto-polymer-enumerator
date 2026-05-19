#!/bin/bash
# Rebuild every downstream artefact that depends on Full-HB outputs after
# Phase 3 has been rerun.
#
# Usage, from the repo root:
#   caffeinate -i -s ./rerun_after_phase3.sh </dev/null
#
# This script intentionally does NOT rerun Phase 3 itself. It first checks
# that the fresh Phase-3 Hilbert-basis files have been projected back to the
# monomer coordinate space, then archives stale downstream leakage outputs,
# and finally regenerates all later figures / leakage analyses / CSV exports.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

if [[ -f ".venv/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source .venv/bin/activate
fi

TS="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$REPO_ROOT/results/logs/post_phase3_$TS"
ARCHIVE_DIR="$REPO_ROOT/backup/results/post_phase3_$TS"
mkdir -p "$LOG_DIR" "$ARCHIVE_DIR"

MASTER="$LOG_DIR/00_master.log"
echo "Post-Phase-3 rebuild started: $(date)" | tee "$MASTER"
echo "Logs:    $LOG_DIR"                     | tee -a "$MASTER"
echo "Archive: $ARCHIVE_DIR"                 | tee -a "$MASTER"

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
    return "$rc"
}

archive_if_present() {
    local rel="$1"
    local src="$REPO_ROOT/$rel"
    local dst="$ARCHIVE_DIR/$rel"
    if [[ -e "$src" ]]; then
        mkdir -p "$(dirname "$dst")"
        mv "$src" "$dst"
        echo "archived $rel" | tee -a "$MASTER"
    fi
}

validate_phase3() {
    python - "$REPO_ROOT" <<'PY'
import json
import os
import sys

repo = sys.argv[1]
path = os.path.join(repo, "results", "experiments", "phase3_full_hb.json")
if not os.path.exists(path):
    raise SystemExit(f"missing {path}; finish Phase 3 before running this script")

payload = json.load(open(path))
bad = []
checked = 0
for rec in payload.get("results", []):
    hb_path = rec.get("hb_path")
    if rec.get("truncated") or not hb_path:
        continue
    if not os.path.exists(hb_path):
        bad.append(f"missing HB file: {hb_path}")
        continue
    expected = rec["n_monomers"]
    width = None
    with open(hb_path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                width = len(line.split())
                break
    if width is None:
        bad.append(f"empty HB file: {hb_path}")
    elif width != expected:
        bad.append(
            f"{os.path.basename(hb_path)} has vector width {width}, expected {expected}"
        )
    checked += 1

if bad:
    print("Phase-3 validation failed:")
    for item in bad:
        print(f"  - {item}")
    raise SystemExit(1)

print(f"Phase-3 validation passed for {checked} completed Full-HB file(s).")
PY
}

validate_downstream_outputs() {
    python - "$REPO_ROOT" <<'PY'
import os
import sys

repo = sys.argv[1]
required = [
    # Paper leakage chain
    "results/leakage/hilbert_basis/hilbert_k25_t3_monomer_n7_incomplete.txt",
    "results/leakage/hilbert_basis/hilbert_k25_t5_monomer_n7_incomplete.txt",
    "results/leakage/hilbert_basis/hilbert_full_p_star_n7_incomplete.txt",
    "results/leakage/coffee/n7_incomplete/full_pstar/coffee_output.txt",
    "results/leakage/coffee/n7_incomplete/full_pstar/coffee_output_sorted.csv",
    "results/leakage/coffee/n7_incomplete/k25_t3/coffee_output.txt",
    "results/leakage/coffee/n7_incomplete/k25_t3/coffee_output_sorted.csv",
    "results/leakage/coffee/n7_incomplete/k25_t5/coffee_output.txt",
    "results/leakage/coffee/n7_incomplete/k25_t5/coffee_output_sorted.csv",
    "results/figures/figure5_leakage_t3.png",
    "results/figures/figure6_leakage_t5.png",

    # Extended leakage-vs-t experiment
    "results/leakage/hilbert_basis/hilbert_k25_t4_monomer_n7_incomplete.txt",
    "results/leakage/hilbert_basis/hilbert_k25_t6_monomer_n7_incomplete.txt",
    "results/leakage/hilbert_basis/hilbert_k25_t7_monomer_n7_incomplete.txt",
    "results/leakage/analysis/vary_t/n7_incomplete/summary.json",
    "results/leakage/analysis/vary_t/n7_incomplete/figure_leakage_vs_t.png",

    # Extended leakage-vs-K experiment
    "results/leakage/hilbert_basis/hilbert_full_p_star_n7_incomplete2.txt",
    "results/leakage/hilbert_basis/hilbert_full_p_star_n7_incomplete3.txt",
    "results/leakage/hilbert_basis/hilbert_full_p_star_n7_incomplete4.txt",
    "results/leakage/hilbert_basis/hilbert_full_p_star_n7_incomplete5.txt",
    "results/leakage/hilbert_basis/hilbert_full_p_star_n7_incomplete6.txt",
    "results/leakage/hilbert_basis/hilbert_full_p_star_n7_incomplete7.txt",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/summary.json",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/figure_leakage_vs_K.png",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/K1/polymers_sorted.csv",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/K2/polymers_sorted.csv",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/K3/polymers_sorted.csv",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/K4/polymers_sorted.csv",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/K5/polymers_sorted.csv",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/K6/polymers_sorted.csv",
    "results/leakage/analysis/vary_removed_inputs/n7_systems_compare/K7/polymers_sorted.csv",

    # External plotting exports
    "results/csv/leakage_vs_t.csv",
    "results/csv/leakage_vs_t_per_polymer.csv",
    "results/csv/leakage_vs_K.csv",
    "results/csv/leakage_vs_K_per_polymer.csv",
]

missing = [rel for rel in required if not os.path.exists(os.path.join(repo, rel))]
if missing:
    print("Downstream rebuild validation failed; missing:")
    for rel in missing:
        print(f"  - {rel}")
    raise SystemExit(1)
print(f"Downstream rebuild validation passed for {len(required)} expected artefacts.")
PY
}

echo | tee -a "$MASTER"
echo "Validating fresh Phase-3 Full-HB outputs before downstream rebuild..." | tee -a "$MASTER"
validate_phase3 | tee "$LOG_DIR/01_validate_phase3.log"

# Archive only artefacts that depend on stale Full-HB files. Covering-design
# Hilbert bases are deliberately kept: their coordinates were already in the
# correct monomer space, and reusing them avoids needless long Normaliz runs.
echo | tee -a "$MASTER"
echo "Archiving stale downstream Full-HB-derived outputs..." | tee -a "$MASTER"
archive_if_present "results/leakage"
archive_if_present "results/figures/figure5_leakage_t3.png"
archive_if_present "results/figures/figure6_leakage_t5.png"
archive_if_present "results/csv/leakage_vs_t.csv"
archive_if_present "results/csv/leakage_vs_t_per_polymer.csv"
archive_if_present "results/csv/leakage_vs_K.csv"
archive_if_present "results/csv/leakage_vs_K_per_polymer.csv"

cd "$REPO_ROOT/src/benchmarks"

# Figures 2/3/4 now read the fresh Phase-3 JSON/HB outputs.
run_step "02_make_plots_234" python -u make_plots.py

# Paper leakage figures (Figures 5/6).
run_step "03_leakage_compute" python -u leakage_compute_all.py --n 7
run_step "04_leakage_coffee"  python -u leakage_coffee.py      --n 7
run_step "05_leakage_plots"   python -u plot_leakage_figures.py --n 7
run_step "06_dump_polymers"   python -u dump_sorted_polymers.py --n 7 --top 50

# Extended leakage experiments added after the paper figures.
run_step "07_leakage_vs_t"     python -u leakage_experiment_t.py      --n 7
run_step "08_leakage_vs_K"     python -u leakage_experiment_inputs.py --n 7 --K-values 1 2 3 4 5 6 7

# Final tidy exports for external plotting.
run_step "09_export_csv"       python -u export_csv.py
echo | tee -a "$MASTER"
echo "Validating rebuilt downstream outputs..." | tee -a "$MASTER"
validate_downstream_outputs | tee "$LOG_DIR/10_validate_outputs.log"

echo                                                | tee -a "$MASTER"
echo "================================================================" | tee -a "$MASTER"
echo "Post-Phase-3 rebuild complete: $(date)"           | tee -a "$MASTER"
echo "================================================================" | tee -a "$MASTER"
