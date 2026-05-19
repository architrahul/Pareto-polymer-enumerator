"""
Paths for the Hilbert basis pipeline.

All paths are derived from this file's location so the repository is
portable: clone it anywhere and everything still resolves.

Normaliz is NOT redistributed with this repo. Set $NORMALIZ_EXE to a
locally installed `normaliz` binary, or place a built copy at
src/Normaliz/source/normaliz (which is gitignored). See README for install
instructions.
"""

import os
import shutil
import sys

REPO_ROOT        = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
EXAMPLE_TBNS_DIR = os.path.join(REPO_ROOT, "example-tbns")
RESULTS_DIR      = os.path.join(REPO_ROOT, "results")
LOGS_DIR         = os.path.join(RESULTS_DIR, "logs")
BENCHMARK_LOGS   = os.path.join(LOGS_DIR, "benchmarking")
LEAKAGE_DIR      = os.path.join(RESULTS_DIR, "leakage")
LEAKAGE_HB_DIR   = os.path.join(LEAKAGE_DIR, "hilbert_basis")
LEAKAGE_COFFEE_DIR = os.path.join(LEAKAGE_DIR, "coffee")
LEAKAGE_ANALYSIS_DIR = os.path.join(LEAKAGE_DIR, "analysis")
LEAKAGE_VARY_T_DIR = os.path.join(LEAKAGE_ANALYSIS_DIR, "vary_t")
LEAKAGE_VARY_INPUTS_DIR = os.path.join(LEAKAGE_ANALYSIS_DIR, "vary_removed_inputs")

DEFAULT_MONOMER_FILE = os.path.join(EXAMPLE_TBNS_DIR, "monomers_cascade_n8.txt")


def _resolve_normaliz_exe() -> str:
    """Return path to a `normaliz` executable, or print a clear error and exit.

    Resolution order:
        1. $NORMALIZ_EXE environment variable (if set and executable)
        2. src/Normaliz/source/normaliz (gitignored vendor location)
        3. `normaliz` on $PATH
    """
    candidate = os.environ.get("NORMALIZ_EXE")
    if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate

    vendored = os.path.join(REPO_ROOT, "src", "Normaliz", "source", "normaliz")
    if os.path.isfile(vendored) and os.access(vendored, os.X_OK):
        return vendored

    on_path = shutil.which("normaliz")
    if on_path:
        return on_path

    sys.stderr.write(
        "ERROR: Normaliz executable not found.\n"
        "  Install Normaliz locally and either:\n"
        "    • set $NORMALIZ_EXE to its absolute path, or\n"
        "    • place the built binary at src/Normaliz/source/normaliz, or\n"
        "    • make sure `normaliz` is on $PATH.\n"
        "  See README.md → 'Installing Normaliz' for instructions.\n"
    )
    sys.exit(1)


NORMALIZ_EXE = _resolve_normaliz_exe()
